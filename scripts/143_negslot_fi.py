#!/usr/bin/env python3
"""Negative slots for MPF_F and MPF_I -- the controls these classes never had.

MPF_T's slot work rests on negative slots (VirB8|B9, VirB9|B10, VirB10|B11).
MPF_F and MPF_I have never had them, so for those classes it is not currently
possible to test whether ANY criterion separates. This builds them.

ANCHOR PAIRS, chosen from MEASURED backbone adjacency, not assumed order:

  MPF_F   TrbC_Ftype -> F_T4SS_TraN  74.6%    TraU -> TrbC_Ftype  54.2%
          (only two usable: the sole other pair >=30% is TraH -> TraG_N at 96.1%,
           which is the target slot's own side)
  MPF_I   traM -> traN 97.1%   traP -> traQ 93.2%   trbA -> trbB 88.7%

`TraG_N` and `traY` appear in NO high-adjacency pair -- their +1 is not normally
another anchor, which is what makes them slots.

The target slot is run through this SAME code for comparison, so any difference is
not an artefact of two implementations.

DECONTAMINATION uses each class's OWN anchor families -- MPF_F's 13, MPF_I's 17 --
not MPF_T's 14. In a colinear operon anchor+1 IS the next anchor, so a positional
nomination set is contaminated by construction.

Reported per slot: pair occupancy, occupant length distribution, named-family
hits, lipobox rate. Lipobox is measured even where the class is not expected to
use it, because it is the one axis comparable ACROSS classes.
"""
import argparse, collections, csv, json, os, re, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
OUT = "/global/scratch/users/kh36969/exclusion_gene/negslot_fi"
LIPO = re.compile(r"[LVI][ASTVIG][GASN]C")
EEX = ["PF10624", "NF033891", "NF033894", "NF041429", "TIGR04359", "PF14729"]

CLASSES = {
  "MPF_F": {"cat": "data/release/v1/catalogue_MPF_F_v1.tsv", "ready": "slot_ready__traG",
            "fams": ["TraU","TraH","TraF","TraC_F_IV","TraG_N","TrbI","P-loop_TraG","TraE",
                     "TrbC_Ftype","F_T4SS_TraN","TrwB_AAD_bind","TraV","TraG-D_C"],
            "hmm": lambda f: [os.path.join(HMM, "mpff_%s.hmm" % f),
                              os.path.join(HMM, "anchor_%s.hmm" % f),
                              os.path.join(HMM, "%s.hmm" % f)],
            "target": "TraG_N"},
  "MPF_I": {"cat": "data/release/v1.1/catalogue_MPF_I_v1.tsv", "ready": "slot_ready__traY",
            "fams": ["traE","traI","traK","traL","traM","traN","traO","traP","traQ","traR",
                     "traT","traU","traV","traW","traY","trbA","trbB"],
            "hmm": lambda f: [os.path.join(CJ, "T4SS_I_%s.hmm" % f)],
            "target": "traY"},
}


def parse(c):
    if not c or c == "none": return None
    f = c.split(":")
    if f[0] == "ga": f = f[1:]
    s, e = f[1].split("..")
    return int(f[0]), int(s), int(e), 1 if f[2] == "+" else -1


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cls", required=True, choices=list(CLASSES))
    ap.add_argument("--anchor", required=True)
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    C = CLASSES[a.cls]
    os.makedirs(OUT, exist_ok=True)
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    cat = [r for r in csv.DictReader(open(os.path.join(PROJ, C["cat"])), delimiter="\t")
           if r.get(C["ready"]) == "1"]
    want = {r["accession"] for r in cat}
    is_target = (a.anchor == C["target"])
    print("%s  anchor %s%s   nomination set %d plasmids"
          % (a.cls, a.anchor, "  (TARGET SLOT)" if is_target else "  (negative)", len(cat)),
          flush=True)

    seqs, owner, aas = [], [], []
    byacc = collections.defaultdict(dict)
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                if acc not in want: continue
                p = line.rstrip("\n")
                byacc[acc][int(f[1])] = len(seqs)
                owner.append(acc); aas.append(p)
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(), sequence=p).digitize(alpha))
    print("  cache subset: %d proteins over %d plasmids" % (len(seqs), len(byacc)), flush=True)

    nom, nocds, noanchor = {}, 0, 0
    for r in cat:
        p = parse(r.get(a.anchor, ""))
        if not p: noanchor += 1; continue
        gi, _, _, sd = p
        j = byacc[r["accession"]].get(gi + sd)
        if j is None: nocds += 1; continue
        nom[r["accession"]] = aas[j]
    occ = 100.0 * len(nom) / len(cat)
    print("  PAIR OCCUPANCY: %d/%d = %.1f%%  (anchor absent %d, no CDS at +1 %d)"
          % (len(nom), len(cat), occ, noanchor, nocds), flush=True)
    uniq = sorted(set(nom.values()))
    print("  records %d | unique %d" % (len(nom), len(uniq)), flush=True)

    # --- decontaminate against this class's OWN anchor families -------------
    members = {}
    for f in C["fams"]:
        path = next((p for p in C["hmm"](f) if os.path.exists(p)), None)
        if not path: continue
        with pyhmmer.plan7.HMMFile(path) as fh: m = next(iter(fh))
        hits = set()
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            for h in top: hits.add(aas[int(dec(h.name))])
        if hits: members[f] = sorted(hits)
    print("  confirmed members: %s"
          % ", ".join("%s %d" % (k, len(v)) for k, v in sorted(members.items())), flush=True)

    dig = lambda L, t: [pyhmmer.easel.TextSequence(
        name=("%s%d" % (t, i)).encode(), sequence=s).digitize(alpha) for i, s in enumerate(L)]
    q = dig(uniq, "q")
    drop = collections.defaultdict(set)
    for fam, ref in members.items():
        bestE = {}
        for top in pyhmmer.phmmer(q, dig(ref, "r"), cpus=a.cpus, E=1000.0):
            i = int(dec(top.query.name)[1:])
            for h in top:
                if i not in bestE or h.evalue < bestE[i]: bestE[i] = h.evalue
        for i, e in bestE.items():
            if e <= 1e-3: drop[fam].add(uniq[i])
    dropped = set().union(*drop.values()) if drop else set()
    kept = [s for s in uniq if s not in dropped]
    print("  decontamination: %d -> %d kept (%d removed, %.1f%%)"
          % (len(uniq), len(kept), len(dropped), 100.0*len(dropped)/max(1, len(uniq))))
    for f, s in sorted(drop.items(), key=lambda x: -len(x[1]))[:6]:
        print("      %-14s %d" % (f, len(s)))
    keptset = set(kept)
    recs = {k: v for k, v in nom.items() if v in keptset}

    # --- characterise --------------------------------------------------------
    L = sorted(len(s) for s in kept)
    lip = sum(1 for s in kept if LIPO.search(s[:40]))
    print("\n  AFTER DECONTAMINATION: records %d | unique %d" % (len(recs), len(kept)))
    if L:
        print("  length: Q1 %d median %d Q3 %d (min %d max %d)"
              % (L[len(L)//4], L[len(L)//2], L[3*len(L)//4], L[0], L[-1]))
        print("  lipobox: %d/%d = %.1f%%" % (lip, len(kept), 100.0*lip/len(kept)))

    named = collections.Counter()
    if kept:
        ks = dig(kept, "k")
        for fam in EEX:
            p = os.path.join(HMM, fam + ".hmm")
            if not os.path.exists(p): continue
            with pyhmmer.plan7.HMMFile(p) as fh: m = next(iter(fh))
            for top in pyhmmer.hmmsearch([m], ks, cpus=a.cpus, bit_cutoffs="gathering"):
                for h in top: named[fam] += 1
    print("  named exclusion family hits on kept occupants: %s"
          % (", ".join("%s %d" % kv for kv in named.most_common()) or "NONE"))

    row = {"class": a.cls, "anchor": a.anchor, "is_target": int(is_target),
           "n_plasmids": len(cat), "pair_occupancy_pct": round(occ, 1),
           "records_nominated": len(nom), "unique_nominated": len(uniq),
           "unique_kept": len(kept), "records_kept": len(recs),
           "pct_removed": round(100.0*len(dropped)/max(1, len(uniq)), 1),
           "len_q1": L[len(L)//4] if L else "", "len_median": L[len(L)//2] if L else "",
           "len_q3": L[3*len(L)//4] if L else "",
           "lipobox_pct": round(100.0*lip/len(kept), 1) if kept else "",
           "named_family_hits": ";".join("%s:%d" % kv for kv in named.most_common()) or ""}
    dest = os.path.join(PROJ, "data", "anchors", "negslot_fi.tsv")
    ex = []
    if os.path.exists(dest):
        ex = [r for r in csv.DictReader(open(dest), delimiter="\t")
              if not (r["class"] == a.cls and r["anchor"] == a.anchor)]
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row.keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in ex + [row]: w.writerow(r)
    with open(os.path.join(OUT, "%s_%s.json" % (a.cls, a.anchor)), "w") as fh:
        json.dump({"records": recs}, fh)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
