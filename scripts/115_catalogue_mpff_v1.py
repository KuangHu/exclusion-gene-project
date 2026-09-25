#!/usr/bin/env python3
"""MPF_F tra gene cluster catalogue, v1. Same shape as MPF_T v1, four divergences.

NOT an operon catalogue: HMM hits and gene positions only.

Entry criterion: TraC_F_IV / PF11130, measured across all 72,558 cached
accessions -- admits 11,120 (15.4%), all five seeds at GA, 93% of admissions
independently called MPF_F by MOB-suite. Chosen as the VirB4 homolog's F-SPECIFIC
domain, not as the narrowest filter; PF19044 is the same protein's generic P-loop
domain and admits 18,223 at 57% purity.

Divergences from MPF_T v1, each because the class differs rather than by choice:

  core_completeness   over the 9 F-SPECIFIC families (>90% MPF_F purity). The four
                      cross-class families (TrbI, P-loop_TraG, TrwB_AAD_bind,
                      TraG-D_C) are genuine MPF_F components but admit the MPF_T
                      set wholesale, so they cannot score class membership.
  layout + order      two columns, not one. Layout is four-valued (contiguous
                      64.8%, split_2 22.1%, split_3 11.8%, split_4+ 1.3%) and a
                      gene-order string cannot express R27's Tra1/Tra2 separation.
  slot_ready__traG    the slot is TraG_N+1, validated on the seeds: F traS,
                      R100 traS, SXT EexR1 are all named exclusion genes at +1.
  traN_present        replaces virB1/virB7 columns; MPF_F has no analogue of
                      either. TraC and TraD are MULTI-DOMAIN, not fused -- two
                      models of one protein, unlike VirB3|VirB4.

Caveat carried into the file: MPF_T's slot is 5' of VirB6 and MPF_F's is 3' of the
VirB6 analogue -- OPPOSITE sides. MPF_F also has no VirB5 analogue in its anchor
set (the cross-class map resolved only VirB4, VirD4, VirB10). So "which side of
VirB6 the slot is on" is not a general rule. Recorded, not explained.
"""
import argparse
import collections
import csv
import os
import re
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
PLSDB = "/global/scratch/users/kh36969/plsdb"
ENTRY = "TraC_F_IV"
FSPEC = ["TraU", "TraH", "TraC_F_IV", "TraG_N", "TraE", "TrbC_Ftype",
         "F_T4SS_TraN", "TraV", "TraF"]
CROSSCLASS = ["TrbI", "P-loop_TraG", "TrwB_AAD_bind", "TraG-D_C"]
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "TraS", "DUF4467"]
LIPO = re.compile(r"[LVI][ASTVIG][GASN]C")
SPLIT_GAP = 20000


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    lines = [l for l in open(os.path.join(PROJ, "data", "anchors",
                                          "anchor_set_MPF_F.tsv"))
             if not l.startswith("#") and l.strip()]
    fams = [(r["pfam_acc"], r["pfam_name"]) for r in csv.DictReader(lines, delimiter="\t")]

    seqs, owner, idx, st, en, strand, aas = [], [], [], [], [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                owner.append(f[0]); idx.append(int(f[1]))
                st.append(int(f[2])); en.append(int(f[3])); strand.append(int(f[4]))
                aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    print("cache: %d proteins over %d accessions" % (len(seqs), len(set(owner))), flush=True)
    bypos = {(owner[i], idx[i]): i for i in range(len(owner))}
    byacc = collections.defaultdict(list)
    for i in range(len(owner)):
        byacc[owner[i]].append(i)

    def scan(fam):
        p = None
        for c in ("mpff_%s.hmm" % fam, "anchor_%s.hmm" % fam, "%s.hmm" % fam):
            if os.path.exists(os.path.join(HMM, c)):
                p = os.path.join(HMM, c); break
        if not p:
            return None
        with pyhmmer.plan7.HMMFile(p) as fh:
            model = next(iter(fh))
        best = {}
        for top in pyhmmer.hmmsearch([model], seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                if owner[i] not in best or h.score > best[owner[i]][1]:
                    best[owner[i]] = (i, h.score)
        return best
    hit = {}
    for acc, name in fams:
        b = scan(name)
        if b is not None:
            hit[name] = b
            print("  %-18s %6d plasmids" % (name, len(b)), flush=True)
    eexhit = {}
    for fam in EEX:
        b = scan(fam)
        if b:
            for accn, (i, sc) in b.items():
                eexhit.setdefault(i, fam)

    mpf, rep = {}, collections.defaultdict(set)
    csv.field_size_limit(10 ** 7)
    for r in csv.DictReader(open(os.path.join(PLSDB, "typing.csv"))):
        mpf[r["NUCCORE_ACC"]] = r.get("mpf_type", "") or "-"
    for r in csv.DictReader(open(os.path.join(PLSDB, "plasmidfinder.csv"))):
        rep[r["NUCCORE_ACC"]].add(r["typing"])

    admitted = sorted(hit.get(ENTRY, {}))
    print("\nadmitted by %s: %d plasmids" % (ENTRY, len(admitted)))
    rows = []
    for acc in admitted:
        pos = {n: hit[n][acc][0] for n in hit if acc in hit[n]}
        nf = sum(1 for n in FSPEC if n in pos)
        ncc = sum(1 for n in CROSSCLASS if n in pos)
        p = sorted(st[i] for i in pos.values())
        gaps = [b - a2 for a2, b in zip(p, p[1:])]
        big = [g for g in gaps if g >= SPLIT_GAP]
        layout = "contiguous" if not big else "split_%d" % (len(big) + 1)
        tier = "A" if nf == len(FSPEC) else "B" if nf >= len(FSPEC) - 1 else \
               "C" if nf >= len(FSPEC) - 3 else "D"
        # the slot: gene immediately 3' of TraG_N
        slot = slotlen = slotfam = slotlipo = ""
        ready = 0
        if "TraG_N" in pos:
            g = pos["TraG_N"]; t = 1 if strand[g] == 1 else -1
            j = bypos.get((acc, idx[g] + t))
            if j is not None and strand[j] == strand[g]:
                ready = 1
                slot = "%d..%d%s" % (st[j], en[j], "+" if strand[j] == 1 else "-")
                slotlen = len(aas[j])
                slotfam = eexhit.get(j, "")
                slotlipo = int(bool(LIPO.search(aas[j][:40])))
        r = {"accession": acc, "mob_mpf_type": mpf.get(acc, "-"),
             "replicons": ";".join(sorted(rep.get(acc, ()))) or "-",
             "n_cds": len(byacc[acc]), "tier": tier,
             "core_completeness_fspecific": "%d/%d" % (nf, len(FSPEC)),
             "n_fspecific": nf, "n_crossclass": ncc,
             "layout": layout, "anchor_span_bp": p[-1] - p[0] if len(p) > 1 else 0,
             "max_internal_gap_bp": max(gaps) if gaps else 0,
             "traN_present": int("F_T4SS_TraN" in pos),
             "slot_ready__traG": ready, "slot_coords": slot,
             "slot_aa_len": slotlen, "slot_eex_family": slotfam,
             "slot_lipobox": slotlipo}
        for n, _ in fams:
            pass
        for acc2, name in fams:
            i = pos.get(name)
            r[name] = ("ga:%d:%d..%d:%s" % (idx[i], st[i], en[i],
                       "+" if strand[i] == 1 else "-")) if i is not None else "none"
        rows.append(r)

    cols = list(rows[0].keys())
    dest = os.path.join(PROJ, "data", "catalogue_MPF_F_regen.tsv")   # A9: see 106
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("wrote %s (%d rows, %d cols)" % (dest, len(rows), len(cols)))

    print("\n=== MPF_F CATALOGUE v1 SUMMARY ===")
    for k in ("tier", "layout"):
        c = collections.Counter(r[k] for r in rows)
        print("  %-8s %s" % (k, ", ".join("%s %d (%.1f%%)" % (a2, b, 100.0 * b / len(rows))
                                          for a2, b in c.most_common())))
    print("  slot_ready__traG : %d (%.1f%%)"
          % (sum(r["slot_ready__traG"] for r in rows),
             100.0 * sum(r["slot_ready__traG"] for r in rows) / len(rows)))
    sl = [r for r in rows if r["slot_ready__traG"]]
    named = [r for r in sl if r["slot_eex_family"]]
    print("  slot occupant hit by a named Eex family: %d/%d (%.1f%%)"
          % (len(named), len(sl), 100.0 * len(named) / max(1, len(sl))))
    fc = collections.Counter(r["slot_eex_family"] for r in named)
    for k, v in fc.most_common():
        print("      %-12s %d" % (k, v))
    L = sorted(r["slot_aa_len"] for r in sl if r["slot_aa_len"])
    if L:
        print("  slot occupant length: median %d  Q1 %d  Q3 %d" % (L[len(L)//2], L[len(L)//4], L[3*len(L)//4]))
    lb = [r["slot_lipobox"] for r in sl if r["slot_lipobox"] != ""]
    print("  slot occupant lipobox-positive: %.1f%%" % (100.0 * sum(lb) / max(1, len(lb))))
    print("  unnamed + lipobox-positive (the candidate set): %d"
          % sum(1 for r in sl if not r["slot_eex_family"] and r["slot_lipobox"] == 1))


if __name__ == "__main__":
    main()
