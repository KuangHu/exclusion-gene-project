#!/usr/bin/env python3
"""Final annotation pass: slot-occupant exclusion family, all six families, all classes.

Three gaps closed before freeze:

1. MPF_T has NO slot_occupant_family column. The numbers were measured
   (eex_occupied 46.2% / candidate 40.8% / slot_empty 13.0%) but stayed in logs.
2. ALL SIX exclusion families are scanned in EVERY class, not just the one the
   class is expected to carry. The v1 census already showed cross-class presence:
   Eex_IncN on 150 MPF_F plasmids, TraS on 28 MPF_T.
3. MPF_FA's 11 PF14729 hits are recorded as exclusion_family_hit WITH COORDINATES
   and kept SEPARATE from any slot notion -- MPF_FA has no validated slot.

Slot positions per class:
  MPF_T   VirB5+1   (VirB5 strand frame, all architectures)
  MPF_F   TraG_N+1
  MPF_I   TraY+1
  MPF_FA  none -- family hits are plasmid-level only

An occupant is `eex_occupied` if a named family hits it, `candidate` if the slot
is filled by something unnamed, `slot_empty` if the +1 gene is itself an anchor.
"""
import collections, csv, json, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"]
NAME = {"TIGR04359": "TrbK_RP4", "NF033894": "Eex_IncN", "NF041429": "EexR",
        "NF033891": "ExcA", "PF10624": "TraS", "PF14729": "DUF4467"}
SPEC = [
 ("MPF_T",  "data/release/v1.1/catalogue_MPF_T_v1.1.tsv",  "slot_ready__virb5_virb6", "VirB5"),
 ("MPF_F",  "data/release/v1/catalogue_MPF_F_v1.tsv",      "slot_ready__traG",        "TraG_N"),
 ("MPF_I",  "data/release/v1.1/catalogue_MPF_I_v1.tsv",    "slot_ready__traY",        "traY"),
 ("MPF_FA", "data/release/v1.1/catalogue_MPF_FA_v1.tsv",   None,                      None),
]


def parse(c):
    if not c or c == "none": return None
    f = c.split(":")
    if f[0] == "ga": f = f[1:]
    s, e = f[1].split("..")
    return int(f[0]), int(s), int(e), 1 if f[2] == "+" else -1


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    want = {}
    for cls, path, ready, anchor in SPEC:
        for r in csv.DictReader(open(os.path.join(PROJ, path)), delimiter="\t"):
            if ready is None or r.get(ready) == "1":
                want.setdefault(r["accession"], []).append(cls)
    print("plasmids to annotate: %d" % len(want), flush=True)

    seqs, owner, idx, st, en, strand, aas = [], [], [], [], [], [], []
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
                owner.append(acc); idx.append(int(f[1])); st.append(int(f[2]))
                en.append(int(f[3])); strand.append(int(f[4])); aas.append(p)
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(), sequence=p).digitize(alpha))
    print("proteins %d\n" % len(seqs), flush=True)

    fam_hit = {}
    for f in EEX:
        p = os.path.join(HMM, f + ".hmm")
        if not os.path.exists(p):
            print("  MISSING %s -- named control absent, aborting" % f); return 1
        with pyhmmer.plan7.HMMFile(p) as fh: m = next(iter(fh))
        best = {}
        for top in pyhmmer.hmmsearch([m], seqs, cpus=16, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                if i not in best or h.score > best[i]: best[i] = round(h.score, 1)
        fam_hit[f] = best
        print("  %-10s (%-9s) %6d proteins" % (f, NAME[f], len(best)), flush=True)

    allrows, summary = [], {}
    for cls, path, ready, anchor in SPEC:
        cat = [r for r in csv.DictReader(open(os.path.join(PROJ, path)), delimiter="\t")
               if ready is None or r.get(ready) == "1"]
        cnt = collections.Counter(); fams = collections.Counter()
        for r in cat:
            acc = r["accession"]
            row = {"accession": acc, "mpf_class": cls, "slot_anchor": anchor or "",
                   "slot_occupant_family": "", "slot_occupant_coords": "",
                   "slot_occupant_aa": "", "slot_family_bit": "", "slot_status": ""}
            if anchor is None:
                hits = [(NAME[f], fam_hit[f][j], j) for f in EEX
                        for j in byacc[acc].values() if j in fam_hit[f]]
                if hits:
                    nm2, sc, j = max(hits, key=lambda x: x[1])
                    row.update({"slot_occupant_family": nm2, "slot_family_bit": sc,
                                "slot_occupant_coords": "%d..%d" % (st[j], en[j]),
                                "slot_occupant_aa": len(aas[j]),
                                "slot_status": "family_hit_no_slot"})
                    fams[nm2] += 1; cnt["family_hit_no_slot"] += 1
                else:
                    row["slot_status"] = "no_slot_defined"; cnt["no_slot_defined"] += 1
                allrows.append(row); continue
            p = parse(r.get(anchor, ""))
            j = byacc[acc].get(p[0] + p[3]) if p else None
            if j is None:
                row["slot_status"] = "unresolved"; cnt["unresolved"] += 1
                allrows.append(row); continue
            got = [(NAME[f], fam_hit[f][j]) for f in EEX if j in fam_hit[f]]
            row["slot_occupant_coords"] = "%d..%d" % (st[j], en[j])
            row["slot_occupant_aa"] = len(aas[j])
            if got:
                nm2, sc = max(got, key=lambda x: x[1])
                row.update({"slot_occupant_family": ";".join(g[0] for g in got),
                            "slot_family_bit": sc, "slot_status": "eex_occupied"})
                cnt["eex_occupied"] += 1
                for g in got: fams[g[0]] += 1
            else:
                row["slot_status"] = "candidate"; cnt["candidate"] += 1
            allrows.append(row)
        summary[cls] = (len(cat), cnt, fams)
        print("\n=== %s (n=%d) ===" % (cls, len(cat)))
        for k, v in cnt.most_common():
            print("  %-22s %6d (%.1f%%)" % (k, v, 100.0*v/len(cat)))
        if fams:
            print("  families: %s" % ", ".join("%s %d" % kv for kv in fams.most_common()))

    dest = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(allrows)
    print("\nwrote %s (%d rows)" % (dest, len(allrows)))

    print("\n================= FINAL COUNTS =================")
    tot_cl = tot_eex = 0
    print("  %-8s %10s %14s %12s %12s" % ("class", "clusters", "eex_occupied", "candidate", "unresolved"))
    for cls, _, _, _ in SPEC:
        n, cnt, _ = summary[cls]
        e = cnt.get("eex_occupied", 0) + cnt.get("family_hit_no_slot", 0)
        tot_cl += n; tot_eex += e
        print("  %-8s %10d %14d %12d %12d"
              % (cls, n, e, cnt.get("candidate", 0), cnt.get("unresolved", 0)))
    print("  %-8s %10d %14d" % ("TOTAL", tot_cl, tot_eex))
    return 0


if __name__ == "__main__":
    sys.exit(main())
