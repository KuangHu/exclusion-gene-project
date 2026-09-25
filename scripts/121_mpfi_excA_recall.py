#!/usr/bin/env python3
"""Tiebreaker for the MPF_I entry criterion: recall on a known-positive set.

Ten T4SS_I_* profiles qualify at >=90% MPF_I purity and fall into two near-nested
groups -- narrow (traT/traV/traE, ~1,850 admitted, 94.8-99.0% pure) and broad
(traL/traP/traQ/trbB/traO/traK/traW, ~3,500, 90.1-94.0%). Cross-group Jaccard is
0.520 against |narrow|/|broad| = 0.532, so narrow is essentially a subset and the
broad profiles add ~1,650 plasmids at ~3 points of purity.

Purity alone cannot choose between them. But MPF_I has an independent positive set
that MPF_T and MPF_F never had: the Eex census found NF033891 (ExcA) on 2,293
plasmids, of which only 28 are admitted by the existing two entry criteria. Those
~2,265 carry a named exclusion family and belong to no built class.

Read the recall THREE ways, decided in advance:

  narrow captures ~all of them   -> narrow is genuinely specific; broad's extra
                                    1,650 are noise
  narrow captures about half     -> narrow is UNDER-DETECTING; broad's extra are
                                    real MPF_I
  both capture ~all of them      -> the extra 1,650 carry no ExcA and this set
                                    cannot separate them; fall back to describing
                                    what those 1,650 are

Second question, one extra column: I_traU misses the 90% threshold by 1.3 points.
What is its 11.3% impurity made of? MPF_F excluded four families that were eating
6,000-9,000 plasmids of another CLASS. If I_traU's impurity is mostly MOB-suite
'-' (untyped), 88.7% means something entirely different from cross-class
contamination -- it may be MPF_I that MOB-suite failed to type.
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
PLSDB = "/global/scratch/users/kh36969/plsdb"
NARROW = ["T4SS_I_traT", "T4SS_I_traV", "T4SS_I_traE"]
BROAD = ["T4SS_I_traL", "T4SS_I_traP", "T4SS_I_traQ", "T4SS_I_trbB",
         "T4SS_I_traO", "T4SS_I_traK", "T4SS_I_traW"]


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner = [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                owner.append(nm.split("|")[0])
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    N = len(set(owner))
    print("full cache: %d proteins over %d accessions" % (len(seqs), N), flush=True)
    mpf = {}
    csv.field_size_limit(10 ** 7)
    for r in csv.DictReader(open(os.path.join(PLSDB, "typing.csv"))):
        mpf[r["NUCCORE_ACC"]] = r.get("mpf_type", "") or "-"

    def admit(path):
        with pyhmmer.plan7.HMMFile(path) as fh:
            m = next(iter(fh))
        kw = {"bit_cutoffs": "gathering"} if m.cutoffs.gathering_available() else {"E": 1e-5}
        s = set()
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, **kw):
            for h in top:
                s.add(owner[int(dec(h.name))])
        return s

    exca = admit(os.path.join(HMM, "NF033891.hmm"))
    mt = admit(os.path.join(HMM, "anchor_CagE_TrbE_VirB.hmm"))
    mf = admit(os.path.join(HMM, "mpff_TraC_F_IV.hmm"))
    target = exca - mt - mf
    print("\nExcA-positive: %d plasmids; outside both existing entries: %d (the target)"
          % (len(exca), len(target)), flush=True)
    tc = collections.Counter(mpf.get(x, "-") for x in target)
    print("  MOB-suite calls the target set: %s" % dict(tc.most_common(5)))

    cand = sorted(f[:-4] for f in os.listdir(CJ) if f.startswith("T4SS_I_") and f.endswith(".hmm"))
    sets = {}
    print("\n%-16s %8s %9s %9s   %s"
          % ("profile", "admits", "ExcA rec", "recall%", "impurity composition (MOB-suite)"))
    rows = []
    for c in cand:
        s = admit(os.path.join(CJ, c + ".hmm"))
        sets[c] = s
        rec = len(s & target)
        cnt = collections.Counter(mpf.get(x, "-") for x in s)
        imp = {k: v for k, v in cnt.items() if k != "MPF_I"}
        tot_imp = sum(imp.values())
        comp = ", ".join("%s %d (%.0f%%)" % (k, v, 100.0 * v / max(1, tot_imp))
                         for k, v in sorted(imp.items(), key=lambda kv: -kv[1])[:3])
        print("%-16s %8d %9d %8.1f%%   %s"
              % (c.replace("T4SS_", ""), len(s), rec, 100.0 * rec / max(1, len(target)), comp),
              flush=True)
        rows.append({"profile": c, "n_admitted": len(s), "excA_recovered": rec,
                     "excA_recall_pct": round(100.0 * rec / max(1, len(target)), 2),
                     "purity_mpfi_pct": round(100.0 * cnt.get("MPF_I", 0) / max(1, len(s)), 2),
                     "impurity_untyped": imp.get("-", 0),
                     "impurity_MPF_T": imp.get("MPF_T", 0),
                     "impurity_MPF_F": imp.get("MPF_F", 0),
                     "impurity_total": tot_imp})

    nb = set().union(*[sets[c] for c in NARROW])
    bb = set().union(*[sets[c] for c in BROAD])
    extra = bb - nb
    print("\n=== what are the %d plasmids broad admits and narrow does not? ===" % len(extra))
    ec = collections.Counter(mpf.get(x, "-") for x in extra)
    for k, v in ec.most_common(5):
        print("  MOB-suite %-12s %5d (%.1f%%)" % (k, v, 100.0 * v / max(1, len(extra))))
    print("  of them, ExcA-positive: %d" % len(extra & target))
    nhit = collections.Counter()
    for x in extra:
        nhit[sum(1 for c in cand if x in sets[c])] += 1
    print("  how many T4SS_I_* profiles each carries: %s"
          % ", ".join("%d:%d" % (k, v) for k, v in sorted(nhit.items(), reverse=True)[:8]))
    print("\n  Many I_* profiles each -> real MPF_I systems narrow is missing.")
    print("  One or two only -> fragments or spurious, and narrow is right.")

    print("\n=== I_traU's impurity, examined (misses 90%% by 1.3 points) ===")
    s = sets.get("T4SS_I_traU", set())
    cnt = collections.Counter(mpf.get(x, "-") for x in s)
    imp = {k: v for k, v in cnt.items() if k != "MPF_I"}
    print("  admitted %d, MPF_I %d (%.1f%%), impure %d"
          % (len(s), cnt.get("MPF_I", 0), 100.0 * cnt.get("MPF_I", 0) / max(1, len(s)),
             sum(imp.values())))
    for k, v in sorted(imp.items(), key=lambda kv: -kv[1]):
        print("    %-12s %5d (%.1f%% of impurity)" % (k, v, 100.0 * v / max(1, sum(imp.values()))))
    print("  untyped impurity is NOT cross-class contamination -- it may be MPF_I")
    print("  that MOB-suite failed to type. MPF_F's excluded families were eating")
    print("  6,000-9,000 plasmids of a DIFFERENT class, which is a different thing.")

    dest = os.path.join(PROJ, "data", "anchors", "mpfi_excA_recall.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
