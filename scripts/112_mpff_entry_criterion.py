#!/usr/bin/env python3
"""Measure the MPF_F entry criterion across the WHOLE database.

MPF_T's entry anchor was not chosen from the literature. VirB4/PF03135 was picked
because it admits all four controls at GA and, run across PLSDB, was the most
conservative family that did so. The entry anchor determines what you can never
see, so it has to be measured, not inferred from seed scores.

Two things this must NOT do:

  * measure inside the VirB4-selected cache. An MPF_F entry criterion computed
    there would be measuring VirB4. Hence the full cache: 72,557 accessions,
    7,725,213 CDS.
  * use MOB-suite's mpf_type as ground truth. It reports 11,714 MPF_F plasmids,
    which is an order-of-magnitude reference only -- it is coarser than CONJscan
    and gives no component coordinates. It is reported alongside for orientation,
    never used to define the set.

Reported per family: plasmids admitted, whether all five seeds are admitted, and
the pairwise overlap structure. The criterion is the most conservative family that
still admits every seed.
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
PLSDB = "/global/scratch/users/kh36969/plsdb"
# PLSDB is RefSeq-based; the seed records are INSDC. None of the five INSDC
# accessions is present in PLSDB, so the first run reported 0/5 seeds for all 13
# families -- a failed baseline, not a biological result. Four seeds have RefSeq
# equivalents in PLSDB. SXT is an ICE rather than a plasmid, so PLSDB legitimately
# lacks it and it is injected into the cache (cds_900.faa), exactly as RP4 (a TPA
# record) was injected into the mash sketch.
SEED_ACC = {"F": "NC_002483.1", "R100": "NC_002134.1", "IncC": "NZ_CP033514.1",
            "SXT": "KJ817376.1", "R27": "NC_002305.1"}


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
    fams = [(r["pfam_acc"], r["pfam_name"], r["role"])
            for r in csv.DictReader(lines, delimiter="\t")]
    print("MPF_F anchor families to test as entry criterion: %d" % len(fams))

    seqs, owner, aas = [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                owner.append(nm.split("|")[0]); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    N = len(set(owner))
    print("full cache: %d proteins over %d accessions" % (len(seqs), N), flush=True)

    mpf = {}
    csv.field_size_limit(10 ** 7)
    tp = os.path.join(PLSDB, "typing.csv")
    if os.path.exists(tp):
        for r in csv.DictReader(open(tp)):
            mpf[r["NUCCORE_ACC"]] = r.get("mpf_type", "") or "-"

    admit = {}
    for acc, name, role in fams:
        p = None
        for cand in ("mpff_%s.hmm" % name, "anchor_%s.hmm" % name, "%s.hmm" % name):
            if os.path.exists(os.path.join(HMM, cand)):
                p = os.path.join(HMM, cand); break
        if not p:
            print("  no HMM for %s (%s)" % (name, acc)); continue
        with pyhmmer.plan7.HMMFile(p) as fh:
            model = next(iter(fh))
        s = set()
        for top in pyhmmer.hmmsearch([model], seqs, cpus=a.cpus,
                                     bit_cutoffs="gathering"):
            for h in top:
                s.add(owner[int(dec(h.name))])
        admit[(acc, name)] = s
        seeds = [k for k, v in SEED_ACC.items() if v in s]
        c = collections.Counter(mpf.get(x, "-") for x in s)
        print("  %-9s %-18s admits %6d plasmids (%.1f%%)  seeds %d/5 [%s]  "
              "MOB-suite: MPF_F %d, MPF_T %d, other %d"
              % (acc, name, len(s), 100.0 * len(s) / N, len(seeds),
                 ",".join(sorted(seeds)), c.get("MPF_F", 0), c.get("MPF_T", 0),
                 len(s) - c.get("MPF_F", 0) - c.get("MPF_T", 0)), flush=True)

    print("\n=== candidates that admit ALL FIVE seeds, ranked most conservative ===")
    ok = [(k, v) for k, v in admit.items()
          if all(x in v for x in SEED_ACC.values())]
    for (acc, name), v in sorted(ok, key=lambda kv: len(kv[1])):
        c = collections.Counter(mpf.get(x, "-") for x in v)
        print("  %-9s %-18s %6d plasmids   MPF_F %5d (%.0f%% of its admissions)"
              % (acc, name, len(v), c.get("MPF_F", 0),
                 100.0 * c.get("MPF_F", 0) / len(v)))
    if not ok:
        print("  NONE -- no single family admits all five seeds at GA.")
        print("  Report which seed each family misses; a disjunction may be needed,")
        print("  as VirB4 needed PF27097 as a secondary check on MPF_T.")
        for (acc, name), v in sorted(admit.items(), key=lambda kv: -len(kv[1])):
            miss = [k for k, x in SEED_ACC.items() if x not in v]
            if miss:
                print("    %-18s misses %s" % (name, ",".join(sorted(miss))))

    print("\n=== pairwise overlap among the all-five families (Jaccard) ===")
    keys = [k for k, _ in sorted(ok, key=lambda kv: len(kv[1]))][:8]
    if keys:
        print("      %s" % " ".join("%-10s" % k[1][:10] for k in keys))
        for x in keys:
            cells = []
            for y in keys:
                A, B = admit[x], admit[y]
                cells.append("%-10.3f" % (len(A & B) / max(1, len(A | B))))
            print("  %-10s %s" % (x[1][:10], " ".join(cells)))

    dest = os.path.join(PROJ, "data", "anchors", "mpff_entry_criterion.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pfam_acc", "pfam_name", "n_admitted", "pct_of_db",
                    "seeds_admitted", "mob_MPF_F", "mob_MPF_T", "mob_other"])
        for (acc, name), v in sorted(admit.items(), key=lambda kv: len(kv[1])):
            c = collections.Counter(mpf.get(x, "-") for x in v)
            w.writerow([acc, name, len(v), round(100.0 * len(v) / N, 2),
                        sum(1 for x in SEED_ACC.values() if x in v),
                        c.get("MPF_F", 0), c.get("MPF_T", 0),
                        len(v) - c.get("MPF_F", 0) - c.get("MPF_T", 0)])
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
