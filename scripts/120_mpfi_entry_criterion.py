#!/usr/bin/env python3
"""MPF_I entry criterion, measured across the whole database. Same procedure as
MPF_T and MPF_F -- no new steps.

All 17 T4SS_I_* profiles are tested, not just the favourite. MPF_F is why: testing
all 13 there is what exposed the four cross-class families (41-57% purity) that
had to be excluded, and the PF19044-vs-PF11130 domain problem -- two domains of
one protein admitting 18,223 at 57% versus 11,120 at 93%. Testing only the
expected winner would have missed both.

I_traU is the favourite on the seed evidence (1455-1461 across all five, GA 413.1)
and it is the VirB4-family ATPase, converging with MPF_T VirB4, MPF_F TraC and
MPF_FA ConE. CONJScan's own T4SS_typeI definition makes T4SS_virb4 mandatory with
T4SS_I_traU as its exchangeable. But favourites get measured like everything else.

Expected to be EXCLUDED on purity, for the same reason TrbI and TrwB_AAD_bind were
excluded from MPF_F: t4cp1/t4cp2 and the MOB relaxases are pan-conjugative
components and will admit the MPF_T and MPF_F sets wholesale.

Scale reference: the Eex census found NF033891 (ExcA) on 2,293 plasmids of which
only 28 are admitted by the existing two entry criteria, so the target population
is ~2,265. A candidate admitting far more than that probably has a purity problem;
far fewer, a detection problem.

MOB-suite mpf_type is a reference, never ground truth -- MPF_T used CONJscan, and
MOB-suite is coarser with no component coordinates.
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node, assert_keys_present

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
PLSDB = "/global/scratch/users/kh36969/plsdb"
SEEDS = {"R64": "NC_005014.1", "R621a": "NC_015965.1", "pEK204": "NC_013120.1",
         "pCVM29188_101": "NC_011077.1", "ColIb-P9": "NC_002122.1"}
MPFT_ENTRY = ("anchor_CagE_TrbE_VirB", HMM)      # PF03135
MPFF_ENTRY = ("mpff_TraC_F_IV", HMM)             # PF11130


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

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
    accs = set(owner)
    N = len(accs)
    print("full cache: %d proteins over %d accessions" % (len(seqs), N), flush=True)
    assert_keys_present(SEEDS, accs, "MPF_I seed accessions", "the full CDS cache")
    print("A12: all five seeds present", flush=True)

    mpf = {}
    csv.field_size_limit(10 ** 7)
    for r in csv.DictReader(open(os.path.join(PLSDB, "typing.csv"))):
        mpf[r["NUCCORE_ACC"]] = r.get("mpf_type", "") or "-"

    def admit(path):
        if not os.path.exists(path):
            return None
        with pyhmmer.plan7.HMMFile(path) as fh:
            m = next(iter(fh))
        kw = {"bit_cutoffs": "gathering"} if m.cutoffs.gathering_available() else {"E": 1e-5}
        s = set()
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, **kw):
            for h in top:
                s.add(owner[int(dec(h.name))])
        return s

    other = {}
    for lab, (mod, d) in (("MPF_T", MPFT_ENTRY), ("MPF_F", MPFF_ENTRY)):
        other[lab] = admit(os.path.join(d, mod + ".hmm")) or set()
        print("  %s entry admits %d" % (lab, len(other[lab])), flush=True)

    cand = sorted(f[:-4] for f in os.listdir(CJ)
                  if f.startswith(("T4SS_I_", "T4SS_t4cp", "T4SS_MOB")) and f.endswith(".hmm"))
    print("\ncandidates to test as entry criterion: %d\n" % len(cand), flush=True)
    print("%-16s %8s %6s %8s %8s %8s %8s"
          % ("profile", "admits", "seeds", "purity", "xMPF_T", "xMPF_F", "%DB"))
    sets, rows = {}, []
    for c in cand:
        s = admit(os.path.join(CJ, c + ".hmm"))
        if s is None:
            continue
        sets[c] = s
        nseed = sum(1 for v in SEEDS.values() if v in s)
        cnt = collections.Counter(mpf.get(x, "-") for x in s)
        pur = 100.0 * cnt.get("MPF_I", 0) / max(1, len(s))
        xt, xf = len(s & other["MPF_T"]), len(s & other["MPF_F"])
        print("%-16s %8d %5d/5 %7.1f%% %8d %8d %7.1f%%"
              % (c.replace("T4SS_", ""), len(s), nseed, pur, xt, xf, 100.0 * len(s) / N),
              flush=True)
        rows.append({"profile": c, "n_admitted": len(s), "seeds_admitted": nseed,
                     "purity_mpfi_pct": round(pur, 2), "cross_MPF_T": xt,
                     "cross_MPF_F": xf, "pct_of_db": round(100.0 * len(s) / N, 2)})

    ok = [r for r in rows if r["seeds_admitted"] == 5 and r["purity_mpfi_pct"] >= 90]
    print("\n=== admit 5/5 seeds AND >=90%% MPF_I purity (the MPF_F threshold) ===")
    for r in sorted(ok, key=lambda x: x["n_admitted"]):
        print("  %-16s %6d admitted, %.1f%% pure, cross MPF_T %d / MPF_F %d"
              % (r["profile"].replace("T4SS_", ""), r["n_admitted"],
                 r["purity_mpfi_pct"], r["cross_MPF_T"], r["cross_MPF_F"]))
    if not ok:
        print("  NONE -- report which criterion each candidate fails")
        for r in sorted(rows, key=lambda x: -x["purity_mpfi_pct"])[:8]:
            print("    %-16s seeds %d/5, purity %.1f%%"
                  % (r["profile"].replace("T4SS_", ""), r["seeds_admitted"],
                     r["purity_mpfi_pct"]))
    keys = [r["profile"] for r in sorted(ok, key=lambda x: x["n_admitted"])][:8]
    if keys:
        print("\n=== pairwise Jaccard among the qualifying candidates ===")
        print("      %s" % " ".join("%-11s" % k.replace("T4SS_", "")[:11] for k in keys))
        for x in keys:
            print("  %-11s %s" % (x.replace("T4SS_", "")[:11],
                  " ".join("%-11.3f" % (len(sets[x] & sets[y]) / max(1, len(sets[x] | sets[y])))
                           for y in keys)))
    print("\n  scale reference: ~2,265 plasmids carry NF033891/ExcA outside both")
    print("  existing entry criteria -- the population MPF_I should recover.")
    dest = os.path.join(PROJ, "data", "anchors", "mpfi_entry_criterion.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
