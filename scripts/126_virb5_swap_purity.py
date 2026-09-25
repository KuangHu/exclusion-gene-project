#!/usr/bin/env python3
"""Purity check before the VirB5 disjunction. Do not adopt T4SS_T_virB5 blind.

The motivation is a correctness failure, not an improvement: RP4 and R751 have
slot_ready = 0 in release v1 although both are eex_occupied. slot_ready requires
PF07996 to detect VirB5, and PF07996 scores their TrbJ at 17.7 and 22.2 against
GA 24.7. T4SS_T_virB5 scores the same proteins at 106 and 122. The two plasmids
whose TrbK DEFINED the slot are excluded by the criterion built from them.

But T4SS_T_virB5's GA is 12.6 against PF07996's 24.7 -- far more permissive across
72,558 accessions. That is the shape of the phmmer recall that had to be
retracted: more sensitive, completeness inflated, and the precision control blind
where it mattered. So two checks first.

  (1) LENGTH BAND. VirB5 members measured at 200-288 aa (n=237 of the GA set,
      95.9% within band). What fraction of the NEW admissions fall in it? This is
      the only control here that does not depend on third-party annotation, and it
      is the same test used on the 454 recovered VirB5 in the candidate pool.

  (2) CROSS-CLASS. How many new admissions are already admitted by the MPF_F entry
      criterion (PF11130)? A VirB5 call landing on an MPF_F plasmid is
      contamination.

The swap, if adopted, is a DISJUNCTION -- `PF07996 (E<=1e-5 AND aa>150) OR
T4SS_T_virB5` -- never a replacement. The ten-seed comparison showed neither
library dominates: CONJScan misses TraG_N on IncC (375) and R27 (103), which is
MPF_F's slot anchor, and has no TraC_F_IV profile at all, which is MPF_F's entry.
Replacing Pfam with CONJScan would break MPF_F.
"""
import argparse
import collections
import csv
import os
import statistics
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
BAND = (200, 288)
CONTROLS = {"RP4": "BN000925.1", "R751": "NC_001735.4", "R388": "NC_028464.1"}


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
    N = len(set(owner))
    print("cache: %d proteins over %d accessions" % (len(seqs), N), flush=True)

    def scan(path, ev=None, minaa=0):
        with pyhmmer.plan7.HMMFile(path) as fh:
            m = next(iter(fh))
        kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
        prot, plas = set(), set()
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, **kw):
            for h in top:
                i = int(dec(h.name))
                if minaa and len(aas[i]) <= minaa:
                    continue
                prot.add(i); plas.add(owner[i])
        return prot, plas

    pf_p, pf_a = scan(os.path.join(HMM, "anchor_T4SS.hmm"), 1e-5, 150)
    cj_p, cj_a = scan(os.path.join(CJ, "T4SS_T_virB5.hmm"))
    mf_p, mf_a = scan(os.path.join(HMM, "mpff_TraC_F_IV.hmm"))
    mt_p, mt_a = scan(os.path.join(HMM, "anchor_CagE_TrbE_VirB.hmm"))
    print("PF07996 (E<=1e-5 AND aa>150): %d plasmids" % len(pf_a))
    print("T4SS_T_virB5 (GA 12.6):       %d plasmids" % len(cj_a))
    print("union (the disjunction):      %d plasmids" % len(pf_a | cj_a), flush=True)

    print("\n=== do the defining controls come back? ===")
    for nm, acc in CONTROLS.items():
        print("  %-6s %-14s PF07996 %-4s  CONJScan %-4s"
              % (nm, acc, "yes" if acc in pf_a else "NO", "yes" if acc in cj_a else "NO"))

    new_p = cj_p - pf_p
    new_a = cj_a - pf_a
    print("\n=== (1) LENGTH BAND, the annotation-free control ===")
    L = sorted(len(aas[i]) for i in new_p)
    old = sorted(len(aas[i]) for i in pf_p)
    print("  PF07996 proteins   n=%d  median %d  in %d-%d aa: %.1f%%"
          % (len(old), statistics.median(old),
             BAND[0], BAND[1], 100.0 * sum(1 for x in old if BAND[0] <= x <= BAND[1]) / len(old)))
    if L:
        inband = sum(1 for x in L if BAND[0] <= x <= BAND[1])
        print("  NEW from CONJScan  n=%d  median %d  in %d-%d aa: %.1f%%"
              % (len(L), statistics.median(L), BAND[0], BAND[1], 100.0 * inband / len(L)))
        h = collections.Counter(min(9, x // 50) for x in L)
        mx = max(h.values())
        for b in range(10):
            lab = "%d-%d" % (b * 50, b * 50 + 49) if b < 9 else ">=450"
            star = "  <- VirB5 band" if b in (4, 5) else ""
            print("    %-10s %6d %s%s" % (lab, h.get(b, 0), "#" * int(40 * h.get(b, 0) / mx), star))
    print("\n=== (2) CROSS-CLASS ===")
    print("  new plasmids also admitted by MPF_F entry (PF11130): %d (%.1f%%)"
          % (len(new_a & mf_a), 100.0 * len(new_a & mf_a) / max(1, len(new_a))))
    print("  new plasmids also admitted by MPF_T entry (PF03135): %d (%.1f%%)"
          % (len(new_a & mt_a), 100.0 * len(new_a & mt_a) / max(1, len(new_a))))
    print("  (a VirB5 call on an MPF_F plasmid is contamination; on an MPF_T")
    print("   plasmid it is exactly what should happen)")

    print("\n=== effect on slot_ready, if adopted as a disjunction ===")
    cat = list(csv.DictReader(open(os.path.join(PROJ, "data", "release", "v1",
                                                "catalogue_MPF_T_v1.tsv")), delimiter="\t"))
    n = len(cat)
    was = sum(1 for r in cat if r["slot_ready__virb5_virb6"] == "1")
    b6 = sum(1 for r in cat if r.get("VirB6", "none") != "none")
    now = sum(1 for r in cat if r["accession"] in (pf_a | cj_a)
              and r.get("VirB6", "none") != "none"
              and r["architecture"] != "uncallable")
    print("  release v1 slot_ready      %d / %d = %.1f%%" % (was, n, 100.0 * was / n))
    print("  with the disjunction       %d / %d = %.1f%%" % (now, n, 100.0 * now / n))
    print("  (VirB6 present on %d; slot_ready also needs a callable architecture)" % b6)


if __name__ == "__main__":
    main()
