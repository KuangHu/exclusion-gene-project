#!/usr/bin/env python3
"""Does the two-tier anchor structure survive the CONJScan disjunction?

VirB5 fell from 27.1% to 3.1% under `PF07996 OR T4SS_T_virB5`, both sides carrying
the same aa>150 conjunct, and it passed the controls the retracted phmmer recall
failed: new hits median 258 aa against the family's own 238 (that recall gave 816
against 238), 88.3% in the measured length band, 6.1% cross-class.

That moves VirB5 from the peripheral tier into the core. So the remaining tier is:

    VirB1  19.5%   SLT -- already shown to be REAL absence: 1,118 of 1,452
                   VirB1-negative plasmids carry all seven core anchors
    VirB2  20.9%   pilin
    VirB3   9.0%

VirB2 and VirB3 are pilus/surface components like VirB5, and CONJScan has profiles
for all nine T4SS_T_*. If they collapse the same way, the two-tier structure
largely does not exist -- and with it the biological reading that peripheral
components are poorly detected because they are under diversifying selection.
That reading is currently in the docs and would have to be withdrawn.

Same controls for every anchor, no exceptions:
    length band   robust band (Q1-1.5*IQR, Q3+1.5*IQR) from the Pfam set itself
    cross-class   new admissions also taken by the MPF_F entry criterion
    conjunct      aa>150 is VirB5-specific and is NOT transferred; each anchor's
                  band is measured, not assumed

VirB1 is expected NOT to move much, since its absence was independently confirmed.
It is tested anyway -- an expectation is not a control.
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
PAIRS = [("VirB1", "anchor_SLT", "T4SS_T_virB1", None, 0),
         ("VirB2", "anchor_TrbC", "T4SS_T_virB2", None, 0),
         ("VirB3", "anchor_VirB3", "T4SS_T_virB3", None, 0),
         ("VirB5", "anchor_T4SS", "T4SS_T_virB5", None, 150),
         ("VirB6", "TrbL", "T4SS_T_virB6", None, 0)]


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
    print("cache: %d proteins" % len(seqs), flush=True)

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

    _, mf_a = scan(os.path.join(HMM, "mpff_TraC_F_IV.hmm"))
    _, mt_a = scan(os.path.join(HMM, "anchor_CagE_TrbE_VirB.hmm"))
    N = len(mt_a)        # the VirB4+ set is the MPF_T denominator
    print("MPF_T set (PF03135): %d ; MPF_F set (PF11130): %d\n" % (N, len(mf_a)), flush=True)

    rows = []
    union_sets, pfam_sets = {}, {}
    print("%-8s %9s %9s %9s %9s %9s %9s %8s"
          % ("anchor", "pfam", "conjscan", "union", "GA fail", "new fail", "in band", "xMPF_F"))
    for name, pf, cj, ev, minaa in PAIRS:
        pp, pa = scan(os.path.join(HMM, pf + ".hmm"), ev, minaa)
        cjp = os.path.join(CJ, cj + ".hmm")
        if not os.path.exists(cjp):
            print("%-8s (no CONJScan profile)" % name); continue
        cp, ca = scan(cjp, None, minaa)          # same conjunct on both sides
        L = sorted(len(aas[i]) for i in pp)
        q1, q3 = L[len(L)//4], L[3*len(L)//4]
        iqr = q3 - q1
        lo, hi = max(1, int(q1 - 1.5 * iqr)), int(q3 + 1.5 * iqr)
        new_p, new_a = cp - pp, ca - pa
        nl = sorted(len(aas[i]) for i in new_p)
        inband = (100.0 * sum(1 for x in nl if lo <= x <= hi) / len(nl)) if nl else float("nan")
        union = pa | ca
        union_sets[name] = union; pfam_sets[name] = pa
        xf = len(new_a & mf_a)
        print("%-8s %9d %9d %9d %8.1f%% %8.1f%% %8.1f%% %7.1f%%"
              % (name, len(pa & mt_a), len(ca & mt_a), len(union & mt_a),
                 100.0 * (N - len(pa & mt_a)) / N, 100.0 * (N - len(union & mt_a)) / N,
                 inband, 100.0 * xf / max(1, len(new_a))), flush=True)
        rows.append({"anchor": name, "pfam_plasmids": len(pa & mt_a),
                     "conjscan_plasmids": len(ca & mt_a), "union_plasmids": len(union & mt_a),
                     "ga_failure_pct": round(100.0 * (N - len(pa & mt_a)) / N, 2),
                     "union_failure_pct": round(100.0 * (N - len(union & mt_a)) / N, 2),
                     "pfam_band": "%d-%d" % (lo, hi),
                     "pfam_median_aa": statistics.median(L),
                     "new_median_aa": statistics.median(nl) if nl else "",
                     "new_in_band_pct": round(inband, 1) if nl else "",
                     "new_crossclass_mpff_pct": round(100.0 * xf / max(1, len(new_a)), 1)})
    print("\n  failure rates are over the %d-plasmid MPF_T set (PF03135-admitted)" % N)
    print("  'in band' uses each anchor's OWN robust band from its Pfam hits --")
    print("  VirB5's aa>150 conjunct is not transferred to the others.")
    print("\n=== does the two-tier structure survive? ===")
    print("  core reference (unchanged): VirD4 4.1, VirB10 4.7, VirB11 4.7, VirB8 4.9, VirB9 4.9")
    print("  A verdict is read ONLY from a disjunction that PASSED both controls.")
    print("  band >=80%% in the anchor's own robust band; cross-class <10%% MPF_F.")
    print("  %-8s %9s %9s %9s  %s" % ("anchor", "GA fail", "union", "admissible", "verdict"))
    for r in rows:
        band = r["new_in_band_pct"]
        xf = r["new_crossclass_mpff_pct"]
        ok_band = isinstance(band, float) and band >= 80.0
        ok_xf = xf < 10.0
        adm = ok_band and ok_xf
        if not adm:
            why = []
            if not ok_band: why.append("band %.1f%%" % band if isinstance(band, float) else "band n/a")
            if not ok_xf: why.append("cross-class %.1f%%" % xf)
            verdict = "REJECTED (%s) -- union NOT read; GA stands at %.1f%%" % (
                "; ".join(why), r["ga_failure_pct"])
            shown = r["ga_failure_pct"]
        else:
            verdict = ("still peripheral" if r["union_failure_pct"] > 8 else "now core-level")
            shown = r["union_failure_pct"]
        r["disjunction_admissible"] = int(adm)
        r["verdict"] = verdict
        r["failure_pct_used"] = round(shown, 2)
        print("  %-8s %8.1f%% %8.1f%% %9s  %s"
              % (r["anchor"], r["ga_failure_pct"], r["union_failure_pct"],
                 "yes" if adm else "NO", verdict))
    # --- effect on slot_ready, using ONLY admissible disjunctions ---
    import csv as _csv
    cat = list(_csv.DictReader(open(os.path.join(PROJ, "data/release/v1/catalogue_MPF_T_v1.tsv")),
                               delimiter="\t"))
    n = len(cat)
    was = sum(1 for r in cat if r.get("slot_ready__virb5_virb6") == "1")
    adm = {r["anchor"]: r for r in rows if r.get("disjunction_admissible")}
    print("\n=== effect on slot_ready (admissible disjunctions only) ===")
    print("  release v1 slot_ready       %d / %d = %.1f%%" % (was, n, 100.0 * was / n))
    if "VirB5" in adm:
        b5 = union_sets["VirB5"] & mt_a
        b6 = union_sets["VirB6"] & mt_a if "VirB6" in adm else pfam_sets["VirB6"] & mt_a
        both = len(b5 & b6)
        print("  VirB5 admissible: union %d ; VirB6 side %s" % (len(b5),
              "union" if "VirB6" in adm else "Pfam GA only (disjunction rejected)"))
        print("  VirB5 AND VirB6 present    %d / %d = %.1f%%" % (both, n, 100.0 * both / n))
        print("  (slot_ready also requires a callable architecture; this is the ceiling)")
    dest = os.path.join(PROJ, "data", "anchors", "peripheral_anchor_swap.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
