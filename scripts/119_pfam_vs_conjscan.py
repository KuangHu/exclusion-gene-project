#!/usr/bin/env python3
"""Pfam anchors vs CONJScan profiles, side by side, on ten seeds. Plus relaxase.

Both frozen catalogues are built entirely on Pfam anchors. CONJScan has been on
disk since MPF_T step one but was only ever used for class assignment -- its
component profiles were never used as anchors. MPF_I exposed the cost: TraU, TrbC
and TraY score NOTHING against Pfam-A at GA, and 1460 / 105 / 951 against
CONJScan. The right protein with the wrong model reads as absence.

So the question for the releases: is the two-tier failure structure partly a
library property? VirB5's 27.1% is load-bearing -- it motivates the
`E<=1e-5 AND aa>150` conjunct, the architecture fallback, every
"conditional on VirB5 detectable" label, and the 5,101 vs 7,014 denominator gap.
If T4SS_T_virB5 simply works, most of that scaffolding was built around PF07996.
If it does not, 27.1% is the real ceiling and the scaffolding stays. Either answer
is clean.

Third block, and the more serious gap: T4SS_MOBB and its nine exchangeables are
`presence="mandatory"` in EVERY CONJScan class definition, and NEITHER catalogue
scores a relaxase. Tier A therefore means "has the channel and motor", not "can
transfer DNA". Measured here at seed level; the database-wide count decides
whether that is a new column or a revised tier.

Scores are each model's own, at its own GA. No threshold is transferred between
libraries.
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"

SEEDS = [
    ("MPF_T", "RP4", "RP4__BN000925.1.gb"),
    ("MPF_T", "R751", "R751__NC_001735.4.gb"),
    ("MPF_T", "pKM101", "pKM101__U09868.1.gb"),
    ("MPF_T", "R388", "R388__NC_028464.1.gb"),
    ("MPF_T", "R46", "R46__AY046276.1.gb"),
    ("MPF_F", "F", "F__AP001918.1.gb"),
    ("MPF_F", "R100", "R100__AP000342.1.gb"),
    ("MPF_F", "IncC", "pVCR94__CP033514.1.gb"),
    ("MPF_F", "SXT", "SXT__KJ817376.1.gb"),
    ("MPF_F", "R27", "R27__AF250878.1.gb"),
]
# anchor -> (pfam model file in HMM, conjscan profile name)
PAIRS_T = [("VirB1", "anchor_SLT", "T4SS_T_virB1"),
           ("VirB2", "anchor_TrbC", "T4SS_T_virB2"),
           ("VirB3", "anchor_VirB3", "T4SS_T_virB3"),
           ("VirB4", "anchor_CagE_TrbE_VirB", "T4SS_virb4"),
           ("VirB5", "anchor_T4SS", "T4SS_T_virB5"),
           ("VirB6", "TrbL", "T4SS_T_virB6"),
           ("VirB8", "anchor_VirB8", "T4SS_T_virB8"),
           ("VirB9", "anchor_CagX", "T4SS_T_virB9"),
           ("VirB10", "anchor_TrbI", "T4SS_T_virB10"),
           ("VirB11", "anchor_T2SSE", "T4SS_T_virB11"),
           ("VirD4", "anchor_T4SS-DNA_transf", "T4SS_t4cp1")]
PAIRS_F = [("TraU", "mpff_TraU", "T4SS_F_traU"),
           ("TraH", "mpff_TraH", "T4SS_F_traH"),
           ("TraF", "mpff_TraF", "T4SS_F_traF"),
           ("TraG_N", "mpff_TraG_N", "T4SS_F_traG"),
           ("TrbI", "mpff_TrbI", "T4SS_F_traB"),
           ("TraE", "mpff_TraE", "T4SS_F_traE"),
           ("TrbC_Ftype", "mpff_TrbC_Ftype", "T4SS_F_trbC"),
           ("F_T4SS_TraN", "mpff_F_T4SS_TraN", "T4SS_F_traN"),
           ("TraV", "mpff_TraV", "T4SS_F_traV"),
           ("TraC_F_IV", "mpff_TraC_F_IV", None),
           ("TrwB_AAD_bind", "mpff_TrwB_AAD_bind", "T4SS_t4cp2")]
MOB = ["T4SS_MOBB", "T4SS_MOBC", "T4SS_MOBF", "T4SS_MOBH", "T4SS_MOBP1",
       "T4SS_MOBP2", "T4SS_MOBP3", "T4SS_MOBQ", "T4SS_MOBT", "T4SS_MOBV"]


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=8)
    a = ap.parse_args()
    import pyhmmer
    from Bio import SeqIO
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    data = {}
    for cls, nm, fn in SEEDS:
        p = os.path.join(PROJ, "data", "seed", "genbank", fn)
        if not os.path.exists(p):
            print("  MISSING seed record: %s" % fn); continue
        rec = next(SeqIO.parse(p, "genbank"))
        cds = []
        for f in rec.features:
            if f.type != "CDS":
                continue
            tr = (f.qualifiers.get("translation") or [""])[0]
            if not tr:
                continue
            g = ((f.qualifiers.get("gene") or [""])[0] or
                 (f.qualifiers.get("locus_tag") or [""])[0])
            cds.append((g, tr))
        data[nm] = cds
        print("  %-8s %-8s %4d CDS" % (cls, nm, len(cds)), flush=True)

    def best(path, seed):
        if not path or not os.path.exists(path):
            return None
        with pyhmmer.plan7.HMMFile(path) as fh:
            m = next(iter(fh))
        ga = m.cutoffs.gathering[0] if m.cutoffs.gathering_available() else None
        seqs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=t).digitize(alpha)
                for i, (g, t) in enumerate(data[seed])]
        kw = {"bit_cutoffs": "gathering"} if ga is not None else {"E": 1e-3}
        top_hit = None
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, **kw):
            for h in top:
                i = int(dec(h.name))
                if top_hit is None or h.score > top_hit[1]:
                    top_hit = (data[seed][i][0], h.score, len(data[seed][i][1]))
        return (top_hit, ga)

    rows = []
    for label, pairs, cls in (("MPF_T", PAIRS_T, "MPF_T"), ("MPF_F", PAIRS_F, "MPF_F")):
        seeds = [n for c, n, _ in SEEDS if c == cls and n in data]
        print("\n" + "=" * 108)
        print("%s anchors: Pfam model vs CONJScan profile, each at its own GA" % label)
        print("=" * 108)
        print("%-14s %-9s %-9s %s" % ("anchor", "pfamGA", "cjGA",
              " ".join("%-17s" % s for s in seeds)))
        for anchor, pf, cj in pairs:
            pfp = os.path.join(HMM, pf + ".hmm")
            cjp = os.path.join(CJ, cj + ".hmm") if cj else None
            pga = cga = None
            cells = []
            for s in seeds:
                rp = best(pfp, s); rc = best(cjp, s)
                if rp: pga = rp[1]
                if rc: cga = rc[1]
                pv = "%.0f" % rp[0][1] if rp and rp[0] else "-"
                cv = "%.0f" % rc[0][1] if rc and rc[0] else ("-" if cjp else "n/a")
                cells.append("%-17s" % ("%s|%s" % (pv, cv)))
                rows.append({"class": label, "anchor": anchor, "seed": s,
                             "pfam_model": pf, "pfam_bit": pv,
                             "conjscan_profile": cj or "", "conjscan_bit": cv})
            print("%-14s %-9s %-9s %s"
                  % (anchor, "%.1f" % pga if pga else "-",
                     "%.1f" % cga if cga else "-", " ".join(cells)))
        print("  cells are  pfam_bit|conjscan_bit   ('-' = no hit at that model's GA)")

    print("\n" + "=" * 108)
    print("RELAXASE -- mandatory in every CONJScan class definition, scored by NEITHER catalogue")
    print("=" * 108)
    print("%-12s %s" % ("seed", " ".join("%-11s" % m.replace("T4SS_", "") for m in MOB)))
    for cls, nm, _ in SEEDS:
        if nm not in data:
            continue
        cells = []
        for m in MOB:
            r = best(os.path.join(CJ, m + ".hmm"), nm)
            cells.append("%-11s" % ("%s:%.0f" % (r[0][0][:5], r[0][1]) if r and r[0] else "-"))
        print("%-12s %s" % (nm, " ".join(cells)))
        rows.append({"class": cls, "anchor": "RELAXASE", "seed": nm,
                     "pfam_model": "", "pfam_bit": "",
                     "conjscan_profile": ";".join(
                         "%s" % m.replace("T4SS_", "") for m in MOB
                         if (lambda r: r and r[0])(best(os.path.join(CJ, m + ".hmm"), nm))),
                     "conjscan_bit": ""})
    dest = os.path.join(PROJ, "data", "anchors", "pfam_vs_conjscan.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["class", "anchor", "seed", "pfam_model",
                                           "pfam_bit", "conjscan_profile", "conjscan_bit"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
