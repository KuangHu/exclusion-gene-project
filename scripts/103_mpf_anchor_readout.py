#!/usr/bin/env python3
"""Read the anchor set for an MPF class off its seeds, empirically.

The MPF_T anchor set was not chosen from the literature -- it was read off the
control plasmids by scanning them against all of Pfam-A and recording which family
actually hits each component. That is the only procedure that survives this
project's repeated finding that family NAMES do not imply mechanism (NF033891 is
called surf_exc_IncI1 and hits IncI1 ExcA, which is entry exclusion).

The same procedure, per class. Nothing is transferred from MPF_T: each class gets
its own anchors, because the slot is defined relative to whatever skeleton that
class actually has.

Known starting points for MPF_F, to be confirmed or corrected by the readout:
    VirB6 analogue  PF07916 (TraG_N)
    TraN            PF06986 (F 268.0, IncC 93.3)
Everything else is unmeasured.

Annotated CDS are used rather than called ones: these are curated seed records, so
the gene names are the point -- the readout maps NAME -> FAMILY, which is what an
anchor set is.
"""
import argparse
import collections
import csv
import os
import re
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PFAM = "/global/scratch/users/kh36969/funcannot_dbs/pfam/Pfam-A.hmm"

CLASSES = {
    "MPF_F": [("F", "F__AP001918.1.gb"), ("R100", "R100__AP000342.1.gb"),
              ("IncC_pVCR94", "pVCR94__CP033514.1.gb"), ("SXT", "SXT__KJ817376.1.gb"),
              ("R27", "R27__AF250878.1.gb")],
    # RefSeq, not INSDC: PLSDB is RefSeq-based and A12 flagged 5 of 9 INSDC
    # accessions as absent. pCVM29188_101 (NC_011077.1) resolved at NCBI but is
    # absent from PLSDB and was injected, as SXT and RP4 were.
    "MPF_I": [("R64", "R64__NC_005014.1.gb"), ("R621a", "R621a__NC_015965.1.gb"),
              ("pEK204", "pEK204__NC_013120.1.gb"),
              ("pCVM29188", "pCVM29188_101__NC_011077.1.gb"),
              ("ColIb-P9", "ColIb-P9__NC_002122.1.gb")],
    "MPF_FA": [("ICEBs1", None), ("pLS20", "pLS20__AB615352.1.gb"),
               ("pCF10", "pCF10__AY855841.2.gb"), ("pAD1", "pAD1__CP046109.1.gb")],
    "MPF_T": [("RP4", "RP4__BN000925.1.gb"), ("R751", "R751__NC_001735.4.gb"),
              ("pKM101", "pKM101__U09868.1.gb"), ("R388", "R388__NC_028464.1.gb")],
}
TRA_RE = re.compile(r"^(tra|trb|trh|trw|vir|dot|icm|pil|prg|pcf|con|cwl|ydd)", re.I)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mpf", default="MPF_F")
    ap.add_argument("--cpus", type=int, default=8)
    a = ap.parse_args()
    import pyhmmer
    from Bio import SeqIO
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, meta = [], []
    for name, fn in CLASSES[a.mpf]:
        if not fn:
            print("  %-14s no seed record in repo -- SKIPPED" % name)
            continue
        p = os.path.join(PROJ, "data", "seed", "genbank", fn)
        if not os.path.exists(p):
            print("  %-14s %s MISSING" % (name, fn))
            continue
        rec = next(SeqIO.parse(p, "genbank"))
        n = 0
        for f in rec.features:
            if f.type != "CDS":
                continue
            tr = (f.qualifiers.get("translation") or [None])[0]
            if not tr:
                continue
            gene = ((f.qualifiers.get("gene") or [""])[0] or
                    (f.qualifiers.get("locus_tag") or [""])[0])
            prod = (f.qualifiers.get("product") or [""])[0]
            meta.append({"seed": name, "gene": gene, "product": prod,
                         "aa": len(tr), "start": int(f.location.start) + 1,
                         "end": int(f.location.end),
                         "strand": 1 if f.location.strand == 1 else -1})
            seqs.append(pyhmmer.easel.TextSequence(
                name=str(len(seqs)).encode(), sequence=tr).digitize(alpha))
            n += 1
        print("  %-14s %-28s %d annotated CDS" % (name, fn, n), flush=True)
    print("total CDS to scan: %d against all of Pfam-A" % len(seqs), flush=True)

    hits = collections.defaultdict(list)
    with pyhmmer.plan7.HMMFile(PFAM) as fh:
        for top in pyhmmer.hmmsearch(fh, seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            q = dec(top.query.name)
            acc = dec(top.query.accession) if top.query.accession else "-"
            for h in top:
                i = int(dec(h.name))
                hits[i].append((q, acc.split(".")[0], h.score))
    print("proteins with >=1 Pfam family at GA: %d/%d" % (len(hits), len(seqs)))

    rows = []
    for i, m in enumerate(meta):
        fams = sorted(hits.get(i, []), key=lambda x: -x[2])
        rows.append(dict(m, n_families=len(fams),
                         families=";".join("%s(%s):%.1f" % f for f in fams[:4])))
    dest = os.path.join(PROJ, "data", "anchors", "readout_%s.tsv" % a.mpf)
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    print("\n=== TRANSFER GENES and the families that actually hit them ===")
    print("%-14s %-10s %5s  %-42s %s" % ("seed", "gene", "aa", "top families (GA)", "product"))
    for r in rows:
        if not TRA_RE.match(r["gene"] or "") and "conjug" not in (r["product"] or "").lower():
            continue
        print("%-14s %-10s %5d  %-42s %s"
              % (r["seed"], (r["gene"] or "-")[:10], r["aa"],
                 (r["families"] or "NO FAMILY AT GA")[:42], (r["product"] or "")[:34]))

    print("\n=== CANDIDATE ANCHORS: families hitting >=3 seeds ===")
    byfam = collections.defaultdict(set)
    for i, m in enumerate(meta):
        for q, acc, sc in hits.get(i, []):
            byfam[(q, acc)].add(m["seed"])
    nseed = len({m["seed"] for m in meta})
    for (q, acc), s in sorted(byfam.items(), key=lambda kv: -len(kv[1])):
        if len(s) >= 3:
            print("  %-24s %-10s %d/%d seeds: %s" % (q, acc, len(s), nseed, ",".join(sorted(s))))
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
