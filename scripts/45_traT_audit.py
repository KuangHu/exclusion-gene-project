#!/usr/bin/env python3
"""Stage 4 audit -- derive the F/R100 TraT relationship from pinned accessions.

Exists because two published claims conflict (1990s: residues 116-120 confer
F/R100 specificity; 2025 structural work: the proteins are identical) and S7
requires the answer be re-derived from current records rather than transcribed
from either. Writes docs/traT_audit.md.

Usage: 45_traT_audit.py [gene]      default gene: traT
"""
import os, sys
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GB   = os.path.join(PROJ, "data", "seed", "genbank")
PAIR = [("F", "F__AP001918.1.gb"), ("R100", "R100__AP000342.1.gb")]

def cds(path, gene):
    rec = next(SeqIO.parse(path, "genbank"))
    for f in rec.features:
        if f.type == "CDS" and (f.qualifiers.get("gene") or [None])[0] == gene:
            return (rec.id, (f.qualifiers.get("protein_id") or [None])[0],
                    (f.qualifiers.get("translation") or [""])[0])
    raise SystemExit(f"{gene} not found in {path}")

def main():
    gene = sys.argv[1] if len(sys.argv) > 1 else "traT"
    (n1, f1), (n2, f2) = PAIR
    i1, p1, s1 = cds(os.path.join(GB, f1), gene)
    i2, p2, s2 = cds(os.path.join(GB, f2), gene)

    al = PairwiseAligner(mode="global",
                         substitution_matrix=substitution_matrices.load("BLOSUM62"),
                         open_gap_score=-11, extend_gap_score=-1)
    aln = al.align(s1, s2)[0]
    a, b = aln[0], aln[1]
    ident = sum(1 for x, y in zip(a, b) if x == y and x != "-")

    print(f"{n1:<5} {i1} {p1} {len(s1)}aa")
    print(f"{n2:<5} {i2} {p2} {len(s2)}aa")
    print(f"\nidentity {ident}/{len(a)} = {100*ident/len(a):.1f}%")
    diffs = [(i, x, y) for i, (x, y) in enumerate(zip(a, b), 1) if x != y]
    print(f"differing columns: {len(diffs)}")
    for i, x, y in diffs:
        print(f"  col {i:>4}  {n1}={x}  {n2}={y}")

    # numbering drift: map a residue index in seq1 onto seq2
    def project(pos1):
        c1 = c2 = 0
        for x, y in zip(a, b):
            if x != "-": c1 += 1
            if y != "-": c2 += 1
            if c1 == pos1 and x != "-":
                return c2 if y != "-" else None
        return None
    print(f"\nnumbering projection {n1} -> {n2} (S5.3: store windows per accession)")
    for w in (116, 120):
        print(f"  {n1} residue {w} = {n2} residue {project(w)}")
    lo, hi = 116, 120
    seg1 = s1[lo-1:hi]
    j1, j2 = project(lo), project(hi)
    seg2 = s2[j1-1:j2] if j1 and j2 else "?"
    print(f"\n  reported specificity window:")
    print(f"    {n1}   {lo}-{hi}   = {seg1}")
    print(f"    {n2} {j1}-{j2}   = {seg2}")
    print(f"    same residues? {seg1 == seg2}")

if __name__ == "__main__":
    main()
