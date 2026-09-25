#!/usr/bin/env python3
"""Build the first negative and ambiguous rows.

These are not padding. docs/heuristic_precision.md measured the discovery
heuristic at 11% precision against 53 in-window decoys, and the two F decoys
below are the sharpest of them: same operon, same strand, inside the +/-6 kb
window, in the length band. Nothing in the implemented filter rejects them.

RP4 trbJ is the ambiguous case and is more informative than any random negative:
same operon, immediately adjacent to the true exclusion gene, partial published
evidence, and an unresolved conflict between reports.
"""
import csv, os, sys
from Bio import SeqIO
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
import ncbi

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(PROJ, "data", "negatives")


def grab(path, genes):
    out = {}
    rec = next(SeqIO.parse(path, "genbank"))
    for f in rec.features:
        if f.type != "CDS":
            continue
        g = ((f.qualifiers.get("gene") or [""])[0]).strip()
        if g not in genes:
            continue
        tr = (f.qualifiers.get("translation") or [""])[0]
        out[g] = {
            "acc": rec.id,
            "pid": (f.qualifiers.get("protein_id") or [""])[0],
            "start": int(f.location.start) + 1,
            "end": int(f.location.end),
            "strand": "+" if f.location.strand == 1 else "-",
            "aa": len(tr),
            "sha": ncbi.aa_sha256(tr),
        }
    return out


def gap(a, b):
    """bp between two features, order-independent."""
    if b["start"] > a["end"]:
        return b["start"] - a["end"] - 1
    return a["start"] - b["end"] - 1


NOTES = {
    "trbF": "In the F tra operon, same strand as traS, 126 aa, inside the +/-6 kb "
            "window around traG. NOT an exclusion gene. Carries a confident "
            "VirB5-family assignment -- which is exactly the filter conjunct that "
            "should reject it and is not yet implemented.",
    "trbH": "In the F tra operon, same strand, 239 aa, inside the window. NOT an "
            "exclusion gene. Routinely omitted from simplified tra gene maps, which "
            "is why it was never considered as a decoy.",
}

TRBJ_NOTE = (
    "THE HARD CASE. Immediately upstream of trbK in the same trb operon. Lessl and "
    "Lyras report low-level entry exclusion activity for IncP-alpha TrbJ; Haase's "
    "data conflict, for reasons not established. UniProt names it gene 'eexA' / "
    "'Entry exclusion protein A' (Q79AS0), which collides with R27 eexA. More "
    "informative than a random negative: same operon, adjacent to the real "
    "exclusion gene, partial evidence, unresolved literature. Any scan that cannot "
    "separate trbJ from trbK has not solved the problem."
)


def main():
    os.makedirs(OUT, exist_ok=True)
    F = grab(os.path.join(PROJ, "data", "seed", "genbank", "F__AP001918.1.gb"),
             {"trbH", "trbF", "traG", "traS"})
    R = grab(os.path.join(PROJ, "data", "seed", "genbank", "RP4__BN000925.1.gb"),
             {"trbJ", "trbK", "trbL"})
    rows = []
    for g in ("trbF", "trbH"):
        x = F[g]
        rows.append({
            "row_id": "F_" + g, "element_id": "F", "gene": g,
            "role": "negative", "row_class": "in_operon_decoy",
            "accession": x["acc"], "protein": x["pid"],
            "coords": "%s:%d..%d" % (x["acc"], x["start"], x["end"]),
            "strand": x["strand"], "aa_len": x["aa"], "aa_sha256": x["sha"],
            "dist_to_virb6_bp": gap(x, F["traG"]),
            "dist_to_true_excl_bp": gap(x, F["traS"]),
            "same_operon": "yes", "evidence_tier": "negative_control",
            "literature_conflict": "no", "notes": NOTES[g],
        })
    x = R["trbJ"]
    rows.append({
        "row_id": "RP4_trbJ", "element_id": "RP4", "gene": "trbJ",
        "role": "ambiguous", "row_class": "disputed_exclusion_activity",
        "accession": x["acc"], "protein": x["pid"],
        "coords": "%s:%d..%d" % (x["acc"], x["start"], x["end"]),
        "strand": x["strand"], "aa_len": x["aa"], "aa_sha256": x["sha"],
        "dist_to_virb6_bp": gap(x, R["trbL"]),
        "dist_to_true_excl_bp": gap(x, R["trbK"]),
        "same_operon": "yes", "evidence_tier": "E4_disputed",
        "literature_conflict": "yes", "notes": TRBJ_NOTE,
    })
    cols = ["row_id", "element_id", "gene", "role", "row_class", "accession",
            "protein", "coords", "strand", "aa_len", "aa_sha256",
            "dist_to_virb6_bp", "dist_to_true_excl_bp", "same_operon",
            "evidence_tier", "literature_conflict", "notes"]
    path = os.path.join(OUT, "negatives_ambiguous.tsv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("  %-10s %-10s %-6s %4daa  d(VirB6)=%6dbp  d(true excl)=%6dbp  %s"
              % (r["row_id"], r["role"], r["gene"], r["aa_len"],
                 r["dist_to_virb6_bp"], r["dist_to_true_excl_bp"],
                 r["evidence_tier"]))
    print("\nwrote %s (%d rows)" % (os.path.relpath(path, PROJ), len(rows)))


if __name__ == "__main__":
    main()
