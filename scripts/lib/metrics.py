#!/usr/bin/env python3
"""A6 -- one value per (gene_a, gene_b, metric) triple, enforced.

The bug this exists to prevent already happened: F-vs-R100 TraS identity was
stored as both 17.7% and 27.9%. Neither was wrong arithmetic --

    identities 36, alignment length 203, ungapped columns 129
    36/203 = 17.7%      gaps counted in the denominator
    36/129 = 27.9%      gaps excluded

-- but they were stored under one unqualified name, "identity", so the pair
looked like a contradiction. The fix is not to pick one: it is to make the
DEFINITION part of the metric name, so two definitions cannot collide, and to
hard-fail on a genuine duplicate.

Consequence worth carrying: `humbert2019`'s "17%" was matched to the 17.7%
figure. If that paper used the ungapped definition its comparator is 27.9% and
the match was coincidence. Recorded as `provenance_caveat` on that metric.
"""
import csv
import json
import os

# Every permitted metric, with its denominator spelled out. Adding a metric means
# adding a row here; an unregistered name is a hard failure.
DEFINITIONS = {
    "aa_identity_pct__over_alignment_length":
        "100 * identities / total alignment columns, gaps INCLUDED in the "
        "denominator. Biopython PairwiseAligner, BLOSUM62, gap open -11 extend -1, global.",
    "aa_identity_pct__over_ungapped_columns":
        "100 * identities / columns where BOTH sequences have a residue. Same aligner.",
    "aa_identity_pct__local_over_alignment_length":
        "as above but mode=local. Used for homolog anchoring (R27).",
    "window_substitution_enrichment_fold":
        "(substitutions in window / window length) / (substitutions overall / "
        "aligned length). Gapped columns excluded from both.",
    "adjacency_bp":
        "intergenic distance in bp between two features on one accession, "
        "order-independent, exclusive of both features.",
    "hmm_bitscore":
        "HMMER3 full-sequence bit score at the stated threshold policy.",
    "hmm_evalue":
        "HMMER3 full-sequence E-value.",
    "mature_length_aa":
        "precursor length minus the signal peptide, cleavage site stated separately.",
}


class MetricConflict(RuntimeError):
    pass


class MetricStore:
    """Append-only store keyed on (gene_a, gene_b, metric). Duplicates fail."""

    def __init__(self):
        self._d = {}

    def put(self, gene_a, gene_b, metric, value, source, accession_a="",
            accession_b="", provenance_caveat=""):
        if metric not in DEFINITIONS:
            raise MetricConflict(
                "unregistered metric %r. Add it to DEFINITIONS with its "
                "denominator spelled out; an unqualified name is how 17.7%% and "
                "27.9%% ended up as the same 'identity'." % metric)
        key = (gene_a, gene_b, metric)
        if key in self._d:
            prev = self._d[key]
            if abs(float(prev["value"]) - float(value)) > 1e-9:
                raise MetricConflict(
                    "A6 VIOLATION: %s already holds %s (from %s); refusing to "
                    "store %s (from %s). One value per triple. If these are "
                    "genuinely different measurements, they need different "
                    "metric names." % (key, prev["value"], prev["source"],
                                       value, source))
            return
        self._d[key] = {"gene_a": gene_a, "gene_b": gene_b, "metric": metric,
                        "value": value, "source": source,
                        "accession_a": accession_a, "accession_b": accession_b,
                        "definition": DEFINITIONS[metric],
                        "provenance_caveat": provenance_caveat}

    def rows(self):
        return sorted(self._d.values(),
                      key=lambda r: (r["gene_a"], r["gene_b"], r["metric"]))

    def write(self, path):
        cols = ["gene_a", "gene_b", "metric", "value", "accession_a",
                "accession_b", "source", "provenance_caveat", "definition"]
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(self.rows())


def seed_store():
    """The measurements taken so far, each under its qualified metric name."""
    s = MetricStore()
    s.put("F/traS", "R100/traS", "aa_identity_pct__over_alignment_length", 17.7,
          "scripts/45_traT_audit.py-style alignment, 2026-09-01",
          "AP001918.1", "AP000342.1",
          "MARGINAL ALIGNMENT -- do not build an argument on this number. 74 of "
          "203 columns are gapped (36%). The value swings with method: 17.7 "
          "(BLOSUM62 -11/-1 over alignment length), 27.9 (same, ungapped), 25.5 "
          "(BLOSUM45), 40.0 (local, 15-residue block). humbert2019's '17%' was "
          "matched to this definition; if that paper used the ungapped "
          "denominator its comparator is 27.9 and the match is coincidence. "
          "UNRESOLVED, and now demoted to a literature-reconciliation item only: "
          "the R100 position-conserved argument uses the PF10624 HMM result "
          "instead (R100 traS 312.5 at GA, F traS undetected), which is "
          "threshold-anchored and alignment-independent.")
    s.put("F/traS", "R100/traS", "aa_identity_pct__over_ungapped_columns", 27.9,
          "scripts/46_window_null.py, 2026-09-02", "AP001918.1", "AP000342.1")
    s.put("F/traG", "R100/traG", "aa_identity_pct__over_alignment_length", 92.8,
          "2026-09-01", "AP001918.1", "AP000342.1")
    s.put("F/traT", "R100/traT", "aa_identity_pct__over_alignment_length", 99.2,
          "scripts/45_traT_audit.py", "AP001918.1", "AP000342.1")
    s.put("F/traG:610-673", "R100/traG", "window_substitution_enrichment_fold",
          6.88, "scripts/46_window_null.py", "AP001918.1", "AP000342.1",
          "empirical p = 0.022 against the max-window null over 45 tra ortholog "
          "pairs; traB reaches 6.75x. Not significant enough to be a primary claim.")
    s.put("RP4/trbK", "R751/trbK", "aa_identity_pct__over_alignment_length", 37.8,
          "slot family analysis 2026-09-02", "BN000925.1", "NC_001735.4")
    s.put("pKM101/eex", "R388/candidate", "aa_identity_pct__over_alignment_length",
          40.7, "slot family analysis 2026-09-02", "U09868.1", "NC_028464.1")
    s.put("R27/eexA", "pAPEC-O1-R/trhZ", "aa_identity_pct__local_over_alignment_length",
          78.1, "scripts/80-era homolog anchoring", "AF250878.1", "DQ517526.1",
          "cross-Inc-group: R27 is IncHI1, pAPEC-O1-R is IncHI2.")
    s.put("R27/trhG", "pAPEC-O1-R/trhG", "aa_identity_pct__local_over_alignment_length",
          76.1, "homolog anchoring", "AF250878.1", "DQ517526.1",
          "cross-Inc-group, as above.")
    return s


if __name__ == "__main__":
    PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    st = seed_store()
    out = os.path.join(PROJ, "data", "metrics.tsv")
    st.write(out)
    print("registered metric definitions: %d" % len(DEFINITIONS))
    print("stored measurements: %d" % len(st.rows()))
    print("\n  %-22s %-22s %-46s %s"
          % ("gene_a", "gene_b", "metric", "value"))
    for r in st.rows():
        print("  %-22s %-22s %-46s %s"
              % (r["gene_a"], r["gene_b"], r["metric"], r["value"]))
    print("\n=== A6 self-test: a genuine duplicate must hard-fail ===")
    try:
        st.put("F/traS", "R100/traS", "aa_identity_pct__over_alignment_length",
               27.9, "a script that used the other denominator")
        print("  NOT CAUGHT -- A6 is not working")
    except MetricConflict as e:
        print("  CAUGHT: " + str(e).split(". ")[0])
    print("\n=== and an unregistered metric name must hard-fail ===")
    try:
        st.put("X", "Y", "identity", 50.0, "sloppy caller")
        print("  NOT CAUGHT -- A6 is not working")
    except MetricConflict as e:
        print("  CAUGHT: " + str(e).split(".")[0])
    print("\nwrote %s" % os.path.relpath(out, PROJ))
