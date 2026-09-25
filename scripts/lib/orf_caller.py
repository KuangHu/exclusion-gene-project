#!/usr/bin/env python3
"""THE ORF/CDS caller for this project. Pyrodigal in anon (metagenomic) mode.

One caller, one configuration. Chosen 2026-09-02 on the evidence in
docs/orf_caller_benchmark/. Do not add alternatives here -- if a second caller is
ever needed, that is a decision to re-open with new evidence, not a parameter.

WHY ANON AND NOT NORMAL
  Normal (trained) mode has a HARD FLOOR, not graceful degradation:
      "sequence must be at least 20000 characters"
  It raises outright on U09868.1 (12,416 bp) and J01566.1 (6,646 bp). Fragment
  records and small plasmids are a routine part of this project's input, so
  normal mode would fail on a subset of any batch run -- silently, if the
  exception were swallowed. `assert_anon_only()` below exists to stop that.
  Normal mode also scored worst on recall in the benchmark (7/10 vs 9/10).

WHY NOT SIX-FRAME AS THE PRIMARY LAYER
  The go/no-go test (docs/slot_hypothesis_GO.md) found the exclusion gene is a
  NORMAL SKELETON CALL in all four MPF_T controls -- Pyrodigal anon found trbK
  (RP4), trbK (R751), eex (pKM101) and the R388 candidate as ordinary genes. The
  canonical slot does not need enumeration.

  KNOWN LIMITS, accepted deliberately:
    * one ORF per stop codon, so a same-stop second start is not expressible.
      R64 excA/excB is the case (shared stop at AP005147.1:76983, starts 219 nt
      apart). Pyrodigal returns excA and cannot return excB.
    * a gene fully nested in another frame inside a called gene is missed.
      ColE1 mob4 inside mob3 is the case.
    * boundary-case recall is ~50% (2/4) on genes that a production annotation
      pipeline missed, versus 90% on genes one had already found.
  These are real and documented. They are accepted because the primary use --
  skeleton definition and slot occupancy -- does not depend on them, and because
  a single deterministic caller is worth more here than marginal recall: on two
  deposits of one identical sequence, every caller tested was byte-identical
  while the GenBank annotations differed by 4/47 CDS.

CIRCULARITY CAVEAT
  Pyrodigal has no notion of a circular replicon. With closed=True a gene
  spanning the origin is not called at all; with closed=False it is reported as
  two partial genes at the two ends. The benchmark used closed=True and the
  default here matches it. For a gene of interest near coordinate 1 or near the
  end of a plasmid, re-run on a rotated sequence rather than trusting this.
"""
from Bio.Seq import Seq

# The frozen configuration. Changing any of these invalidates
# docs/orf_caller_benchmark/ and must be justified there.
ANON = True          # pyrodigal -p anon / meta mode
MIN_GENE = 90        # nt; Prodigal's own default. 60 gave identical results
                     # on the benchmark set, so 90 is kept as the documented one.
CLOSED = True        # see CIRCULARITY CAVEAT
TRANSL_TABLE = 11


class CallerMisconfigured(RuntimeError):
    pass


def assert_anon_only(meta):
    """Guard: normal mode raises on records under 20 kb, which this project has."""
    if not meta:
        raise CallerMisconfigured(
            "trained/normal mode is not permitted: it raises "
            "'sequence must be at least 20000 characters' on fragment records "
            "and small plasmids, which are routine input here. Use anon mode.")


def call(seq, min_gene=MIN_GENE, closed=CLOSED):
    """Call CDS on one nucleotide sequence.

    Returns a list of dicts, 1-based inclusive genomic coordinates INCLUDING the
    stop codon (GenBank CDS convention):

        {start, end, strand, aa_len, stop_coord, start_codon, aa}

    stop_coord is the genome coordinate of the last base of the stop codon. It is
    the stable identity of a called gene: two callers that disagree about the
    start still agree about the stop.
    """
    import pyrodigal
    assert_anon_only(ANON)
    s = str(seq).upper()
    gf = pyrodigal.GeneFinder(meta=ANON, closed=closed, min_gene=min_gene)
    out = []
    for g in gf.find_genes(s):
        st = 1 if g.strand == 1 else -1
        a, b = g.begin, g.end
        sub = s[a - 1:b]
        if st == -1:
            sub = str(Seq(sub).reverse_complement())
        aa = str(Seq(sub).translate(table=TRANSL_TABLE))
        aa = aa[:-1] if aa.endswith("*") else aa
        out.append({
            "start": a, "end": b, "strand": st,
            "aa_len": len(aa),
            "stop_coord": b if st == 1 else a,
            "start_codon": sub[:3],
            "aa": aa,
        })
    return sorted(out, key=lambda c: c["start"])


def call_record(path):
    """Convenience: call on the first record of a GenBank/FASTA file."""
    from Bio import SeqIO
    fmt = "genbank" if path.endswith((".gb", ".gbk", ".gbff")) else "fasta"
    rec = next(SeqIO.parse(path, fmt))
    return rec, call(str(rec.seq))


if __name__ == "__main__":
    import os
    import sys
    PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    GB = os.path.join(PROJ, "exclusion_gene_project", "data", "seed", "genbank") \
        if not os.path.isdir(os.path.join(PROJ, "data", "seed", "genbank")) \
        else os.path.join(PROJ, "data", "seed", "genbank")
    # regression sentinel: the four MPF_T slot genes must all be called exactly
    EXPECT = [
        ("RP4__BN000925.1.gb",  27450, 27659,  1, 69,  "trbK"),
        ("pKM101__U09868.1.gb",  6533,  6760,  1, 75,  "eex"),
        ("R751__NC_001735.4.gb", 26340, 26567, 1, 75,  "trbK"),
        ("R388__NC_028464.1.gb",  5883,  6113, 1, 76,  "R388 candidate"),
    ]
    print("Pyrodigal anon -- regression sentinel on the four MPF_T slot genes\n")
    print("  min_gene=%d  closed=%s  table=%d\n" % (MIN_GENE, CLOSED, TRANSL_TABLE))
    bad = 0
    for fn, s, e, st, aa, label in EXPECT:
        p = os.path.join(GB, fn)
        if not os.path.exists(p):
            print("  %-24s RECORD MISSING" % label); bad += 1; continue
        rec, calls = call_record(p)
        hit = [c for c in calls
               if c["start"] == s and c["end"] == e and c["strand"] == st]
        ok = bool(hit) and hit[0]["aa_len"] == aa
        print("  %-24s %-16s %s  (%d calls in %s)"
              % (label, "%d..%d" % (s, e), "OK" if ok else "FAIL",
                 len(calls), rec.id))
        if not ok:
            bad += 1
    print("\n%s" % ("PASS" if bad == 0 else "FAIL: %d sentinel(s) broken" % bad))
    sys.exit(1 if bad else 0)
