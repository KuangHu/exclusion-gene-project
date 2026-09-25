#!/usr/bin/env python3
"""Six-frame stop-to-stop enumeration with FULL multi-start retention.

This exists because orfipy, StORF-Reporter and every Prodigal variant emit at
most ONE ORF per stop codon. That is a model-expressiveness limit, not a
precision setting: R64 excA and excB share a stop codon at AP005147.1:76983 and
differ only in which in-frame start is used (219 nt apart). No one-per-stop
caller can emit both, at any threshold.

Algorithm, per frame, per strand:
  1. split the frame into stop-to-stop intervals
  2. within each interval, enumerate EVERY in-frame candidate start
  3. emit one record per (stop, start) pair above the length floor

Two modes:
  first_start   one ORF per stop, using the 5'-most start (comparable to orfipy)
  all_starts    every in-frame start per stop (the mode that can express excA+excB)
"""
from Bio.Seq import Seq

STOPS = {"TAA", "TAG", "TGA"}
DEFAULT_STARTS = ("ATG", "GTG", "TTG")


def _frame_orfs(seq, frame, starts, min_aa, mode):
    """Yield (start0, end0, aa_len) in seq coordinates for one forward frame."""
    n = len(seq)
    codon_start = frame
    interval_begin = frame          # first base of the current stop-to-stop interval
    out = []
    i = codon_start
    while i + 3 <= n:
        cod = seq[i:i + 3]
        if cod in STOPS:
            # interval is [interval_begin, i+3) with the stop included at the end
            cands = []
            j = interval_begin
            while j + 3 <= i:
                if seq[j:j + 3] in starts:
                    aa = (i - j) // 3          # coding aa, stop excluded
                    if aa >= min_aa:
                        cands.append((j, i + 3, aa))
                    if mode == "first_start":
                        break
                j += 3
            out.extend(cands)
            interval_begin = i + 3
        i += 3
    return out


def call(seq_str, starts=DEFAULT_STARTS, min_aa=60, mode="all_starts"):
    """Return list of dicts with 1-based inclusive genomic coords and strand.

    `stop_coord` is the last base of the stop codon on the genome, which is the
    key the 3'-recall test matches on.
    """
    s = seq_str.upper()
    n = len(s)
    rc = str(Seq(s).reverse_complement())
    res = []
    for frame in (0, 1, 2):
        for st, work in ((1, s), (-1, rc)):
            for a, b, aa in _frame_orfs(work, frame, set(starts), min_aa, mode):
                if st == 1:
                    g_start, g_end = a + 1, b            # 1-based inclusive
                else:
                    g_start, g_end = n - b + 1, n - a
                res.append({
                    "start": g_start, "end": g_end, "strand": st,
                    "aa_len": aa,
                    "stop_coord": g_end if st == 1 else g_start,
                    "start_codon": work[a:a + 3],
                    "frame": frame,
                })
    return res


if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from Bio import SeqIO
    import testset

    # sanity: can all_starts express R64 excA AND excB from the same stop?
    p = testset.rec_path("AP005147.1")
    rec = next(SeqIO.parse(p, "genbank"))
    sub = str(rec.seq)[76000:78200]                # 1-based 76001..78200
    off = 76000
    for mode in ("first_start", "all_starts"):
        orfs = call(sub, min_aa=60, mode=mode)
        # excA 76983..77645(-)  excB 76983..77426(-)  -> stop at 76983
        want_stop = 76983 - off
        same = [o for o in orfs if o["strand"] == -1 and o["stop_coord"] == want_stop]
        print("mode=%-12s ORFs=%-5d sharing the excA/excB stop: %d"
              % (mode, len(orfs), len(same)))
        for o in sorted(same, key=lambda x: -x["aa_len"]):
            print("    aa=%-4d genomic %d..%d (%s) start=%s"
                  % (o["aa_len"], o["start"] + off, o["end"] + off,
                     "+" if o["strand"] == 1 else "-", o["start_codon"]))
