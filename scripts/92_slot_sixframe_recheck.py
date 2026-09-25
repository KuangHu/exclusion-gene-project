#!/usr/bin/env python3
"""Is the empty slot really empty, or did Pyrodigal just not call it?

The VirB7 retraction's payload was not "PF20898 exists". It was: pKM101 traN is a
real, 48 aa, GA-scoring gene, and Pyrodigal -p anon does not call it. Every number
downstream of "Pyrodigal called nothing there" is therefore suspect in one
direction only -- absence -- and under v2 6.1 the occupancy figure decides whether
the whole presence/absence direction survives.

Two earlier versions of this script were wrong, both caught by controls:

  1. Searching the whole inter-anchor window returned a "novel" 58 aa ORF in 11 of
     13 occupied slots. Those windows are 250-280 bp and already hold an 85 aa
     gene: the extras were frame shadows of the called gene, every time.

  2. Restricting the search to non-coding space killed the shadows, but pKM101
     traN OVERLAPS the next called CDS by 11 bp -- that overlap is precisely why
     Pyrodigal misses it. The fix would have excluded the one gene known to be
     missed, which is the failure this script exists to detect.

So the filter is an overlap FRACTION, not a geometric exclusion: an ORF counts as
novel if its stop is not a Pyrodigal stop AND at most half of it lies inside a
called CDS. traN is 92% outside and passes; a frame shadow is ~100% inside and
does not. A gene wholly nested inside a called gene remains undetectable, which is
the same limit v2 8 already records for StORF-Reporter. Stated, not hidden.

Controls, both hard gates:
  POSITIVE -- pKM101 traN 7899..8045 must be recovered from U09868.1. This is a
              real gene that the caller really misses; it is the only control here
              that can fail for an interesting reason.
  BASELINE -- six-frame must re-find every Pyrodigal occupant. Stop-to-stop
              enumeration is a strict superset of any one-ORF-per-stop caller, so
              anything below 100% means the windowing is wrong. It was: the window
              ran anchor-end+1 to anchor-start-1, and occupants routinely overlap
              an anchor, putting their stop codon outside it. 249/2972 missed.

Scope note: six-frame is used ONLY inside an already-anchored window, the exact
scope docs/orf_caller_decision.md reserved for it. The canonical genome-wide
caller is unchanged.

Red line (v2 1.3): window geometry only. No exclusion-family HMM enters this file.
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
import orf_caller
import sixframe
from Bio.Seq import Seq

HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
PLSDB = "/global/scratch/users/kh36969/plsdb/sequences.fasta"
EXTRA = "/global/scratch/users/kh36969/exclusion_gene/mash/extra_controls.fa"

VB5 = ("anchor_T4SS", 1e-5, 150)
VB6 = ("TrbL", None, 0)
MIN_AA = 40                       # pKM101 traN is 48 aa; 40 clears it with margin
MIN_BP = 3 * MIN_AA + 3           # 123 bp of ORF
MAX_CDS_OVERLAP = 0.50            # above this, an ORF is a shadow of a called gene

FORBIDDEN = {"TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"}


def _aa(genome, o):
    sub = Seq(genome[o["start"] - 1:o["end"]])
    if o["strand"] == -1:
        sub = sub.reverse_complement()
    return str(sub.translate(table=11)).rstrip("*")


def _cds_overlap_bp(o, calls):
    """bp of ORF o that lie inside any Pyrodigal CDS."""
    tot = 0
    for c in calls:
        lo, hi = max(o["start"], c["start"]), min(o["end"], c["end"])
        if hi >= lo:
            tot += hi - lo + 1
    return tot


def scan(genome, lo_bp, hi_bp, calls, py_stops):
    """Novel ORFs in [lo_bp, hi_bp]: longest start per stop, shadows removed."""
    if hi_bp - lo_bp + 1 < MIN_BP:
        return []
    out = []
    for o in sixframe.call(genome[lo_bp - 1:hi_bp], min_aa=MIN_AA,
                           mode="first_start"):     # 5'-most start = longest per stop
        o["start"] += lo_bp - 1
        o["end"] += lo_bp - 1
        o["stop_coord"] = o["end"] if o["strand"] == 1 else o["start"]
        if (o["strand"], o["stop_coord"]) in py_stops:
            continue
        span = o["end"] - o["start"] + 1
        if _cds_overlap_bp(o, calls) > MAX_CDS_OVERLAP * span:
            continue
        out.append(o)
    return sorted(out, key=lambda x: -x["aa_len"])


def positive_control():
    """pKM101 traN: a real gene Pyrodigal misses. Must be recovered."""
    from Bio import SeqIO
    p = os.path.join(PROJ, "data", "seed", "genbank", "pKM101__U09868.1.gb")
    rec = next(SeqIO.parse(p, "genbank"))
    seq = str(rec.seq).upper()
    calls = orf_caller.call(seq)
    py_stops = {(c["strand"], c["stop_coord"]) for c in calls}
    want = (1, 8045)                                  # traN 7899..8045 (+)
    if any(c["start"] == 7899 and c["end"] == 8045 for c in calls):
        raise SystemExit("positive control invalid: Pyrodigal now calls traN "
                         "exactly, so it no longer tests a real caller miss")
    got = scan(seq, 7500, 8400, calls, py_stops)
    ok = any((o["strand"], o["stop_coord"]) == want for o in got)
    print("POSITIVE CONTROL  pKM101 traN 7899..8045 (48 aa, missed by Pyrodigal): "
          "%s" % ("RECOVERED" if ok else "*** NOT RECOVERED ***"))
    if not ok:
        raise SystemExit("positive control failed -- the method cannot find the "
                         "one gene it is meant to find. Nothing else may be read.")
    return True


def load_targets(src):
    out = {}
    for fn in sorted(os.listdir(src)):
        if not fn.startswith("virb5_virb6_gap") or not fn.endswith(".tsv"):
            continue
        for r in csv.DictReader(open(os.path.join(src, fn)), delimiter="\t"):
            if r["verdict"] in ("EMPTY", "OCCUPIED"):
                out[r["accession"]] = r["verdict"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(PROJ, "data", "positional"))
    ap.add_argument("--tag", default="v2")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    a = ap.parse_args()
    for h in (VB5[0], VB6[0]):
        if h in FORBIDDEN:
            raise SystemExit("RED LINE: exclusion-family HMM %r as an anchor" % h)
    orf_caller.assert_anon_only(orf_caller.ANON)
    positive_control()

    import hashlib
    import pyhmmer
    from Bio import SeqIO
    alpha = pyhmmer.easel.Alphabet.amino()

    def load(fn):
        with pyhmmer.plan7.HMMFile(os.path.join(HMM, fn + ".hmm")) as fh:
            return next(iter(fh))
    m5, m6 = load(VB5[0]), load(VB6[0])

    want = load_targets(a.src)
    _h = lambda x: int(hashlib.md5(x.encode()).hexdigest(), 16)
    if a.nshards > 1:
        want = {k: v for k, v in want.items() if _h(k) % a.nshards == a.shard}
    sys.stderr.write("targets: %d\n" % len(want))

    rows, stat = [], collections.Counter()
    for path in (PLSDB, EXTRA):
        if not os.path.exists(path):
            continue
        for rec in SeqIO.parse(path, "fasta"):
            if rec.id not in want:
                continue
            seq = str(rec.seq).upper()
            calls = orf_caller.call(seq)
            if not calls:
                stat["no_calls"] += 1
                continue
            digs = [pyhmmer.easel.TextSequence(name=str(i).encode(),
                                               sequence=c["aa"]).digitize(alpha)
                    for i, c in enumerate(calls)]

            def hits(model, ev, minaa):
                out = {}
                kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
                for top in pyhmmer.hmmsearch([model], digs, cpus=1, **kw):
                    for h in top:
                        nm = h.name
                        i = int(nm.decode() if isinstance(nm, bytes) else nm)
                        if minaa and calls[i]["aa_len"] <= minaa:
                            continue
                        if i not in out or h.score > out[i]:
                            out[i] = h.score
                return out
            h5, h6 = hits(m5, VB5[1], VB5[2]), hits(m6, VB6[1], VB6[2])
            best = None
            for i5 in h5:
                for i6 in h6:
                    if calls[i5]["strand"] != calls[i6]["strand"]:
                        continue
                    d = abs(i6 - i5)
                    if d == 0 or d > 6:
                        continue
                    if best is None or d < best[0]:
                        best = (d, i5, i6)
            if best is None:
                stat["anchor_pair_not_reproduced"] += 1
                continue
            d, i5, i6 = best
            lo, hi = (i5, i6) if i5 < i6 else (i6, i5)
            occ = [calls[j] for j in range(lo + 1, hi)]
            inter_bp = calls[hi]["start"] - calls[lo]["end"] - 1

            # Discovery window spans the anchors themselves. Occupants overlapping
            # an anchor are then fully contained (that was the 249-miss bug), and
            # the anchors and their shadows are removed by the overlap filter, not
            # by cutting them out of the search space.
            w_lo, w_hi = calls[lo]["start"], calls[hi]["end"]
            py_stops = {(c["strand"], c["stop_coord"]) for c in calls}

            # BASELINE: recover every Pyrodigal occupant, matched on stop.
            sf_all = sixframe.call(seq[w_lo - 1:w_hi], min_aa=MIN_AA,
                                   mode="first_start")
            sf_stops = set()
            for o in sf_all:
                s, e = o["start"] + w_lo - 1, o["end"] + w_lo - 1
                sf_stops.add((o["strand"], e if o["strand"] == 1 else s))
            recovered = missed = 0
            for c in occ:
                if c["aa_len"] < MIN_AA:
                    continue
                if (c["strand"], c["stop_coord"]) in sf_stops:
                    recovered += 1
                else:
                    missed += 1
                    stat["BASELINE_MISS"] += 1
            stat["baseline_recovered"] += recovered

            nov = scan(seq, w_lo, w_hi, calls, py_stops)

            # Background: the same procedure over the whole plasmid, so the slot
            # rate is compared against this plasmid's own ORF density.
            bg = scan(seq, 1, len(seq), calls, py_stops)
            bg_rate = 1000.0 * len(bg) / len(seq)
            expected = bg_rate * (w_hi - w_lo + 1) / 1000.0

            verdict = want[rec.id]
            stat[verdict] += 1
            if verdict == "EMPTY":
                if inter_bp < MIN_BP:
                    stat["EMPTY_gap_too_short"] += 1
                elif nov:
                    stat["EMPTY_now_occupied"] += 1
                else:
                    stat["EMPTY_confirmed"] += 1
            elif nov:
                stat["OCCUPIED_extra_orf"] += 1

            rows.append({
                "accession": rec.id, "verdict_pyrodigal": verdict,
                "intergenic_bp": inter_bp,
                "window": "%d..%d" % (w_lo, w_hi), "window_bp": w_hi - w_lo + 1,
                "n_pyrodigal_in_window": len(occ),
                "pyrodigal_lengths": ";".join(str(c["aa_len"]) for c in occ),
                "baseline_recovered": recovered, "baseline_missed": missed,
                "bg_orf_per_kb": round(bg_rate, 3),
                "expected_by_chance": round(expected, 3),
                "n_sixframe_novel": len(nov),
                "sixframe_novel_lengths": ";".join(str(o["aa_len"]) for o in nov),
                "sixframe_novel_coords": ";".join(
                    "%d..%d%s" % (o["start"], o["end"],
                                  "+" if o["strand"] == 1 else "-") for o in nov),
                "sixframe_novel_starts": ";".join(o["start_codon"] for o in nov),
                "sixframe_novel_seqs": ";".join(_aa(seq, o) for o in nov),
            })

    out = os.path.join(PROJ, "data", "positional")
    os.makedirs(out, exist_ok=True)
    sfx = "" if a.nshards == 1 else "_%03d" % a.shard
    dest = os.path.join(out, "slot_sixframe_%s%s.tsv" % (a.tag, sfx))
    cols = ["accession", "verdict_pyrodigal", "intergenic_bp", "window",
            "window_bp", "n_pyrodigal_in_window", "pyrodigal_lengths",
            "baseline_recovered", "baseline_missed", "bg_orf_per_kb",
            "expected_by_chance", "n_sixframe_novel", "sixframe_novel_lengths",
            "sixframe_novel_coords", "sixframe_novel_starts", "sixframe_novel_seqs"]
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    for k, v in sorted(stat.items()):
        print("%-32s %d" % (k, v))
    if stat["BASELINE_MISS"]:
        print("*** BASELINE FAILED (%d misses) -- do not read any other number ***"
              % stat["BASELINE_MISS"])
    print("wrote %s (%d rows)" % (dest, len(rows)))


if __name__ == "__main__":
    main()
