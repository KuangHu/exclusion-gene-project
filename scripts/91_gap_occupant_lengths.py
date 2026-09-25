#!/usr/bin/env python3
"""What is actually sitting in the VirB5|VirB6 gap? Length distribution.

The LOO rejection counts already contain a preliminary occupancy figure, and it
lands in the early-kill zone: for VirB5, 105 gaps held an extra gene against 6
that did not, i.e. ~95%. Under pipeline_v2 §6.1, occupancy near 100% means
presence/absence carries no signal and the criterion must change.

But that 95% is ambiguous, and the ambiguity is decisive:

    extra gene 69-76 aa  -> plausibly the exclusion gene      -> real occupancy
    extra gene ~48 aa    -> an unrecognised canonical component -> Pfam blind spot

VirB7 is a known blind spot (pKM101 traN, 48 aa, no Pfam family at all). ANY
unrecognised canonical component between two GA anchors triggers the same count
mismatch. So the question is what the length distribution looks like:

    unimodal near 70  -> genuinely occupied; occupancy ~95%; change direction
    bimodal 50 + 70   -> partly blind-spot artefact; true occupancy much lower
    unimodal near 50  -> mostly VirB7-class artefact, not exclusion at all

This is measurable now, and it outranks the step-1a three numbers because it
decides whether the whole presence/absence direction survives.

Reported for the VirB5|VirB6 gap specifically -- the slot -- not for all gaps.
"""
import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import orf_caller

HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
STEP1A = "/global/scratch/users/kh36969/exclusion_gene/step1a"
PLSDB = "/global/scratch/users/kh36969/plsdb/sequences.fasta"
EXTRA = "/global/scratch/users/kh36969/exclusion_gene/mash/extra_controls.fa"

# Only the two anchors that bound the slot. VirB5 needs its full rule
# (E<=1e-5 AND aa>150) because GA misses both IncP controls.
VB5 = ("anchor_T4SS", 1e-5, 150)
VB6 = ("TrbL", None, 0)          # GA


def virb4_positive():
    accs = set()
    for fn in sorted(os.listdir(STEP1A)):
        if not fn.startswith("plasmids_"):
            continue
        for line in open(os.path.join(STEP1A, fn)):
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[0] == "accession":
                continue
            if "VirB4" in f[4]:
                accs.add(f[0])
    return accs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    a = ap.parse_args()
    import pyhmmer
    from Bio import SeqIO
    alpha = pyhmmer.easel.Alphabet.amino()

    def load(fn):
        with pyhmmer.plan7.HMMFile(os.path.join(HMM, fn + ".hmm")) as fh:
            return next(iter(fh))
    m5, m6 = load(VB5[0]), load(VB6[0])

    want = virb4_positive()
    import hashlib
    _h = lambda x: int(hashlib.md5(x.encode()).hexdigest(), 16)
    if a.limit:
        # NOT sorted()[:n]. That took the first 800 by accession string, which on
        # PLSDB is a contiguous submission block: 694 CP + 95 AP and ZERO NZ_,
        # against a full set that is 80% NZ_. Deterministic hash sampling instead.
        keep = sorted(want, key=_h)[:a.limit]
        want = set(keep)
    if a.nshards > 1:
        want = {x for x in want if _h(x) % a.nshards == a.shard}
    sys.stderr.write("VirB4+ plasmids to examine: %d\n" % len(want))

    occ_lengths = []          # aa length of every gene in an occupied VirB5|VirB6 gap
    rej = collections.Counter()   # why a VirB4+ plasmid yielded no measurable gap
    n_gap_found = n_empty = n_occupied = 0
    multi = 0
    rows = []
    n = 0
    for path in (PLSDB, EXTRA):
        if not os.path.exists(path):
            continue
        for rec in SeqIO.parse(path, "fasta"):
            if rec.id not in want:
                continue
            calls = orf_caller.call(str(rec.seq).upper())
            if not calls:
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
            h5 = hits(m5, VB5[1], VB5[2])
            h6 = hits(m6, VB6[1], VB6[2])
            n += 1
            if not h5 and not h6:
                rej["no_VirB5_and_no_VirB6"] += 1; continue
            if not h5:
                rej["no_VirB5_only"] += 1; continue
            if not h6:
                rej["no_VirB6_only"] += 1; continue
            # the VirB5 and VirB6 that are nearest each other, same strand
            best = None
            for i5 in h5:
                for i6 in h6:
                    if calls[i5]["strand"] != calls[i6]["strand"]:
                        continue
                    d = abs(i6 - i5)
                    if d == 0 or d > 6:
                        continue
                    # NO presupposed order. pipeline_v2 4.4: name the slot by the
                    # LOCALLY OBSERVED adjacent pair, never by forcing global VirB
                    # numbering. Requiring VirB5-before-VirB6 excluded 78.3% of the
                    # "not pairable" class -- ~688 plasmids where VirB6 simply
                    # precedes VirB5. Those are skeleton rearrangements (Guglielmini
                    # 2013: virB5/virB6 sometimes follow virB10), the slot is intact,
                    # and they are ENRICHED for divergent systems -- exactly the
                    # subpopulation most likely to carry a novel family.
                    if best is None or d < best[0]:
                        best = (d, i5, i6)
            if best is None:
                wide = None
                for i5 in h5:
                    for i6 in h6:
                        if calls[i5]["strand"] != calls[i6]["strand"]:
                            continue
                        dd = abs(i6 - i5)
                        if dd and (wide is None or dd < wide[0]):
                            wide = (dd, i5, i6)
                if wide:
                    rej["gap_gt_6_genes"] += 1
                    rows.append({"accession": rec.id, "n_between": wide[0] - 1,
                                 "lengths": "", "verdict": "GAP_TOO_WIDE",
                                 "coords": "", "seqs": "", "order": ""})
                else:
                    rej["opposite_strand_only"] += 1
                continue
            d, i5, i6 = best
            fwd = calls[i5]["strand"] == 1
            canonical = (i5 < i6) if fwd else (i5 > i6)
            n_gap_found += 1
            lo, hi = (i5, i6) if i5 < i6 else (i6, i5)
            between = [calls[j] for j in range(lo + 1, hi)]
            if not between:
                n_empty += 1
                rows.append({"accession": rec.id, "n_between": 0, "lengths": "",
                             "verdict": "EMPTY", "coords": "", "seqs": "",
                             "order": "canonical" if canonical else "inverted"})
            else:
                n_occupied += 1
                if len(between) > 1:
                    multi += 1
                for c in between:
                    occ_lengths.append(c["aa_len"])
                rows.append({"accession": rec.id, "n_between": len(between),
                             "lengths": ";".join(str(c["aa_len"]) for c in between),
                             "verdict": "OCCUPIED",
                             "coords": ";".join("%d..%d" % (c["start"], c["end"])
                                                for c in between),
                             "seqs": ";".join(c["aa"] for c in between),
                             "order": "canonical" if canonical else "inverted"})
            if n_gap_found % 250 == 0:
                sys.stderr.write("  %d gaps\n" % n_gap_found)

    PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = os.path.join(PROJ, "data", "positional")
    os.makedirs(out, exist_ok=True)
    sfx = "" if a.nshards == 1 else "_%03d" % a.shard
    with open(os.path.join(out, "virb5_virb6_gap%s.tsv" % sfx), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "n_between", "lengths",
                                           "verdict", "order", "coords", "seqs"],
                           delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    print("=" * 78)
    print("VirB5 | VirB6 GAP -- what occupies it?")
    print("=" * 78)
    print("  plasmids examined            %d" % n)
    print("  with BOTH anchors, same strand, <=6 genes apart:  %d" % n_gap_found)
    print("\n  WHY THE REST WERE EXCLUDED. v2 3.1's logic applies here too: VirB5")
    print("  is the WEAK anchor, so these exclusions ENRICH for divergent systems")
    print("  and the occupancy figure is conditional on VirB5 being detectable.")
    for k, v in rej.most_common():
        print("    %-34s %6d  (%.1f%% of examined)" % (k, v, 100.0 * v / max(1, n)))
    if not n_gap_found:
        return
    print("  gap EMPTY (adjacent)         %d  (%.1f%%)"
          % (n_empty, 100.0 * n_empty / n_gap_found))
    print("  gap OCCUPIED                 %d  (%.1f%%)"
          % (n_occupied, 100.0 * n_occupied / n_gap_found))
    print("    of those, >1 gene in gap   %d" % multi)

    print("\n  LENGTH DISTRIBUTION of gap occupants (n=%d)" % len(occ_lengths))
    print("  %-14s %6s  %s" % ("aa range", "count", "bar"))
    h = collections.Counter(min(20, l // 25) for l in occ_lengths)
    for b in range(0, 21):
        c = h.get(b, 0)
        lab = "%d-%d" % (b * 25, b * 25 + 24) if b < 20 else ">=500"
        bar = "#" * max(0, int(c / max(1, max(h.values())) * 50))
        star = ""
        if b == 2:
            star = "   <- 50-74: VirB7-class blind spot (traN 48 aa)"
        if b == 3:
            star = "   <- 75-99: exclusion-gene band (69-76 aa straddles 2/3)"
        print("  %-14s %6d  %s%s" % (lab, c, bar, star))

    import statistics
    if occ_lengths:
        s = sorted(occ_lengths)
        print("\n  median %d   Q1 %d   Q3 %d   min %d   max %d"
              % (statistics.median(s), s[len(s)//4], s[3*len(s)//4], s[0], s[-1]))
        band = sum(1 for l in s if 60 <= l <= 100)
        small = sum(1 for l in s if l < 60)
        print("  in the 60-100 aa exclusion band : %d (%.1f%%)"
              % (band, 100.0 * band / len(s)))
        print("  below 60 aa (VirB7-class)       : %d (%.1f%%)"
              % (small, 100.0 * small / len(s)))
    print("\nwrote %s" % os.path.join(out, "virb5_virb6_gap%s.tsv" % sfx))


if __name__ == "__main__":
    main()
