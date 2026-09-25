#!/usr/bin/env python3
"""MPF_I architecture: descriptive first, no presupposed number of classes.

Entry is T4SS_I_traU (3,732 plasmids). Two dimensions measured together, because
MPF_F showed they are not the same thing -- there, layout is four-valued and 35%
split, which no gene-order string can express:

    gene ORDER   anchor tokens in transcription frame, canonicalised for circular
                 rotation and for multi-domain proteins contributing two tokens
    LAYOUT       contiguous vs split, from anchor-block spacing (>=20 kb)

Reference points, not expectations: MPF_T resolved to three architectures
(canonical / IncI2 / pEC4115, the third found only by chasing an anomaly); MPF_F
to layout four-valued x 298 canonical orders with one dominant at 60.5%.

Frame is the entry anchor's strand. Off-strand anchors are REPORTED, not dropped:
the first MPF_T survey took its frame from VirB4 while the class labels came from
the VirB5/VirB6 strand, and every order string on a discordant plasmid came out
reversed.

I_traR is excluded from the anchor set: 179 cross-class admissions (4.8%), 166 of
them MPF_T -- the only genuine contamination among the seventeen.
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
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
ENTRY = "T4SS_I_traU"
EXCLUDE = {"T4SS_I_traR"}
SPLIT_GAP = 20000


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner, idx, st, strand = [], [], [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                owner.append(f[0]); idx.append(int(f[1]))
                st.append(int(f[2])); strand.append(int(f[4]))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    print("cache: %d proteins over %d accessions" % (len(seqs), len(set(owner))), flush=True)

    profiles = sorted(f[:-4] for f in os.listdir(CJ)
                      if f.startswith("T4SS_I_") and f.endswith(".hmm")
                      and f[:-4] not in EXCLUDE)
    hit = {}
    for p in profiles:
        with pyhmmer.plan7.HMMFile(os.path.join(CJ, p + ".hmm")) as fh:
            m = next(iter(fh))
        best = {}
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                if owner[i] not in best or h.score > best[owner[i]][1]:
                    best[owner[i]] = (i, h.score)
        hit[p] = best
        print("  %-16s %6d" % (p.replace("T4SS_", ""), len(best)), flush=True)

    admitted = sorted(hit[ENTRY])
    print("\nadmitted by %s: %d" % (ENTRY, len(admitted)))
    short = lambda p: p.replace("T4SS_I_", "")

    def canon(toks):
        if not toks:
            return ""
        return min(" ".join(toks[i:] + toks[:i]) for i in range(len(toks)))

    orders, layouts, offs, rows = collections.Counter(), collections.Counter(), collections.Counter(), []
    for acc in admitted:
        e = hit[ENTRY][acc][0]
        s0 = strand[e]
        t = lambda i: idx[i] * (1 if s0 == 1 else -1)
        present = {short(p): hit[p][acc][0] for p in profiles if acc in hit[p]}
        same = {n: i for n, i in present.items() if strand[i] == s0}
        offs[len(present) - len(same)] += 1
        if len(same) < 4:
            continue
        toks = [n for n, _ in sorted(same.items(), key=lambda kv: t(kv[1]))]
        orders[canon(toks)] += 1
        pos = sorted(st[i] for i in present.values())
        gaps = [b - a2 for a2, b in zip(pos, pos[1:])]
        big = [g for g in gaps if g >= SPLIT_GAP]
        layout = "contiguous" if not big else "split_%d" % (len(big) + 1)
        layouts[layout] += 1
        rows.append({"accession": acc, "n_anchors": len(present),
                     "n_same_strand": len(same), "layout": layout,
                     "anchor_span_bp": pos[-1] - pos[0],
                     "max_internal_gap_bp": max(gaps) if gaps else 0,
                     "order_canonical": canon(toks)})

    print("\n=== LAYOUT (blocks separated by >=%d bp) ===" % SPLIT_GAP)
    tot = max(1, sum(layouts.values()))
    for k, c in layouts.most_common():
        print("  %-14s %6d (%.1f%%)" % (k, c, 100.0 * c / tot))
    sp = sorted(r["anchor_span_bp"] for r in rows)
    print("  anchor span bp: median %d  Q1 %d  Q3 %d  max %d"
          % (statistics.median(sp), sp[len(sp)//4], sp[3*len(sp)//4], sp[-1]))

    print("\n=== CANONICAL GENE ORDERS ===")
    to = sum(orders.values())
    for o, c in orders.most_common(10):
        print("  %5d (%5.1f%%)  %s" % (c, 100.0 * c / to, o[:96]))
    print("  distinct canonical orders: %d over %d plasmids" % (len(orders), to))
    cum = 0
    for i, (o, c) in enumerate(orders.most_common(), 1):
        cum += c
        if i in (1, 2, 3, 5, 10):
            print("  top %-3d cumulative %5.1f%%" % (i, 100.0 * cum / to))

    print("\n=== layout x dominant order ===")
    top = orders.most_common(1)[0][0]
    x = collections.Counter((r["layout"], "dominant" if r["order_canonical"] == top else "other")
                            for r in rows)
    for k in sorted(x):
        print("  %-12s %-10s %5d" % (k[0], k[1], x[k]))

    print("\n=== off-strand anchors (reported, not dropped) ===")
    for k, c in sorted(offs.items())[:8]:
        print("  %d off-strand: %5d plasmids" % (k, c))

    dest = os.path.join(PROJ, "data", "anchors", "mpfi_architecture_survey.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))


if __name__ == "__main__":
    main()
