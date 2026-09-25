#!/usr/bin/env python3
"""MPF_F architecture: descriptive survey first, discriminator later.

The MPF_T template does NOT transfer as-is. There the sequence was:

    an anomaly appeared (554 plasmids with VirB6 before VirB5)
      -> that DEFINED the architecture classes
      -> ground truth = the observable VirB5/VirB6 order where both are detected
      -> a discriminator that avoids VirB5 (sign of VirB6-VirB4) was built and
         validated at 98.1% against that ground truth

For MPF_F we do not yet know what the classes ARE. R27's Tra1/Tra2 split is a
hypothesis about what varies, not a measured class. Building a discriminator now
would mean assuming the answer, and defining "ground truth" for it would be
circular.

So: measure modal orders across the TraC_F_IV-admitted set and let the classes
declare themselves. Two properties are recorded because MPF_F may vary in a way
MPF_T did not -- R27 splits its transfer genes into two genomic REGIONS, which is
a different kind of rearrangement than one component relocating:

    gene ORDER      the token sequence of anchors in transcription frame
    genomic LAYOUT  contiguous vs split, from the anchor span and internal gaps

Frame: the entry anchor's strand (TraC_F_IV), which is present by construction.
An earlier MPF_T survey used VirB4's strand while the class labels came from the
VirB5/VirB6 strand, and every order string on a discordant plasmid came out
reversed. Anchors on the opposite strand are reported, not silently dropped.
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
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
ENTRY = ("TraC_F_IV", "mpff_TraC_F_IV")
SPLIT_GAP = 20000          # bp between anchor blocks that counts as a split layout


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    lines = [l for l in open(os.path.join(PROJ, "data", "anchors",
                                          "anchor_set_MPF_F.tsv"))
             if not l.startswith("#") and l.strip()]
    fams = [(r["pfam_acc"], r["pfam_name"]) for r in csv.DictReader(lines, delimiter="\t")]

    seqs, owner, idx, st, en, strand, aas = [], [], [], [], [], [], []
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
                st.append(int(f[2])); en.append(int(f[3])); strand.append(int(f[4]))
                aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    print("cache: %d proteins over %d accessions" % (len(seqs), len(set(owner))), flush=True)

    hit = {}
    for acc, name in fams:
        p = None
        for c in ("mpff_%s.hmm" % name, "anchor_%s.hmm" % name, "%s.hmm" % name):
            if os.path.exists(os.path.join(HMM, c)):
                p = os.path.join(HMM, c); break
        if not p:
            continue
        with pyhmmer.plan7.HMMFile(p) as fh:
            model = next(iter(fh))
        best = {}
        for top in pyhmmer.hmmsearch([model], seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                if owner[i] not in best or h.score > best[owner[i]][1]:
                    best[owner[i]] = (i, h.score)
        hit[name] = best
        print("  %-18s %6d plasmids" % (name, len(best)), flush=True)

    admitted = set(hit.get(ENTRY[0], {}))
    print("\nadmitted by the entry criterion %s: %d plasmids" % (ENTRY[0], len(admitted)))

    orders = collections.Counter()
    layouts = collections.Counter()
    offstrand = collections.Counter()
    rows = []
    for acc in admitted:
        e = hit[ENTRY[0]][acc][0]
        s0 = strand[e]
        t = lambda i: idx[i] * (1 if s0 == 1 else -1)
        present = {n: hit[n][acc][0] for n in hit if acc in hit[n]}
        same = {n: i for n, i in present.items() if strand[i] == s0}
        offstrand[len(present) - len(same)] += 1
        if len(same) < 4:
            continue
        seq_order = " ".join(n for n, _ in sorted(same.items(), key=lambda kv: t(kv[1])))
        orders[seq_order] += 1
        pos = sorted(st[i] for i in present.values())
        gaps = [b - a2 for a2, b in zip(pos, pos[1:])]
        big = [g for g in gaps if g >= SPLIT_GAP]
        layout = "contiguous" if not big else "split_%d" % (len(big) + 1)
        layouts[layout] += 1
        rows.append({"accession": acc, "n_anchors": len(present),
                     "n_same_strand": len(same), "layout": layout,
                     "anchor_span_bp": pos[-1] - pos[0],
                     "max_internal_gap_bp": max(gaps) if gaps else 0,
                     "order": seq_order})

    print("\n=== GENOMIC LAYOUT (anchor blocks separated by >=%d bp) ===" % SPLIT_GAP)
    for k, c in layouts.most_common():
        print("  %-14s %6d (%.1f%%)" % (k, c, 100.0 * c / max(1, sum(layouts.values()))))
    sp = sorted(r["anchor_span_bp"] for r in rows)
    print("  anchor span bp: median %d  Q1 %d  Q3 %d  max %d"
          % (statistics.median(sp), sp[len(sp)//4], sp[3*len(sp)//4], sp[-1]))

    print("\n=== MODAL GENE ORDERS (>=4 anchors on the entry anchor's strand) ===")
    tot = sum(orders.values())
    for o, c in orders.most_common(12):
        print("  %5d (%4.1f%%)  %s" % (c, 100.0 * c / tot, o))
    print("  distinct order strings: %d over %d plasmids" % (len(orders), tot))

    print("\n=== anchors on the OPPOSITE strand from the entry anchor ===")
    for k, c in sorted(offstrand.items())[:8]:
        print("  %d off-strand anchors: %5d plasmids" % (k, c))
    print("  (reported, not dropped -- a discordant strand is what reversed every")
    print("   order string in the first MPF_T survey)")

    dest = os.path.join(PROJ, "data", "anchors", "mpff_architecture_survey.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))


if __name__ == "__main__":
    main()
