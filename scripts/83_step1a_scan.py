#!/usr/bin/env python3
"""Pipeline v2 step 0 + 1a, fused. One shard of PLSDB.

Per plasmid: call CDS with Pyrodigal anon, scan the proteins against the anchor
set, and emit one row per plasmid plus one row per anchor hit. Protein FASTAs are
never written -- 72k of them would be pointless I/O.

The entry criterion is applied at analysis time, not here: this records every
anchor hit so that the VirB4-negative remainder stays inspectable. Filtering here
would define away the denominator, which is the whole reason VirB4 rather than
VirB6 is the entry criterion (v2 §3.1).
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import orf_caller

HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
PLSDB = "/global/scratch/users/kh36969/plsdb/sequences.fasta"
EXTRA = "/global/scratch/users/kh36969/exclusion_gene/mash/extra_controls.fa"

# (label, model file, threshold mode, min_aa floor) -- from data/anchors/anchor_set.tsv
# Exclusion-gene families are deliberately absent (v2 1.3 red line).
ANCHORS = [
    ("VirB4", "anchor_CagE_TrbE_VirB", "ga",     0),    # PF03135 -- ENTRY CRITERION
    ("VirB6", "TrbL",                  "ga",     0),    # PF04610 -- slot boundary
    ("VirB5", "anchor_T4SS",           "e1e-5", 150),   # PF07996 -- empty-slot test
    ("VirB8", "anchor_VirB8",          "ga",     0),    # PF04335 -- annotation
    ("VirB9", "anchor_CagX",           "ga",     0),    # PF03524 -- annotation
]


def load_models():
    import pyhmmer
    out = []
    for label, fn, mode, minaa in ANCHORS:
        path = os.path.join(HMM, fn + ".hmm")
        if not os.path.exists(path):
            sys.stderr.write("missing model: %s\n" % path)
            continue
        with pyhmmer.plan7.HMMFile(path) as fh:
            for hmm in fh:
                out.append((label, hmm, mode, minaa))
    return out


def shard_records(shard, nshards):
    """Stream only this shard's records. Round-robin so shards are size-balanced."""
    from Bio import SeqIO
    i = 0
    for path in (PLSDB, EXTRA):
        if not os.path.exists(path):
            continue
        for rec in SeqIO.parse(path, "fasta"):
            if i % nshards == shard:
                yield rec
            i += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--nshards", type=int, required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    models = load_models()
    sys.stderr.write("anchors loaded: %s\n"
                     % ", ".join("%s(%s)" % (l, m or "-") for l, _, m, _ in models))

    plas_p = os.path.join(a.out, "plasmids_%03d.tsv" % a.shard)
    hits_p = os.path.join(a.out, "anchor_hits_%03d.tsv" % a.shard)
    cds_p = os.path.join(a.out, "cds_%03d.tsv" % a.shard)
    fp = open(plas_p, "w", newline="")
    fh_ = open(hits_p, "w", newline="")
    fc = open(cds_p, "w", newline="")
    # lineterminator="\n" is NOT optional: csv.writer defaults to \r\n, which
    # puts a stray \r in the last field of every row. The first smoke test read
    # 146/146 plasmids as "has an anchor" because the empty anchor field held \r.
    kw = {"delimiter": "\t", "lineterminator": "\n"}
    wp = csv.writer(fp, **kw)
    wh = csv.writer(fh_, **kw)
    wc = csv.writer(fc, **kw)
    wp.writerow(["accession", "length_bp", "n_cds", "n_anchor_hits", "anchors_found"])
    wh.writerow(["accession", "anchor", "cds_index", "start", "end", "strand",
                 "aa_len", "bitscore", "evalue"])
    wc.writerow(["accession", "cds_index", "start", "end", "strand", "aa_len",
                 "start_codon"])

    n = 0
    for rec in shard_records(a.shard, a.nshards):
        seq = str(rec.seq).upper()
        if len(seq) < 1000:
            continue
        try:
            calls = orf_caller.call(seq)
        except Exception as exc:
            sys.stderr.write("CALLER FAILED on %s: %s\n" % (rec.id, exc))
            continue
        if not calls:
            wp.writerow([rec.id, len(seq), 0, 0, ""])
            continue
        digs = []
        for i, c in enumerate(calls):
            try:
                digs.append(pyhmmer.easel.TextSequence(
                    name=str(i).encode(), sequence=c["aa"]).digitize(alpha))
            except Exception:
                pass
        found = {}
        for label, hmm, mode, minaa in models:
            kw = {"bit_cutoffs": "gathering"} if mode == "ga" else {"E": 1e-5}
            try:
                for top in pyhmmer.hmmsearch([hmm], digs, cpus=1, **kw):
                    for h in top:
                        nm = h.name
                        i = int(nm.decode() if isinstance(nm, bytes) else nm)
                        c = calls[i]
                        if minaa and c["aa_len"] <= minaa:
                            continue
                        prev = found.get((label, i))
                        if prev is None or h.score > prev[0]:
                            found[(label, i)] = (h.score, h.evalue)
            except Exception as exc:
                # A4: never let a scan failure look like "no hits". A broken scan
                # and an empty scan are indistinguishable downstream, and that is
                # how the first smoke test produced a plausible 0/146.
                raise RuntimeError("anchor %s failed on %s: %s -- refusing to "
                                   "emit a result that would read as 'no hits'"
                                   % (label, rec.id, exc))
        for (label, i), (sc, ev) in sorted(found.items()):
            c = calls[i]
            wh.writerow([rec.id, label, i, c["start"], c["end"], c["strand"],
                         c["aa_len"], round(sc, 1), "%.2g" % ev])
        labs = sorted({l for l, _ in found})
        wp.writerow([rec.id, len(seq), len(calls), len(found), ";".join(labs)])
        # CDS table only for plasmids that hit the entry criterion, to bound size
        if "VirB4" in labs:
            for i, c in enumerate(calls):
                wc.writerow([rec.id, i, c["start"], c["end"], c["strand"],
                             c["aa_len"], c["start_codon"]])
        n += 1
        if n % 200 == 0:
            sys.stderr.write("  %s: %d plasmids\n" % (a.shard, n))
            fp.flush(); fh_.flush(); fc.flush()
    for f in (fp, fh_, fc):
        f.close()
    sys.stderr.write("shard %d done: %d plasmids\n" % (a.shard, n))


if __name__ == "__main__":
    main()
