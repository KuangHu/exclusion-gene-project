#!/usr/bin/env python3
"""CDS cache for ALL of PLSDB, not just the VirB4+ subset.

The existing cache holds the 7,436 plasmids that pass the MPF_T entry criterion
(VirB4/PF03135). That is fine for MPF_T and useless for anything else: an MPF_F
entry criterion measured inside a VirB4-selected set would be measuring VirB4.

MPF_T's entry anchor was chosen by running candidates across the whole database
and taking the most conservative. MPF_F, MPF_I and MPF_FA need the same, and
PLSDB holds 11,714 MPF_F and 3,613 MPF_I plasmids by MOB-suite typing that the
current cache never saw.

Same frozen caller (anon, min_gene 90, closed, table 11) and the same header
format, so everything downstream reads it unchanged.
"""
import argparse
import hashlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import orf_caller
from assertions import require_compute_node

PLSDB = "/global/scratch/users/kh36969/plsdb/sequences.fasta"
EXTRA = "/global/scratch/users/kh36969/exclusion_gene/mash/extra_controls.fa"
CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    a = ap.parse_args()
    orf_caller.assert_anon_only(orf_caller.ANON)
    from Bio import SeqIO
    _h = lambda x: int(hashlib.md5(x.encode()).hexdigest(), 16)
    os.makedirs(CACHE, exist_ok=True)
    sfx = "%03d" % a.shard
    n_acc = n_cds = 0
    with open(os.path.join(CACHE, "cds_%s.faa" % sfx), "w") as fa, \
         open(os.path.join(CACHE, "index_%s.tsv" % sfx), "w") as ix:
        ix.write("accession\tlength_bp\tn_cds\n")
        for path in (PLSDB, EXTRA):
            if not os.path.exists(path):
                continue
            for rec in SeqIO.parse(path, "fasta"):
                if a.nshards > 1 and _h(rec.id) % a.nshards != a.shard:
                    continue
                seq = str(rec.seq).upper()
                calls = orf_caller.call(seq)
                n_acc += 1
                n_cds += len(calls)
                ix.write("%s\t%d\t%d\n" % (rec.id, len(seq), len(calls)))
                for i, c in enumerate(calls):
                    fa.write(">%s|%d|%d|%d|%d|%d\n%s\n"
                             % (rec.id, i, c["start"], c["end"], c["strand"],
                                c["aa_len"], c["aa"]))
                if n_acc % 500 == 0:
                    sys.stderr.write("  %d accessions\n" % n_acc)
    print("shard %s: %d accessions, %d CDS" % (sfx, n_acc, n_cds))


if __name__ == "__main__":
    main()
