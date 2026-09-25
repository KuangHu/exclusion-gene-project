#!/usr/bin/env python3
"""Call CDS once for the VirB4+ set and cache the proteins.

Every scan so far has re-run Pyrodigal over the same 7,436 plasmids: the step-1a
scan, the gap measurement, the six-frame recheck, and now the anchor ladder. The
calls are deterministic (orf_caller is frozen: anon, min_gene 90, closed, table
11), so re-deriving them is pure waste and, worse, a chance for two scans to
disagree about what a gene is.

Cache format, one protein per record:

    >{accession}|{index}|{start}|{end}|{strand}|{aa_len}

`index` is the position in orf_caller's output order for that accession, which is
the key the slot logic uses (adjacency is |i6 - i5|). Anything reading this cache
gets exactly the genes the pipeline's other steps saw.
"""
import argparse
import hashlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import orf_caller

STEP1A = "/global/scratch/users/kh36969/exclusion_gene/step1a"
PLSDB = "/global/scratch/users/kh36969/plsdb/sequences.fasta"
EXTRA = "/global/scratch/users/kh36969/exclusion_gene/mash/extra_controls.fa"
CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"


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
    from assertions import require_compute_node
    require_compute_node()          # A11: no heavy scans on a login node
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    a = ap.parse_args()
    orf_caller.assert_anon_only(orf_caller.ANON)
    from Bio import SeqIO

    _h = lambda x: int(hashlib.md5(x.encode()).hexdigest(), 16)
    want = virb4_positive()
    if a.nshards > 1:
        want = {x for x in want if _h(x) % a.nshards == a.shard}
    sys.stderr.write("accessions in this shard: %d\n" % len(want))

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
                if rec.id not in want:
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
                if n_acc % 100 == 0:
                    sys.stderr.write("  %d accessions\n" % n_acc)
    print("shard %s: %d accessions, %d CDS" % (sfx, n_acc, n_cds))


if __name__ == "__main__":
    main()
