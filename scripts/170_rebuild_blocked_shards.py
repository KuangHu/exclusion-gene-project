#!/usr/bin/env python3
"""Rebuild the cds_cache_full shards made unreadable by the Lustre OST 106 fault.

Lustre OST 106 (lr6-OST006a) has been returning ESHUTDOWN since 2026-09-18. The
failure is per-FILE: any file whose first stripe component landed on OST 106
cannot be read at all. In cds_cache_full that is 10 of 34 `cds_*.faa` and 9 of 34
`index_*.tsv` -- roughly 18,200 of 72,556 plasmids, which blocks every scan that
needs the whole database.

`plsdb/sequences.fasta` is blocked too, but `sequences.fasta.bz2` is intact, so
the source is not lost. Sharding in script 109 is `md5(accession) % 32`, and the
ORF caller is frozen (anon, min_gene 90, closed, table 11), so a shard rebuilt
from the same records is reproducible byte for byte.

  THE CONTROL IS THE POINT. This script refuses to rebuild anything until it has
  rebuilt a shard that IS readable and byte-compared it against the original. If
  that control does not match exactly, the rebuild is not reproducing the
  pipeline's genes and every downstream coordinate would silently shift --
  gene INDEX is the slot key, so a one-gene drift moves every slot on the shard.

Output goes to a staging directory. Nothing is written into cds_cache_full by
this script; installation is a separate, explicit step once the control passes.
"""
import argparse
import hashlib
import os
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
import orf_caller
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
EXTRA = "/global/scratch/users/kh36969/exclusion_gene/mash/extra_controls.fa"
STAGE = "/global/scratch/users/kh36969/exclusion_gene/recovery/cds_cache_rebuild"
NSHARDS = 32

_h = lambda x: int(hashlib.md5(x.encode()).hexdigest(), 16)


def readable(path):
    try:
        with open(path, "rb") as fh:
            fh.read(32)
        return True
    except OSError:
        return False


def build(records, shard, nshards, outdir, sfx=None):
    """Byte-identical to script 109's writer. Do not reformat."""
    from Bio import SeqIO
    sfx = sfx if sfx is not None else "%03d" % shard
    fa_p = os.path.join(outdir, "cds_%s.faa" % sfx)
    ix_p = os.path.join(outdir, "index_%s.tsv" % sfx)
    n_acc = n_cds = 0
    with open(fa_p, "w") as fa, open(ix_p, "w") as ix:
        ix.write("accession\tlength_bp\tn_cds\n")
        for path in records:
            if not os.path.exists(path):
                raise SystemExit("input missing: %s -- refusing to write a short shard" % path)
            for rec in SeqIO.parse(path, "fasta"):
                if nshards > 1 and _h(rec.id) % nshards != shard:
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
    return fa_p, ix_p, n_acc, n_cds


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--plsdb", required=True, help="rebuilt sequences.fasta")
    ap.add_argument("--control", type=int, default=0,
                    help="a READABLE shard to reproduce and byte-compare")
    ap.add_argument("--shards", default="", help="comma list to rebuild, after the control passes")
    a = ap.parse_args()
    orf_caller.assert_anon_only(orf_caller.ANON)
    os.makedirs(STAGE, exist_ok=True)
    srcs = [a.plsdb, EXTRA]

    # ---- the control: reproduce a shard we can still read -------------------
    ctl = "%03d" % a.control
    orig = os.path.join(CACHE, "cds_%s.faa" % ctl)
    if not readable(orig):
        raise SystemExit("control shard %s is itself unreadable -- pick another" % ctl)
    print("=== CONTROL: rebuilding readable shard %s ===" % ctl, flush=True)
    cdir = os.path.join(STAGE, "control")
    os.makedirs(cdir, exist_ok=True)
    fa_p, ix_p, n_acc, n_cds = build(srcs, a.control, NSHARDS, cdir)
    print("  rebuilt %d accessions, %d CDS" % (n_acc, n_cds), flush=True)
    s_new, s_old = sha(fa_p), sha(orig)
    print("  rebuilt sha256 %s" % s_new)
    print("  original sha256 %s" % s_old)
    if s_new != s_old:
        print("\n  CONTROL FAILED -- the rebuild does not reproduce the cache.")
        print("  STOPPING. Nothing else in this run may be read: gene INDEX is the")
        print("  slot key, so any drift silently relocates every slot on the shard.")
        d = subprocess.run(["diff", fa_p, orig], capture_output=True, text=True)
        print("  first differing lines:\n%s" % "\n".join(d.stdout.splitlines()[:8]))
        return 1
    print("  CONTROL PASSED -- byte identical.\n", flush=True)

    if not a.shards:
        print("no --shards given; control only.")
        return 0

    # ---- rebuild the blocked shards ----------------------------------------
    for tok in a.shards.split(","):
        tok = tok.strip()
        if not tok:
            continue
        shard = int(tok)
        sfx = "%03d" % shard
        if shard >= 900:
            # 900/901 came from the EXTRA pass, not the 32-way PLSDB split
            print("skipping %s: EXTRA-pass shard, rebuild separately" % sfx)
            continue
        print("=== rebuilding shard %s ===" % sfx, flush=True)
        fa_p, ix_p, n_acc, n_cds = build(srcs, shard, NSHARDS, STAGE)
        print("  %s: %d accessions, %d CDS, %d bytes"
              % (sfx, n_acc, n_cds, os.path.getsize(fa_p)), flush=True)
        if n_acc == 0:
            print("  REFUSING: shard is empty. Input is wrong.")
            return 1
    print("\nstaged in %s -- install is a separate step." % STAGE)


if __name__ == "__main__":
    sys.exit(main() or 0)
