#!/usr/bin/env python3
"""Is VirB5's threshold database-size dependent?

All GA anchors reproduced EXACTLY between a 1.1M-protein cache and a 7.7M-protein
cache (VirB1 5984, VirB2 5884, VirB3 6769). VirB5 alone moved, 5424 -> 5341.
VirB5 is also the only anchor whose threshold is an E-value rather than GA.

hmmsearch reports E = P * N. Enlarging the database 7x inflates every E-value 7x,
so a fixed E<=1e-5 cutoff admits FEWER proteins in the larger database -- with no
change to any sequence.

PREDICTION (stated before the run): scanning PF07996 at E<=1e-5 AND aa>150 over
only the 7,436 catalogue accessions, drawn out of cds_cache_full, reproduces 5424
exactly. If it does, the threshold is db-dependent and must be restated as a
bitscore, which is db-independent.

The script also reports the bitscore actually corresponding to E=1e-5 in each
database, and the plasmid counts at fixed bitscores, so a replacement threshold
can be chosen from measurement rather than convention.
"""
import os, sys, collections
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
FULL = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
SUB  = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
HMM  = "/global/scratch/users/kh36969/exclusion_gene/hmm/anchor_T4SS.hmm"

def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    cat = set()
    for i, l in enumerate(open(os.path.join(PROJ, "data/release/v1/catalogue_MPF_T_v1.tsv"))):
        if i: cat.add(l.split("\t")[0])

    def load(d, restrict=None):
        seqs, owner, aas = [], [], []
        nm = None
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".faa"): continue
            for line in open(os.path.join(d, fn)):
                if line[0] == ">": nm = line[1:].rstrip("\n")
                else:
                    acc = nm.split("|")[0]
                    if restrict and acc not in restrict: continue
                    owner.append(acc); aas.append(line.rstrip("\n"))
                    seqs.append(pyhmmer.easel.TextSequence(
                        name=str(len(seqs)).encode(), sequence=line.rstrip("\n")).digitize(alpha))
        return seqs, owner, aas

    with pyhmmer.plan7.HMMFile(HMM) as fh:
        model = next(iter(fh))

    def run(seqs, owner, aas, label):
        # collect ALL hits with a permissive cutoff, then threshold locally
        hits = []
        for top in pyhmmer.hmmsearch([model], seqs, cpus=16, E=10.0):
            for h in top:
                i = int(dec(h.name))
                if len(aas[i]) <= 150: continue
                hits.append((h.evalue, h.score, owner[i]))
        e5 = {o for e, s, o in hits if e <= 1e-5}
        print("\n--- %s : %d proteins searched ---" % (label, len(seqs)))
        print("  E<=1e-5 AND aa>150 : %d plasmids" % len(e5))
        near = sorted((abs(e - 1e-5), s) for e, s, o in hits)[:40]
        if near:
            bs = sorted(s for _, s in near)
            print("  bitscore at the E=1e-5 boundary: median %.1f  range %.1f-%.1f"
                  % (bs[len(bs)//2], bs[0], bs[-1]))
        for cut in (18, 20, 22, 24, 24.7, 26, 28):
            n = len({o for e, s, o in hits if s >= cut})
            print("    bitscore >= %-5s : %d plasmids" % (cut, n))
        return e5

    s, o, a = load(FULL, cat)
    got = run(s, o, a, "cds_cache_full RESTRICTED to the 7,436 catalogue accessions")
    print("\n  PREDICTION was 5424 (the figure in anchor_set.tsv / release v1)")
    print("  OBSERVED   %d   -> %s" % (len(got), "CONFIRMED" if len(got) == 5424 else "NOT confirmed"))
    del s, o, a
    s, o, a = load(SUB)
    run(s, o, a, "cds_cache (the original 1.1M cache, unrestricted)")

if __name__ == "__main__":
    main()
