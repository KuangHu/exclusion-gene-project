#!/usr/bin/env python3
"""The MPF_I candidate container: what is actually at TraY+1 where ExcA is absent?

ExcA presence turned out to be strongly architecture-dependent -- 99.7% in the
stratum all five seeds occupy, 18.5% in split_2, 0% in order#2 and split_3+. So
~1,460 TraY-positive plasmids carry no named exclusion family. Before calling that
a candidate pool, two things have to be checked.

(1) order#2, 132 plasmids, zero ExcA, and its token string is missing traE, traT,
    traV, trbA and trbB relative to order#1 -- a CONSISTENT absence of the same
    five genes, not random loss. What sits at TraY+1? If one unnamed protein
    recurs, that is the same shape as the 84 aa IncI2 lipoprotein found earlier.
    Counted by MASH CLUSTER, not records: IncI1 is highly clonal in clinical
    isolates and 132 records could be a handful of lineages.

(2) split_2, 974 ExcA-negative plasmids. A split layout means anchor blocks
    separated by >=20 kb, so "TraY+1" may straddle a breakpoint. If TraY is the
    last gene of its block, its +1 neighbour is not part of the same cluster and
    the position is meaningless there. Measured as the bp distance from TraY to
    its +1 neighbour.

Also reports whether the order#2 plasmids retain a relaxase and a coupling protein
-- a system missing five tra genes AND a relaxase is degenerate, and its empty slot
is not evidence of anything.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import os
import re
import statistics
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
MASH = "/global/scratch/users/kh36969/exclusion_gene/mash"
LIPO = re.compile(r"[LVI][ASTVIG][GASN]C")
MOB = ["T4SS_MOBB", "T4SS_MOBC", "T4SS_MOBF", "T4SS_MOBH", "T4SS_MOBP1",
       "T4SS_MOBP2", "T4SS_MOBP3", "T4SS_MOBQ", "T4SS_MOBT", "T4SS_MOBV"]


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    csv.field_size_limit(10 ** 8)

    surv = {r["accession"]: r for r in csv.DictReader(
        open(os.path.join(PROJ, "data", "anchors", "mpfi_architecture_survey.tsv")),
        delimiter="\t")}
    oc = collections.Counter(r["order_canonical"] for r in surv.values())
    rank = {o: i + 1 for i, (o, _) in enumerate(oc.most_common())}

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
    bypos = {(owner[i], idx[i]): i for i in range(len(owner))}
    print("cache: %d proteins" % len(seqs), flush=True)

    def scan(path):
        with pyhmmer.plan7.HMMFile(path) as fh:
            m = next(iter(fh))
        out = collections.defaultdict(list)
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                out[owner[i]].append((idx[i], strand[i], h.score, i))
        return out
    tray = scan(os.path.join(CJ, "T4SS_I_traY.hmm"))
    exca = scan(os.path.join(HMM, "NF033891.hmm"))
    t4cp = scan(os.path.join(CJ, "T4SS_t4cp2.hmm"))
    relax = set()
    for m in MOB:
        for acc in scan(os.path.join(CJ, m + ".hmm")):
            relax.add(acc)
    print("TraY %d, ExcA %d, t4cp2 %d, any relaxase %d"
          % (len(tray), len(exca), len(t4cp), len(relax)), flush=True)
    best = lambda d, acc: max(d[acc], key=lambda x: x[2]) if acc in d else None

    ids = [l.strip() for l in open(os.path.join(MASH, "all_ids.txt")) if l.strip()]
    ix = {a2: i for i, a2 in enumerate(ids)}
    p = list(range(len(ids)))

    def find(x):
        while p[x] != x:
            p[x] = p[p[x]]; x = p[x]
        return x
    with gzip.open(os.path.join(MASH, "pairs_d010.tsv.gz"), "rt") as fh:
        for line in fh:
            f = line.split("\t", 3)
            if float(f[2]) > 0.007:
                continue
            u, v = ix.get(f[0]), ix.get(f[1])
            if u is None or v is None:
                continue
            ru, rv = find(u), find(v)
            if ru != rv:
                p[rv] = ru

    o2 = [acc for acc, r in surv.items()
          if r["layout"] == "contiguous" and rank.get(r["order_canonical"]) == 2]
    sp2 = [acc for acc, r in surv.items()
           if r["layout"] == "split_2" and acc not in exca]
    print("\norder#2 contiguous: %d ; split_2 ExcA-negative: %d" % (len(o2), len(sp2)))

    def occupant(acc):
        ty = best(tray, acc)
        if not ty:
            return None
        t = 1 if ty[1] == 1 else -1
        j = bypos.get((acc, ty[0] + t))
        if j is None or strand[j] != ty[1]:
            return None
        gap = (st[j] - en[bypos[(acc, ty[0])]]) if t == 1 else (st[bypos[(acc, ty[0])]] - en[j])
        return (j, abs(gap))

    print("\n=== (1) order#2: what is at TraY+1? ===")
    seqs2, gaps2 = [], []
    for acc in o2:
        r = occupant(acc)
        if r:
            seqs2.append((acc, aas[r[0]])); gaps2.append(r[1])
    print("  TraY+1 resolved on %d/%d" % (len(seqs2), len(o2)))
    if seqs2:
        uq = collections.Counter(s for _, s in seqs2)
        print("  unique sequences: %d ; top multiplicities %s"
              % (len(uq), [n for _, n in uq.most_common(5)]))
        L = sorted(len(s) for _, s in seqs2)
        print("  length: median %d Q1 %d Q3 %d" % (statistics.median(L), L[len(L)//4], L[3*len(L)//4]))
        lb = sum(1 for _, s in seqs2 if LIPO.search(s[:40]))
        print("  lipobox-positive: %d/%d = %.1f%%" % (lb, len(seqs2), 100.0 * lb / len(seqs2)))
        for s, n in uq.most_common(3):
            accs = [a2 for a2, x in seqs2 if x == s]
            cl = len({find(ix[a2]) for a2 in accs if a2 in ix})
            m = LIPO.search(s[:40])
            print("    n=%-4d %3d aa  %2d Mash clusters  lipobox %-8s %s"
                  % (n, len(s), cl, m.group() if m else "no", s[:34]))
        allcl = len({find(ix[a2]) for a2, _ in seqs2 if a2 in ix})
        print("  the %d records span %d Mash clusters" % (len(seqs2), allcl))
    nrel = sum(1 for acc in o2 if acc in relax)
    ncp = sum(1 for acc in o2 if acc in t4cp)
    print("  degeneracy check: relaxase on %d/%d (%.0f%%), t4cp2 on %d/%d (%.0f%%)"
          % (nrel, len(o2), 100.0 * nrel / len(o2), ncp, len(o2), 100.0 * ncp / len(o2)))

    print("\n=== (2) split_2: is TraY+1 inside the same block? ===")
    ok = far = unres = 0
    gaps = []
    for acc in sp2:
        r = occupant(acc)
        if not r:
            unres += 1; continue
        gaps.append(r[1])
        if r[1] < 20000:
            ok += 1
        else:
            far += 1
    print("  TraY+1 unresolved      %d" % unres)
    print("  within 20 kb of TraY   %d" % ok)
    print("  ACROSS a >=20 kb gap   %d  <- position meaningless there" % far)
    if gaps:
        g = sorted(gaps)
        print("  TraY-to-neighbour bp: median %d Q1 %d Q3 %d max %d"
              % (statistics.median(g), g[len(g)//4], g[3*len(g)//4], g[-1]))


if __name__ == "__main__":
    main()
