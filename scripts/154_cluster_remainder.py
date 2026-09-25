#!/usr/bin/env python3
"""Cluster the unnamed remainder and apply the recurrence criterion.

Operational definition of "new", as set: a candidate must
  * have no HMM hit to any known exclusion family (E > 1e-3)   -- already applied
  * occupy a slot that RECURS across independent elements
  * and those occupants must cluster into an internally consistent family

The last two are what this script measures. They also dispose of categories 3 and
4 without any new tool: an ORF-calling artefact or a one-off pseudogene does not
recur on independent elements, and a pilus-accessory or regulatory protein that
does recur will at least form a coherent family that can then be named.

**An isolated ORF on a single element is not a discovery.** Same argument as
covariation needing a population: one sequence carries no evidence of anything.

INDEPENDENCE IS COUNTED AT THE MASH-CLUSTER LEVEL, not per record. PLSDB is
heavily over-sampled for clinical Enterobacteriaceae -- 8,515 named records
collapse to 801 unique sequences (10.6x). Counting records would let one clonal
expansion look like recurrence.

Clustering: phmmer all-vs-all, E<=1e-5, Leiden CPM resolution 0.01 -- the
configuration that passed the seed control in the MPF_T slot iteration. NOT the
0.05 originally specified, which splits the Eex_IncN seed pair.
"""
import collections, csv, gzip, math, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
PAIRS = "/global/scratch/users/kh36969/exclusion_gene/mash/pairs_d010.tsv.gz"
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"]
D_MASH = 0.007


def main():
    require_compute_node()
    import pyhmmer, igraph
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    rows = [r for r in csv.DictReader(open(SOF), delimiter="\t")
            if r["slot_status"] == "candidate" and r["slot_occupant_coords"]]
    want = collections.defaultdict(set)
    for r in rows: want[r["accession"]].add(r["slot_occupant_coords"])
    seq_of = {}
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                if acc in want:
                    k = "%s..%s" % (f[2], f[3])
                    if k in want[acc]: seq_of[(acc, k)] = line.rstrip("\n")
    acc_of = collections.defaultdict(set)
    for r in rows:
        s = seq_of.get((r["accession"], r["slot_occupant_coords"]))
        if s: acc_of[s].add(r["accession"])
    uniq = sorted(acc_of)
    print("unnamed candidates: %d records -> %d unique sequences"
          % (len(rows), len(uniq)), flush=True)

    # drop the 279 already shown to be distant members of known families
    dig = lambda L: [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
                     for i, s in enumerate(L)]
    U = dig(uniq)
    known = set()
    for f in EEX:
        p = os.path.join(HMM, f + ".hmm")
        if not os.path.exists(p): continue
        with pyhmmer.plan7.HMMFile(p) as fh: m = next(iter(fh))
        for top in pyhmmer.hmmsearch([m], U, cpus=8, E=1e-3):
            for h in top: known.add(int(dec(h.name)))
    rem = [uniq[i] for i in range(len(uniq)) if i not in known]
    print("  minus %d distant known-family members -> %d remainder"
          % (len(known), len(rem)), flush=True)

    R = dig(rem)
    best = {}
    for top in pyhmmer.phmmer(R, R, cpus=8, E=10.0):
        i = int(dec(top.query.name))
        for h in top:
            j = int(dec(h.name))
            if i == j: continue
            k = (min(i, j), max(i, j))
            if k not in best or h.evalue < best[k]: best[k] = h.evalue
    g = igraph.Graph(n=len(rem), edges=[(u, v) for (u, v), e in best.items() if e <= 1e-5])
    part = g.community_leiden(objective_function="CPM", resolution=0.01, n_iterations=-1)
    cl = collections.defaultdict(list)
    for i, c in enumerate(part.membership): cl[c].append(i)
    print("  clusters: %d" % len(cl), flush=True)

    need = set()
    for s in rem: need |= acc_of[s]
    par = {a: a for a in need}
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    with gzip.open(PAIRS, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3: continue
            try:
                if float(f[2]) > D_MASH: continue
            except ValueError: continue
            if f[0] in par and f[1] in par:
                ra, rb = find(f[0]), find(f[1])
                if ra != rb: par[ra] = rb
    out = []
    for c, idxs in cl.items():
        seqs = [rem[i] for i in idxs]
        accs = set()
        for s in seqs: accs |= acc_of[s]
        mash = len({find(a) for a in accs})
        out.append({"cluster": c, "n_unique": len(seqs), "n_records": len(accs),
                    "n_mash_clusters": mash,
                    "median_aa": sorted(len(s) for s in seqs)[len(seqs)//2],
                    "medoid_aa": len(max(seqs, key=len)),
                    "recurs": int(mash >= 2), "seq": max(seqs, key=len)[:60]})
    out.sort(key=lambda r: (-r["n_mash_clusters"], -r["n_unique"]))
    dest = os.path.join(PROJ, "data", "anchors", "remainder_clusters.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(out)

    rec = [r for r in out if r["recurs"]]
    iso = [r for r in out if not r["recurs"]]
    print("\n=== THE RECURRENCE FILTER ===")
    print("  clusters spanning >=2 independent Mash clusters : %d" % len(rec))
    print("    unique sequences in them                     : %d" % sum(r["n_unique"] for r in rec))
    print("  clusters confined to ONE Mash cluster          : %d  <- ELIMINATED" % len(iso))
    print("    unique sequences in them                     : %d" % sum(r["n_unique"] for r in iso))
    sing = [r for r in out if r["n_unique"] == 1 and not r["recurs"]]
    print("  of the eliminated, singletons on one element    : %d" % len(sing))
    print("\n=== surviving clusters, ranked by independence ===")
    print("  %-8s %7s %8s %9s %9s" % ("cluster", "uniq", "records", "Mash cl", "median_aa"))
    for r in rec[:20]:
        print("  %-8d %7d %8d %9d %9d"
              % (r["cluster"], r["n_unique"], r["n_records"], r["n_mash_clusters"], r["median_aa"]))
    print("\n  candidate pool after recurrence filter: %d clusters / %d unique sequences"
          % (len(rec), sum(r["n_unique"] for r in rec)))
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
