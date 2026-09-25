#!/usr/bin/env python3
"""Size bands on the CORRECTED candidate pool (1,734), not the contaminated one.

The existing bands (large 894 / medium 460 / small 172) were computed before the
anchor contamination was removed, so they are stale. Cluster 21 in particular --
193 unique over 1,196 Mash clusters, median 185 aa -- is almost certainly VirB6:
TrbL was both the largest Pfam hit (324) and the largest contaminant (460).

This rebuilds the pool exactly as script 158 defines it (decontamination is part
of the slot-occupant definition, and the exclusion family set is the corrected
nine) and then bands what is left.

WHY BANDS. The pool needs a magnitude, not another binary cut:

  large  (>=31)  prevalence too high for an exclusion gene. An exclusion gene is
                 not the DEFAULT occupant of its own slot, and TrbK's mature form
                 is 47 aa, Eex ~75 aa -- the large clusters run much bigger.
  medium (5-30)  the plausible magnitude for a real family. Known families sit
                 near here: TrbK 218 in plasmids, DUF4467 98, EexR 38.
  small  (2-4)   past the isolation filter but still consistent with artefact.

Only the medium band is worth individual inspection, and only it is affordable
for the functional-constraint screens.

Independence is counted at the Mash-cluster level (d=0.007), never per record:
PLSDB over-samples clinical Enterobacteriaceae and one clonal expansion would
otherwise look like recurrence.
"""
import collections, csv, gzip, os, sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import decontaminate as DC

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
PAIRS = "/global/scratch/users/kh36969/exclusion_gene/mash/pairs_d010.tsv.gz"
# the corrected nine
ENTRY = ["TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729", "PF20084"]
SURFACE = ["PF05818", "PrgA_Sea1"]
D_MASH = 0.007


def main():
    require_compute_node()
    import pyhmmer, igraph
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    DC.assert_no_evalue_anchor()

    rows = [r for r in csv.DictReader(open(SOF), delimiter="\t")
            if r["slot_status"] in ("candidate", "eex_occupied") and r["slot_occupant_coords"]]
    want = collections.defaultdict(set)
    for r in rows:
        want[r["accession"]].add(r["slot_occupant_coords"])

    seq_of, allseq, allaa = {}, [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]; p = line.rstrip("\n")
                if acc in want:
                    k = "%s..%s" % (f[2], f[3])
                    if k in want[acc]:
                        seq_of[(acc, k)] = p
                    allaa.append(p)
                    allseq.append(pyhmmer.easel.TextSequence(
                        name=str(len(allaa) - 1).encode(), sequence=p).digitize(alpha))
    occ = collections.defaultdict(set)
    for r in rows:
        s = seq_of.get((r["accession"], r["slot_occupant_coords"]))
        if s:
            occ[s].add(r["accession"])
    uniq = sorted(occ)
    print("slot occupants: %d unique" % len(uniq), flush=True)

    dig = lambda L: [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
                     for i, s in enumerate(L)]
    U = dig(uniq)

    def scan(acc, ev=None):
        p = os.path.join(HMM, acc + ".hmm")
        if not os.path.exists(p):
            return set()
        with pyhmmer.plan7.HMMFile(p) as fh:
            m = next(iter(fh))
        got = set()
        try:
            kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
            for top in pyhmmer.hmmsearch([m], U, cpus=8, **kw):
                for h in top:
                    got.add(int(dec(h.name)))
        except pyhmmer.errors.MissingCutoffs:
            for top in pyhmmer.hmmsearch([m], U, cpus=8, E=1e-5):
                for h in top:
                    got.add(int(dec(h.name)))
        return got

    named = set()
    for a in ENTRY:
        named |= scan(a)
    for a in SURFACE:
        named |= scan(a, 1e-5)
    print("  named by the corrected nine: %d" % len(named), flush=True)

    members = DC.confirmed_members(allseq, allaa, cpus=8)
    kept, _ = DC.decontaminate(uniq, members, cpus=8, verbose=False)
    anchor = {i for i, s in enumerate(uniq) if s not in set(kept)}
    print("  anchor-family occupants: %d" % len(anchor), flush=True)

    pool = [i for i in range(len(uniq)) if i not in named and i not in anchor]
    seqs = [uniq[i] for i in pool]
    print("  CORRECTED CANDIDATE POOL: %d unique\n" % len(seqs), flush=True)

    R = dig(seqs)
    best = {}
    for top in pyhmmer.phmmer(R, R, cpus=8, E=10.0):
        i = int(dec(top.query.name))
        for h in top:
            j = int(dec(h.name))
            if i == j:
                continue
            k = (min(i, j), max(i, j))
            if k not in best or h.evalue < best[k]:
                best[k] = h.evalue
    g = igraph.Graph(n=len(seqs), edges=[(u, v) for (u, v), e in best.items() if e <= 1e-5])
    part = g.community_leiden(objective_function="CPM", resolution=0.01, n_iterations=-1)
    cl = collections.defaultdict(list)
    for i, c in enumerate(part.membership):
        cl[c].append(i)

    need = set()
    for s in seqs:
        need |= occ[s]
    par = {a: a for a in need}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x

    with gzip.open(PAIRS, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3:
                continue
            try:
                if float(f[2]) > D_MASH:
                    continue
            except ValueError:
                continue
            if f[0] in par and f[1] in par:
                ra, rb = find(f[0]), find(f[1])
                if ra != rb:
                    par[ra] = rb

    out = []
    for c, idxs in cl.items():
        accs = set()
        for i in idxs:
            accs |= occ[seqs[i]]
        mash = len({find(a) for a in accs})
        L = sorted(len(seqs[i]) for i in idxs)
        n = len(idxs)
        band = ("large (>=31)" if n >= 31 else "medium (5-30)" if n >= 5
                else "small (2-4)" if n >= 2 else "singleton")
        out.append({"band": band, "cluster": c, "n_unique": n, "n_records": len(accs),
                    "n_mash_clusters": mash, "median_aa": L[len(L) // 2],
                    "recurs": int(mash >= 2)})
    out.sort(key=lambda r: (-r["n_mash_clusters"], -r["n_unique"]))

    print("=== SIZE BANDS on the corrected pool ===")
    print("  %-14s %9s %9s %11s" % ("band", "clusters", "unique", "Mash cl"))
    for b in ("large (>=31)", "medium (5-30)", "small (2-4)", "singleton"):
        v = [r for r in out if r["band"] == b]
        print("  %-14s %9d %9d %11d"
              % (b, len(v), sum(r["n_unique"] for r in v), sum(r["n_mash_clusters"] for r in v)))

    med = [r for r in out if r["band"] == "medium (5-30)" and r["recurs"]]
    print("\n=== MEDIUM BAND (5-30, recurring) -- the inspection set ===")
    print("  %-8s %7s %9s %10s %10s" % ("cluster", "unique", "records", "Mash cl", "median_aa"))
    for r in med[:30]:
        print("  %-8s %7d %9d %10d %10d"
              % (r["cluster"], r["n_unique"], r["n_records"], r["n_mash_clusters"], r["median_aa"]))
    print("\n  MEDIUM BAND: %d clusters / %d unique sequences" %
          (len(med), sum(r["n_unique"] for r in med)))
    print("  This is the set the functional-constraint screens can afford.")

    lg = [r for r in out if r["band"] == "large (>=31)"]
    if lg:
        print("\n=== LARGE BAND -- prevalence wrong for an exclusion gene ===")
        for r in lg[:8]:
            print("  cluster %-6s unique %4d  records %5d  Mash %5d  median %d aa"
                  % (r["cluster"], r["n_unique"], r["n_records"], r["n_mash_clusters"],
                     r["median_aa"]))

    dest = os.path.join(PROJ, "data", "anchors", "bands_corrected.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
