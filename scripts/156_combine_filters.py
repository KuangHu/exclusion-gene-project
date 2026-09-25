#!/usr/bin/env python3
"""Intersect the two filters, and stratify the survivors by cluster size.

TWO QUESTIONS, both cheap:

1. DO THE FILTERS OVERLAP? The isolated-ORF filter (161) was computed on the
   2,297 remainder BEFORE the structural filter (612), so 2,297-612-161=1,524 is
   an UPPER BOUND, not a count. If the overlap is large, the isolated singletons
   were mostly distant known-family members anyway and the recurrence filter adds
   little independent information.

2. SIZE STRATA. The survivors are not one population:
     large  (>=31)  prevalence too high for an exclusion gene -> category 3,
                    the cluster-22 situation
     medium (5-30)  the band where real families live: TrbK is 218 in plasmids
                    and 10 in ICEs; this is the size an undescribed family has
     small  (2-4)   past the isolation filter but still consistent with artefact

   Only the medium band plausibly needs individual inspection, and the three
   functional-constraint screens are affordable on that band and not on 1,524.
"""
import collections, csv, gzip, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
PAIRS = "/global/scratch/users/kh36969/exclusion_gene/mash/pairs_d010.tsv.gz"
FSK = os.path.join(PROJ, "data", "anchors", "foldseek_remainder.tsv")
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"]


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
    dig = lambda L: [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
                     for i, s in enumerate(L)]
    U = dig(uniq); known = set()
    for f in EEX:
        p = os.path.join(HMM, f + ".hmm")
        if not os.path.exists(p): continue
        with pyhmmer.plan7.HMMFile(p) as fh: m = next(iter(fh))
        for top in pyhmmer.hmmsearch([m], U, cpus=8, E=1e-3):
            for h in top: known.add(int(dec(h.name)))
    rem = [uniq[i] for i in range(len(uniq)) if i not in known]
    print("remainder: %d" % len(rem), flush=True)

    fs = list(csv.DictReader(open(FSK), delimiter="\t"))
    assert len(fs) == len(rem), "index mismatch: foldseek %d vs rem %d" % (len(fs), len(rem))
    struct = {int(r["query_idx"]) for r in fs if r["assigned"] == "1"}
    fam_of = {int(r["query_idx"]): r["best_named_family"] for r in fs if r["assigned"] == "1"}
    print("  structural known-family: %d" % len(struct), flush=True)

    R = dig(rem); best = {}
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
                if float(f[2]) > 0.007: continue
            except ValueError: continue
            if f[0] in par and f[1] in par:
                ra, rb = find(f[0]), find(f[1])
                if ra != rb: par[ra] = rb

    iso = set()
    info = {}
    for c, idxs in cl.items():
        accs = set()
        for i in idxs: accs |= acc_of[rem[i]]
        mash = len({find(a) for a in accs})
        info[c] = (len(idxs), len(accs), mash)
        if mash < 2: iso |= set(idxs)

    both = struct & iso
    print("\n=== 1. DO THE FILTERS OVERLAP? ===")
    print("  structural known-family      %5d" % len(struct))
    print("  isolated on one element      %5d" % len(iso))
    print("  BOTH                         %5d" % len(both))
    print("  union                        %5d" % len(struct | iso))
    print("  survivors = %d - %d           %5d" % (len(rem), len(struct | iso), len(rem)-len(struct|iso)))
    print("  -> the recurrence filter adds %d sequences the structural filter did not"
          % (len(iso) - len(both)))
    if len(iso):
        print("  -> %.0f%% of isolated singletons were ALSO structurally known-family"
              % (100.0*len(both)/len(iso)))

    surv = [i for i in range(len(rem)) if i not in struct and i not in iso]
    print("\n=== 2. SIZE STRATA of the %d survivors ===" % len(surv))
    bands = {"large (>=31)": [], "medium (5-30)": [], "small (2-4)": [], "singleton": []}
    for c, idxs in cl.items():
        keep = [i for i in idxs if i in set(surv)]
        if not keep: continue
        n = len(keep)
        b = ("large (>=31)" if n >= 31 else "medium (5-30)" if n >= 5
             else "small (2-4)" if n >= 2 else "singleton")
        bands[b].append((c, n, info[c][1], info[c][2]))
    print("  %-14s %8s %9s %11s" % ("band", "clusters", "unique", "Mash cl (sum)"))
    for b in ("large (>=31)", "medium (5-30)", "small (2-4)", "singleton"):
        v = bands[b]
        print("  %-14s %8d %9d %11d" % (b, len(v), sum(x[1] for x in v), sum(x[3] for x in v)))
    print("\n  === MEDIUM BAND (5-30) -- the band where a real family would sit ===")
    print("  %-8s %7s %9s %10s %10s" % ("cluster", "unique", "records", "Mash cl", "median_aa"))
    med = sorted(bands["medium (5-30)"], key=lambda x: -x[3])
    for c, n, nrec, mash in med[:25]:
        idxs = [i for i in cl[c] if i in set(surv)]
        L = sorted(len(rem[i]) for i in idxs)
        print("  %-8d %7d %9d %10d %10d" % (c, n, nrec, mash, L[len(L)//2]))
    print("\n  medium band: %d clusters, %d unique sequences -- THIS is the set that"
          % (len(med), sum(x[1] for x in med)))
    print("  the functional-constraint screens can afford to run on.")
    print("\n  === LARGE BAND -- prevalence too high for an exclusion gene ===")
    for c, n, nrec, mash in sorted(bands["large (>=31)"], key=lambda x: -x[3])[:8]:
        idxs = [i for i in cl[c] if i in set(surv)]
        L = sorted(len(rem[i]) for i in idxs)
        print("  cluster %-6d unique %4d  records %5d  Mash %5d  median %d aa"
              % (c, n, nrec, mash, L[len(L)//2]))
    dest = os.path.join(PROJ, "data", "anchors", "candidate_bands.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["band", "cluster", "n_unique_surviving", "n_records", "n_mash_clusters", "median_aa"])
        for b, v in bands.items():
            for c, n, nrec, mash in v:
                idxs = [i for i in cl[c] if i in set(surv)]
                L = sorted(len(rem[i]) for i in idxs)
                w.writerow([b, c, n, nrec, mash, L[len(L)//2]])
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
