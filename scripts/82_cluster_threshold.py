#!/usr/bin/env python3
"""Pipeline v2 §0.1 -- choose the ANI clustering threshold from the data.

Two questions, per §0.1 and the review:

  1. Where is the valley in the pairwise ANI distribution? The threshold goes
     there, rather than at the conventional 95%.
  2. At each candidate threshold, do the positive controls land in SEPARATE
     clusters? RP4 (IncP-alpha) and R751 (IncP-beta) merging is DISQUALIFYING:
     their TrbK proteins are only 37.8% identical, so they are plainly
     independent evolutionary events and must not be counted as one observation.

Clustering is single-linkage (connected components) on the sparse graph of pairs
below the threshold. Single linkage is the permissive choice, which is the
conservative direction here: it merges more, so a threshold that keeps the
controls separate under single linkage keeps them separate under any stricter
linkage.

Note pKM101 cannot be tested -- it has no complete plasmid sequence, so it is
not in PLSDB and cannot be sketched. RP4 is a TPA record, also absent from
PLSDB, and was injected into the sketch explicitly.
"""
import gzip
import os
import sys
from collections import Counter

PAIRS = "/global/scratch/users/kh36969/exclusion_gene/mash/pairs_d010.tsv.gz"
IDS = "/global/scratch/users/kh36969/exclusion_gene/mash/all_ids.txt"

# The controls, and what a correct threshold must do with them.
CONTROLS = {
    "BN000925.1":  ("RP4",  "IncP-alpha"),
    "NC_001735.4": ("R751", "IncP-beta"),
    "NC_028464.1": ("R388", "IncW"),
}
THRESHOLDS = [0.001, 0.002, 0.003, 0.005, 0.007, 0.010, 0.020, 0.030, 0.050]
LOAD_MAX = 0.050        # only pairs this close are ever needed


class UF:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def main():
    ids = [l.strip() for l in open(IDS) if l.strip()]
    idx = {a: i for i, a in enumerate(ids)}
    n = len(ids)
    print("sketches: %d" % n)
    for acc, (name, inc) in CONTROLS.items():
        print("  control %-8s %-12s %s" % (name, inc,
              "in sketch" if acc in idx else "*** NOT IN SKETCH ***"))

    print("\nloading pairs with d <= %.3f ..." % LOAD_MAX)
    A, B, D = [], [], []
    kept = seen = 0
    with gzip.open(PAIRS, "rt") as fh:
        for line in fh:
            seen += 1
            f = line.split("\t", 3)
            d = float(f[2])
            if d > LOAD_MAX:
                continue
            a, b = idx.get(f[0]), idx.get(f[1])
            if a is None or b is None:
                continue
            A.append(a); B.append(b); D.append(d)
            kept += 1
    print("  read %d pairs, kept %d (d <= %.3f)" % (seen, kept, LOAD_MAX))

    print("\n%-8s %-9s %9s %9s %9s %9s   %s"
          % ("d", "ANI %", "clusters", "singleton", "largest", "top-10 %", "controls"))
    print("-" * 96)
    rows = []
    for t in THRESHOLDS:
        uf = UF(n)
        for a, b, d in zip(A, B, D):
            if d <= t:
                uf.union(a, b)
        comp = Counter(uf.find(i) for i in range(n))
        sizes = sorted(comp.values(), reverse=True)
        nclust = len(sizes)
        single = sum(1 for s in sizes if s == 1)
        top10 = 100.0 * sum(sizes[:10]) / n
        # control separation
        croots = {}
        for acc, (name, inc) in CONTROLS.items():
            if acc in idx:
                croots.setdefault(uf.find(idx[acc]), []).append(name)
        merged = [v for v in croots.values() if len(v) > 1]
        if merged:
            cstat = "MERGED: " + " + ".join("/".join(v) for v in merged)
        else:
            cstat = "separate (%d/%d)" % (len(croots), len(CONTROLS))
        print("%-8.3f %-9.1f %9d %9d %9d %8.1f%%   %s"
              % (t, (1 - t) * 100, nclust, single, sizes[0], top10, cstat))
        rows.append({"d": t, "ani": (1 - t) * 100, "clusters": nclust,
                     "singletons": single, "largest": sizes[0],
                     "top10_pct": round(top10, 2), "controls": cstat,
                     "sizes": sizes})

    print("\n=== cluster size distribution at each threshold ===")
    for r in rows:
        s = r["sizes"]
        q = lambda p: s[min(len(s) - 1, int(len(s) * p))]
        big = sum(1 for x in s if x >= 100)
        print("  d=%.3f (ANI %.1f%%)  clusters=%-6d  largest=%-5d  "
              ">=100 members: %-4d  median=%d"
              % (r["d"], r["ani"], r["clusters"], r["largest"], big,
                 s[len(s) // 2]))

    import csv
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "data", "anchors", "cluster_threshold_scan.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["d", "ani", "clusters", "singletons",
                                           "largest", "top10_pct", "controls"],
                           delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("\nwrote %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
