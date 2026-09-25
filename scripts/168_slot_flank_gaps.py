#!/usr/bin/env python3
"""How much non-coding DNA actually flanks a slot occupant?

The pilot scored only 7 upstream and 9 downstream spacers across 120 plasmids.
That is consistent with TWO OPPOSITE explanations and the pilot cannot tell them
apart:

  (a) the gaps are genuinely tiny -- genes are packed, there is no room for a
      terminator. That is a REAL result about operon architecture.
  (b) the gene-index matching failed and nothing was measured. That is a bug.

This measures the distance directly, with no folding and no scoring, so the two
are separated before any conclusion is drawn about detectability.

Reported per slot occupant:
  gap_up    bases between the end of the upstream gene and the start of the slot
  gap_down  bases between the end of the slot and the start of the downstream gene
  and the same for anchor-anchor pairs as the internal reference.

A Rho-independent terminator needs roughly 30-50 nt (stem + loop + poly-U). If
the flanking gaps are mostly below that, "no terminator detected" is not a
statement about detection sensitivity -- it is a statement that the DNA is not
there.
"""
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CAT = os.path.join(PROJ, "data", "release", "v1.1", "catalogue_MPF_T_v1.1.tsv")
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
ANCH = ["VirB1", "VirB2", "VirB3", "VirB4", "VirB5", "VirB6", "VirB8",
        "VirB9", "VirB10", "VirB11", "VirD4"]
TERM_MIN = 30      # minimum room a Rho-independent terminator needs


def parse_cell(c):
    if not c or c == "none":
        return None
    f = c.split(":")
    if f[0] == "ga":
        f = f[1:]
    s, e = f[1].split("..")
    return int(f[0]), int(s), int(e), 1 if f[2] == "+" else -1


def main():
    require_compute_node()
    cat = {r["accession"]: r for r in csv.DictReader(open(CAT), delimiter="\t")
           if r["slot_ready__virb5_virb6"] == "1"}
    sof = {r["accession"]: r for r in csv.DictReader(open(SOF), delimiter="\t")
           if r["slot_status"] in ("eex_occupied", "candidate") and r["slot_occupant_coords"]}
    accs = [a for a in sof if a in cat]
    print("plasmids with a resolved slot: %d" % len(accs), flush=True)

    # gene coordinates straight from the cache, indexed the SAME way the
    # catalogue indexes them -- no re-calling, so no index drift
    genes = collections.defaultdict(dict)
    want = set(accs)
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                if f[0] in want:
                    genes[f[0]][int(f[1])] = (int(f[2]), int(f[3]), int(f[4]))
    print("  gene coordinates loaded for %d plasmids\n" % len(genes), flush=True)

    rows = []
    miss = collections.Counter()
    for acc in accs:
        g = genes.get(acc)
        if not g:
            miss["no_genes"] += 1
            continue
        p5 = parse_cell(cat[acc].get("VirB5", ""))
        if not p5:
            miss["no_virb5"] += 1
            continue
        gi5, _, _, sd = p5
        slot_gi = gi5 + sd
        if slot_gi not in g:
            miss["no_slot_gene"] += 1
            continue
        ss, se, _ = g[slot_gi]
        # neighbours in COORDINATE order, which is what matters for spacing
        prev_gi, next_gi = slot_gi - 1, slot_gi + 1
        gu = gd = None
        if prev_gi in g:
            gu = ss - g[prev_gi][1] - 1
        if next_gi in g:
            gd = g[next_gi][0] - se - 1
        # internal reference: gaps between adjacent anchors on this plasmid
        idx = []
        for a in ANCH:
            c = parse_cell(cat[acc].get(a, ""))
            if c:
                idx.append(c[0])
        idx.sort()
        anchgaps = []
        for j in range(len(idx) - 1):
            if idx[j + 1] - idx[j] == 1 and idx[j] in g and idx[j + 1] in g:
                anchgaps.append(g[idx[j + 1]][0] - g[idx[j]][1] - 1)
        rows.append({"accession": acc, "slot_family": sof[acc]["slot_occupant_family"],
                     "slot_status": sof[acc]["slot_status"],
                     "gap_up": gu if gu is not None else "",
                     "gap_down": gd if gd is not None else "",
                     "anchor_gap_median": (sorted(anchgaps)[len(anchgaps) // 2]
                                           if anchgaps else ""),
                     "n_anchor_gaps": len(anchgaps)})
    print("=== resolution ===")
    print("  slot gene resolved on %d of %d plasmids" % (len(rows), len(accs)))
    for k, v in miss.most_common():
        print("    unresolved: %-14s %d" % (k, v))

    def band(v):
        v = sorted(x for x in v if x != "" and x is not None)
        if not v:
            return None
        n = len(v)
        return (n, v[n // 4], v[n // 2], v[3 * n // 4],
                100.0 * sum(1 for x in v if x >= TERM_MIN) / n,
                100.0 * sum(1 for x in v if x <= 0) / n)

    print("\n=== FLANKING NON-CODING DISTANCE (bp) ===")
    print("  %-20s %7s %7s %8s %7s %11s %10s"
          % ("region", "n", "Q1", "median", "Q3", ">=%dbp" % TERM_MIN, "overlap<=0"))
    for lab, key in (("upstream of slot", "gap_up"), ("downstream of slot", "gap_down"),
                     ("anchor-anchor (ref)", "anchor_gap_median")):
        b = band([r[key] for r in rows])
        if not b:
            print("  %-20s %7d" % (lab, 0)); continue
        n, q1, md, q3, ge, ov = b
        print("  %-20s %7d %7d %8d %7d %10.1f%% %9.1f%%" % (lab, n, q1, md, q3, ge, ov))

    print("\n=== READING ===")
    bu = band([r["gap_up"] for r in rows])
    bd = band([r["gap_down"] for r in rows])
    if bu and bd:
        print("  A Rho-independent terminator needs ~%d nt of room." % TERM_MIN)
        print("  upstream   %.1f%% of slots have that much" % bu[4])
        print("  downstream %.1f%% of slots have that much" % bd[4])
        if bu[4] < 25 and bd[4] < 25:
            print("\n  -> The flanking DNA IS NOT THERE on most plasmids.")
            print("     'No terminator detected' would be a statement about")
            print("     ARCHITECTURE, not about detection sensitivity.")
        else:
            print("\n  -> Enough plasmids have room; a detection test is meaningful")
            print("     and the pilot's low n was an indexing failure, not biology.")

    # named vs unnamed, in case the architecture differs
    print("\n=== by slot status ===")
    for st in ("eex_occupied", "candidate"):
        sel = [r for r in rows if r["slot_status"] == st]
        b1 = band([r["gap_up"] for r in sel]); b2 = band([r["gap_down"] for r in sel])
        if b1 and b2:
            print("  %-14s n=%4d  up median %4d (%.0f%% >=%d)  down median %4d (%.0f%% >=%d)"
                  % (st, len(sel), b1[2], b1[4], TERM_MIN, b2[2], b2[4], TERM_MIN))

    dest = os.path.join(PROJ, "data", "anchors", "slot_flank_gaps.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))


if __name__ == "__main__":
    sys.exit(main() or 0)
