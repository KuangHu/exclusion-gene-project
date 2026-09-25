#!/usr/bin/env python3
"""Calibrate the in-slot/total threshold against known families and a null.

Round 1's payload is: build a model per candidate family, search the whole
database, and ask what fraction of its homologs sit in the slot.

    in-slot / total ~ 1   slot-specific -> a genuinely position-defined family
    clearly < 1           a general small-protein family; position caught coincidence

That ratio is meaningless without two calibration points, because neither end of
the scale is known a priori:

  UPPER  the named exclusion families (NF033894 Eex_IncN, TIGR04359 TrbK_RP4).
         These are real slot genes, so their ratio is what "slot-specific" looks
         like in this data -- including however much leakage is normal.

  LOWER  families nominated the same way at VirB10, the negative anchor that
         produced 62 families at 0.7% lipobox. Their ratio is what a positionally
         nominated NON-family looks like. This is the null the threshold sits above.

Without the lower point a ratio of 0.6 is unreadable; with it, 0.6 is either well
above or indistinguishable from what an arbitrary anchor produces.
"""
import argparse
import collections
import csv
import glob
import os
import statistics
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
E_CLUSTER = 1e-05
RESOLUTION = 0.05
E_SEARCH = 1e-05          # same threshold as clustering, so "homolog" means one thing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--families", type=int, default=8)
    a = ap.parse_args()
    import igraph
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner, idx, strand, aas = [], [], [], [], []
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
                strand.append(int(f[4])); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    print("cached proteins %d over %d accessions" % (len(seqs), len(set(owner))), flush=True)
    bypos = {(owner[i], idx[i]): i for i in range(len(owner))}

    def scan(fam, ev=None, minaa=0):
        with pyhmmer.plan7.HMMFile(os.path.join(HMM, fam + ".hmm")) as fh:
            m = next(iter(fh))
        kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
        out = collections.defaultdict(list)
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, **kw):
            for h in top:
                i = int(dec(h.name))
                if minaa and len(aas[i]) <= minaa:
                    continue
                out[owner[i]].append((idx[i], strand[i], h.score, i))
        return out
    B4, B5, B6 = scan("anchor_CagE_TrbE_VirB"), scan("anchor_T4SS", 1e-5, 150), scan("TrbL")
    B10 = scan("anchor_TrbI")
    best = lambda d, acc: max(d[acc], key=lambda x: x[2]) if acc in d else None

    # which protein index is "the slot" on each plasmid, and which is "VirB10+1"
    slot_i, neg_i = {}, {}
    for acc in B6:
        b6, b4 = best(B6, acc), best(B4, acc)
        b5 = best(B5, acc)
        t = 1 if b6[1] == 1 else -1
        if b5 and b5[1] == b6[1]:
            j = bypos.get((acc, b5[0] + t))          # +1 from VirB5
        elif b4 and b4[1] == b6[1]:
            arch = "canonical" if b4[0] * t < b6[0] * t else "rearranged"
            j = bypos.get((acc, b6[0] + (-1 if arch == "canonical" else 2) * t))
        else:
            j = None
        if j is not None and strand[j] == b6[1]:
            slot_i[acc] = j
    for acc in B10:
        b = best(B10, acc)
        j = bypos.get((acc, b[0] + (1 if b[1] == 1 else -1)))
        if j is not None and strand[j] == b[1]:
            neg_i[acc] = j
    print("slot position resolved on %d plasmids; VirB10+1 on %d"
          % (len(slot_i), len(neg_i)), flush=True)

    rows = []

    def ratio(hit_idx, position_map, label, note):
        inpos = sum(1 for i in hit_idx if position_map.get(owner[i]) == i)
        tot = len(hit_idx)
        r = inpos / tot if tot else float("nan")
        print("  %-34s in-position %5d / %6d = %.3f" % (label, inpos, tot, r))
        rows.append({"family": label, "kind": note, "in_position": inpos,
                     "total_hits": tot, "in_slot_ratio": round(r, 4)})
        return r

    print("\n=== UPPER calibration: the named exclusion families ===")
    for fam in ("NF033894", "TIGR04359"):
        h = scan(fam)
        allidx = [gi for v in h.values() for _, _, _, gi in v]
        ratio(allidx, slot_i, fam, "known_exclusion_family")

    print("\n=== LOWER calibration: families nominated at VirB10 (the null) ===")
    nom = {}
    for acc, j in neg_i.items():
        nom.setdefault(aas[j], []).append(acc)
    uniq = list(nom)
    du = [pyhmmer.easel.TextSequence(name=str(k).encode(), sequence=s).digitize(alpha)
          for k, s in enumerate(uniq)]
    edges = {}
    for top in pyhmmer.phmmer(du, du, cpus=a.cpus, E=1.0):
        q = int(dec(top.query.name))
        for h in top:
            t2 = int(dec(h.name))
            if t2 == q:
                continue
            k = (q, t2) if q < t2 else (t2, q)
            if k not in edges or h.evalue < edges[k]:
                edges[k] = h.evalue
    g = igraph.Graph(n=len(uniq), edges=[k for k, w in edges.items() if w <= E_CLUSTER])
    mem = g.community_leiden(objective_function="CPM", resolution=RESOLUTION, n_iterations=10).membership
    byc = collections.defaultdict(list)
    for i, m in enumerate(mem):
        byc[m].append(i)
    fams = sorted([v for v in byc.values() if len(v) >= 2], key=len, reverse=True)[:a.families]
    print("  %d VirB10 families >=2 members; calibrating on the largest %d"
          % (sum(1 for v in byc.values() if len(v) >= 2), len(fams)), flush=True)
    nulls = []
    for n, v in enumerate(fams):
        reps = [uniq[i] for i in sorted(v, key=lambda i: -len(uniq[i]))[:a.reps]]
        q = [pyhmmer.easel.TextSequence(name=("r%d" % k).encode(), sequence=s).digitize(alpha)
             for k, s in enumerate(reps)]
        hit = set()
        for top in pyhmmer.phmmer(q, seqs, cpus=a.cpus, E=E_SEARCH):
            for h in top:
                hit.add(int(dec(h.name)))
        r = ratio(sorted(hit), neg_i, "VirB10_family_%d (n=%d)" % (n, len(v)), "null_negative_anchor")
        if r == r:
            nulls.append(r)
    if nulls:
        print("\n  NULL distribution: median %.3f  max %.3f  (n=%d families)"
              % (statistics.median(nulls), max(nulls), len(nulls)))
        print("  A candidate family is only slot-specific if its ratio exceeds this,")
        print("  not merely if it is 'high'.")
    dest = os.path.join(PROJ, "data", "anchors", "inslot_calibration.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
