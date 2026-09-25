#!/usr/bin/env python3
"""Parallel negative slots: what does this pipeline produce at an arbitrary anchor?

Round 0 output is 63 families with no named Eex member. That number is
uninterpretable on its own -- position-nominate the gene next to ANY conserved
anchor, cluster it, and some number of families falls out. The question is whether
63 is above what the machinery yields anywhere.

So run the identical procedure -- nominate by position, dedupe, phmmer all-vs-all,
Leiden at the seed-validated settings -- at anchors that have no exclusion-gene
hypothesis attached:

    TARGET    gene immediately 3' of VirB5   (the slot)
    NEGATIVE  gene immediately 3' of VirB8, VirB9, VirB10

VirB8/9/10 are chosen because their detection failure (4.9%, 5.0%, 4.7%) is close
to VirB5's role in the target rule and they all have GA-level Pfam families, so
the nomination step is equally reliable. They are structural channel components
with no reported exclusion gene adjacent to them.

This is the only thing that gives 63 a denominator, and it is harder evidence than
any statistical test, because it uses the same code path end to end.

Second output: lipobox rate per family. Slot occupants run 67.8% lipobox against
0.5% in matched control gaps, so it is the strongest discriminator available.
Families should separate into a lipoprotein mode and a non-lipoprotein mode; only
the former are exclusion-gene candidates.

Red line (v2 1.3): exclusion families are used ONLY to label a protein that
position has already nominated. They never choose an anchor or bound a window.
"""
import argparse
import collections
import csv
import hashlib
import os
import re
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
OUT = "/global/scratch/users/kh36969/exclusion_gene/negctrl"

ANCHORS = {                      # label -> (hmm, evalue, min_aa)
    "VirB5":  ("anchor_T4SS", 1e-5, 150),
    "VirB8":  ("anchor_VirB8", None, 0),
    "VirB9":  ("anchor_CagX", None, 0),
    "VirB10": ("anchor_TrbI", None, 0),
}
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "TraS", "DUF4467"]
LIPOBOX = re.compile(r"[LVI][ASTVIG][GASN]C")
E_CLUSTER = 1e-05                # seed-control validated window is 1e-03..1e-05
RESOLUTION = 0.05                # PASS across 0.005-0.1; 0.05 is mid-range


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import igraph
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    os.makedirs(OUT, exist_ok=True)
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    # ---- load the frozen CDS cache -----------------------------------------
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
    print("cached proteins %d over %d accessions"
          % (len(seqs), len(set(owner))), flush=True)
    bypos = {}
    for i in range(len(owner)):
        bypos[(owner[i], idx[i])] = i

    def scan(fam, ev, minaa):
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

    eexhit = {}
    for fam in EEX:
        for acc, v in scan(fam, None, 0).items():
            for _, _, _, i in v:
                eexhit[i] = fam
    print("Eex-family labelled proteins: %d" % len(eexhit), flush=True)

    rows = []
    for lab, (fam, ev, minaa) in ANCHORS.items():
        hits = scan(fam, ev, minaa)
        nominated = {}
        for acc, v in hits.items():
            ai, ast, _, _ = max(v, key=lambda x: x[2])
            nxt = ai + 1 if ast == 1 else ai - 1     # immediately 3' in transcription
            j = bypos.get((acc, nxt))
            if j is None or strand[j] != ast:
                continue
            nominated.setdefault(aas[j], []).append(j)
        uniq = list(nominated)
        print("\n=== %s : %d plasmids anchored, %d unique nominated proteins ==="
              % (lab, len(hits), len(uniq)), flush=True)
        if len(uniq) < 10:
            continue
        du = [pyhmmer.easel.TextSequence(name=str(k).encode(), sequence=s).digitize(alpha)
              for k, s in enumerate(uniq)]
        edges = {}
        for top in pyhmmer.phmmer(du, du, cpus=a.cpus, E=1.0):
            q = int(dec(top.query.name))
            for h in top:
                t = int(dec(h.name))
                if t == q:
                    continue
                k = (q, t) if q < t else (t, q)
                if k not in edges or h.evalue < edges[k]:
                    edges[k] = h.evalue
        el = [k for k, w in edges.items() if w <= E_CLUSTER]
        g = igraph.Graph(n=len(uniq), edges=el)
        part = g.community_leiden(objective_function="CPM",
                                  resolution=RESOLUTION, n_iterations=10)
        mem = part.membership
        byc = collections.defaultdict(list)
        for i, m in enumerate(mem):
            byc[m].append(i)
        named = {i for i, s in enumerate(uniq)
                 if any(j in eexhit for j in nominated[s])}
        fams = [v for v in byc.values() if len(v) >= 2]
        novel = [v for v in fams if not (set(v) & named)]
        lipo_all = sum(1 for s in uniq if LIPOBOX.search(s[:40]))
        # per-family lipobox rate
        rates = []
        for v in novel:
            r = sum(1 for i in v if LIPOBOX.search(uniq[i][:40])) / len(v)
            rates.append((r, len(v), sorted(len(uniq[i]) for i in v)[len(v) // 2]))
        hi = [x for x in rates if x[0] >= 0.5]
        print("  unique nominated      %d" % len(uniq))
        print("  named-Eex members     %d" % len(named))
        print("  clusters (>=2)        %d" % len(fams))
        print("  families with NO named member  %d" % len(novel))
        print("  of those, lipobox rate >=50%%  %d" % len(hi))
        print("  overall lipobox rate  %.1f%%" % (100.0 * lipo_all / len(uniq)))
        rows.append({"anchor": lab, "role": "TARGET" if lab == "VirB5" else "negative",
                     "plasmids_anchored": len(hits), "unique_nominated": len(uniq),
                     "named_eex_members": len(named), "clusters_ge2": len(fams),
                     "families_no_named_member": len(novel),
                     "families_lipobox_ge50pct": len(hi),
                     "overall_lipobox_pct": round(100.0 * lipo_all / len(uniq), 1)})
        with open(os.path.join(OUT, "families_%s.tsv" % lab), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["family_size", "lipobox_rate", "median_aa"])
            for r, n, med in sorted(rates, reverse=True):
                w.writerow([n, round(r, 3), med])

    print("\n" + "=" * 78)
    print("NEGATIVE-SLOT CONTROL -- does the pipeline yield families anywhere?")
    print("=" * 78)
    print("%-8s %-9s %10s %8s %12s %14s %8s"
          % ("anchor", "role", "unique", "named", "families", "lipobox>=50%", "lipo%"))
    for r in rows:
        print("%-8s %-9s %10d %8d %12d %14d %8.1f"
              % (r["anchor"], r["role"], r["unique_nominated"], r["named_eex_members"],
                 r["families_no_named_member"], r["families_lipobox_ge50pct"],
                 r["overall_lipobox_pct"]))
    dest = os.path.join(PROJ, "data", "anchors", "negative_slot_control.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
