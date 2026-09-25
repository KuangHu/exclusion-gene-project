#!/usr/bin/env python3
"""Cluster the slot proteins and ask whether the candidates form families.

Input is the gene immediately 5' of VirB6, deduplicated to unique sequences, with
the 186 phmmer-confirmed unrecognised VirB5 proteins REMOVED. Leaving them in
would let a block of 237 aa VirB5 dominate the family structure.

Method choices, each already forced by a control:

  phmmer, not MMseqs2   MMseqs2 returned zero cross-hits on the four seeds at the
                        proposed parameters. phmmer recovers both families.
  Leiden, not connected components
                        single linkage percolates: at d=0.007 one component held
                        53% of the database. Connected components on a sequence
                        graph does the same thing.
  E-value, not bit score
                        bit score scales with length, and these proteins span
                        40-300 aa.

SEED CONTROL, and it gates everything: the four controls are two families --
TrbK_RP4 (RP4/R751 trbK) and Eex_IncN (pKM101 eex / R388 candidate). A usable
threshold must put each pair together AND keep the two pairs apart. Their
separation (within-family 1.3e-09 and 1.1e-11; across-family 1.9e-04) was
measured on a 4,410-sequence background; the background here is different, so it
is re-established rather than assumed.

Three branches are reported for the candidate pool:
  joins_known_family   clusters with an eex_occupied member -> divergent member
  novel_family         >=2 candidates, no eex member        -> candidate family
  singleton            alone                                -> no family structure
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CL = "/global/scratch/users/kh36969/exclusion_gene/cluster"
THRESHOLDS = [1e-03, 1e-05, 1e-07, 1e-09]
SEED_FAMS = {"SEED_RP4_trbK": "TrbK_RP4", "SEED_R751_trbK": "TrbK_RP4",
             "SEED_pKM101_eex": "Eex_IncN", "SEED_R388_cand": "Eex_IncN"}


def parse_hdr(h):
    f = h.split("|")
    d = {"id": f[0], "verdict": f[1] if len(f) > 1 else "?"}
    for x in f[2:]:
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import igraph
    import pyhmmer
    from Bio import SeqIO

    recs = list(SeqIO.parse(os.path.join(CL, "cluster_input.faa"), "fasta"))
    meta = {}
    for r in recs:
        d = parse_hdr(r.description)
        meta[d["id"]] = d
    names = [parse_hdr(r.description)["id"] for r in recs]
    pos = {n: i for i, n in enumerate(names)}
    print("sequences: %d" % len(recs))
    vc = collections.Counter(meta[n]["verdict"] for n in names)
    for k, v in vc.most_common():
        print("  %-16s %d" % (k, v))

    alpha = pyhmmer.easel.Alphabet.amino()
    seqs = [pyhmmer.easel.TextSequence(name=n.encode(), sequence=str(r.seq)).digitize(alpha)
            for n, r in zip(names, recs)]

    print("\nphmmer all-vs-all ...")
    edges = {}
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    for top in pyhmmer.phmmer(seqs, seqs, cpus=a.cpus, E=1.0):
        q = dec(top.query.name)
        for h in top:
            t = dec(h.name)
            if t == q:
                continue
            k = (q, t) if q < t else (t, q)
            if k not in edges or h.evalue < edges[k]:
                edges[k] = h.evalue
    print("  edges at E<=1: %d" % len(edges))

    # ---- seed separation, measured on THIS background -----------------------
    print("\nSEED SEPARATION on this input (%d sequences)" % len(recs))
    def ev(x, y):
        k = (x, y) if x < y else (y, x)
        return edges.get(k)
    for x, y, lab in (("SEED_RP4_trbK", "SEED_R751_trbK", "within TrbK_RP4"),
                      ("SEED_pKM101_eex", "SEED_R388_cand", "within Eex_IncN"),
                      ("SEED_RP4_trbK", "SEED_pKM101_eex", "ACROSS families"),
                      ("SEED_R751_trbK", "SEED_R388_cand", "ACROSS families")):
        e = ev(x, y)
        print("  %-16s %-16s %-18s E = %s"
              % (x.replace("SEED_", ""), y.replace("SEED_", ""), lab,
                 ("%.1e" % e) if e else "no edge at E<=1"))

    rows = []
    print("\n%-9s %8s %9s %8s %9s   %-26s %s"
          % ("E", "edges", "clusters", "largest", "singleton", "SEED CONTROL", "seed coverage"))
    print("-" * 118)
    for t in THRESHOLDS:
        el = [(pos[u], pos[v], w) for (u, v), w in edges.items() if w <= t]
        g = igraph.Graph(n=len(names), edges=[(u, v) for u, v, _ in el])
        g.vs["name"] = names
        part = g.community_leiden(objective_function="CPM", resolution=0.05, n_iterations=-1)
        mem = part.membership
        sizes = collections.Counter(mem)
        # seed control
        cl = {s: mem[pos[s]] for s in SEED_FAMS if s in pos}
        trbk = {cl["SEED_RP4_trbK"], cl["SEED_R751_trbK"]}
        eexn = {cl["SEED_pKM101_eex"], cl["SEED_R388_cand"]}
        together = len(trbk) == 1 and len(eexn) == 1
        apart = not (trbk & eexn)
        seedstat = ("PASS" if (together and apart) else
                    "FAIL: " + ("pair split" if not together else "families merged"))
        # seed coverage: how many non-seed sequences share a seed's cluster
        cov = sum(1 for i, m in enumerate(mem)
                  if m in set(cl.values()) and not names[i].startswith("SEED_"))
        rows.append({"E": t, "edges": len(el), "clusters": len(sizes),
                     "largest": max(sizes.values()), "singletons": sum(1 for v in sizes.values() if v == 1),
                     "seed_control": seedstat, "seed_cluster_members": cov,
                     "_mem": mem})
        print("%-9.0e %8d %9d %8d %9d   %-26s %d"
              % (t, len(el), len(sizes), max(sizes.values()),
                 sum(1 for v in sizes.values() if v == 1), seedstat, cov))

    print("\n=== THREE BRANCHES for the candidate pool, per threshold ===")
    print("%-9s %10s %10s %10s   %s"
          % ("E", "joins_known", "novel_fam", "singleton", "novel families (n>=2 cands, no eex)"))
    for r in rows:
        mem = r["_mem"]
        byc = collections.defaultdict(list)
        for i, m in enumerate(mem):
            byc[m].append(names[i])
        j = nf = sg = 0
        nfam = 0
        for m, mem_names in byc.items():
            cands = [x for x in mem_names if meta[x]["verdict"] == "candidate"]
            if not cands:
                continue
            has_eex = any(meta[x]["verdict"] == "eex_occupied" or x.startswith("SEED_")
                          for x in mem_names)
            if has_eex:
                j += len(cands)
            elif len(mem_names) == 1:
                sg += 1
            else:
                nf += len(cands); nfam += 1
        print("%-9.0e %10d %10d %10d   %d" % (r["E"], j, nf, sg, nfam))

    print("\n=== the 39 inverted-only candidates, tracked separately ===")
    inv = [n for n in names if meta[n].get("inv") == "1"]
    print("  tagged inverted-only sequences: %d" % len(inv))
    for r in rows:
        mem = r["_mem"]
        cl = collections.Counter(mem[pos[n]] for n in inv)
        withother = 0
        byc = collections.defaultdict(list)
        for i, m in enumerate(mem):
            byc[m].append(names[i])
        mixed = sum(1 for m in cl if any(not names[pos[x]].startswith("SEED_") and
                                         meta[x].get("inv") != "1" for x in byc[m]))
        witheex = sum(1 for m in cl if any(meta[x]["verdict"] == "eex_occupied"
                                           for x in byc[m]))
        print("  E<=%.0e : %d clusters, largest %d, %d clusters also contain "
              "non-inverted members, %d contain an eex member"
              % (r["E"], len(cl), max(cl.values()), mixed, witheex))

    for r in rows:
        r.pop("_mem")
    dest = os.path.join(PROJ, "data", "anchors", "slot_cluster_scan.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
