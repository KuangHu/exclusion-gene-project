#!/usr/bin/env python3
"""Resolution sweep for the slot clustering, gated on the seed control.

The largest cluster held ~21% of the input at the chosen resolution (0.05), and
this project has already been burned once by percolation: single-linkage on the
Mash graph put 53% of the database in one component at 95% ANI. So the resolution
has to be shown not to be doing the same thing.

Cluster count and largest-cluster size are weak evidence for that. The seed
control is not: the four controls are two families, TrbK_RP4 (RP4/R751 trbK) and
Eex_IncN (pKM101 eex / R388 candidate), and a resolution that over-merges will
fuse them. That is ground truth, so it is the primary criterion here and the
size statistics are reported alongside rather than relied on.

Connected components is computed at the same threshold purely as the percolation
comparison -- the method NOT used.
"""
import collections
import csv
import gzip
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CL = "/global/scratch/users/kh36969/exclusion_gene/cluster"
SEEDS = ["SEED_RP4_trbK", "SEED_R751_trbK", "SEED_pKM101_eex", "SEED_R388_cand"]
RESOLUTIONS = [0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
E_PASSING = [1e-03, 1e-05]      # 1e-07 and 1e-09 split the TrbK pair


def hdr(h):
    f = h.split("|")
    d = {"id": f[0], "verdict": f[1]}
    for x in f[2:]:
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
    return d


def main():
    import igraph
    from Bio import SeqIO
    recs = list(SeqIO.parse(os.path.join(CL, "cluster_input.faa"), "fasta"))
    meta, names = {}, []
    for r in recs:
        d = hdr(r.description)
        meta[d["id"]] = d
        names.append(d["id"])
    pos = {n: i for i, n in enumerate(names)}
    edges = []
    with gzip.open(os.path.join(CL, "edges.tsv.gz"), "rt") as fh:
        for line in fh:
            u, v, w = line.rstrip("\n").split("\t")
            edges.append((u, v, float(w)))
    print("sequences %d, cached edges %d" % (len(names), len(edges)))

    rows = []
    for E in E_PASSING:
        el = [(pos[u], pos[v]) for u, v, w in edges if w <= E]
        g = igraph.Graph(n=len(names), edges=el)
        cc = g.connected_components()
        szc = collections.Counter(cc.membership)
        print("\n=== E <= %.0e : %d edges ===" % (E, len(el)))
        print("  connected components (the method NOT used): %d components, "
              "largest %d = %.1f%% of input" %
              (len(szc), max(szc.values()), 100.0 * max(szc.values()) / len(names)))
        print("  %-8s %9s %8s %9s  %-30s %s"
              % ("res", "clusters", "largest", "singleton", "SEED CONTROL", "largest-cluster composition"))
        for res in RESOLUTIONS:
            sys.stdout.flush()
            part = g.community_leiden(objective_function="CPM", resolution=res,
                                      n_iterations=10)
            mem = part.membership
            sz = collections.Counter(mem)
            cl = {s: mem[pos[s]] for s in SEEDS}
            trbk = {cl["SEED_RP4_trbK"], cl["SEED_R751_trbK"]}
            eexn = {cl["SEED_pKM101_eex"], cl["SEED_R388_cand"]}
            if len(trbk) > 1 or len(eexn) > 1:
                stat = "FAIL: a known pair is split"
            elif trbk & eexn:
                stat = "FAIL: the two families MERGED"
            else:
                stat = "PASS"
            big = sz.most_common(1)[0][0]
            comp = collections.Counter(meta[names[i]]["verdict"]
                                       for i, m in enumerate(mem) if m == big)
            print("  %-8s %9d %8d %9d  %-30s %s"
                  % (res, len(sz), sz.most_common(1)[0][1],
                     sum(1 for v in sz.values() if v == 1), stat,
                     dict(comp)))
            # family extension = candidates sharing a cluster with a named member
            byc = collections.defaultdict(list)
            for i, m in enumerate(mem):
                byc[m].append(names[i])
            ext = nf = sg = nfam = 0
            for m, mm in byc.items():
                cands = [x for x in mm if meta[x]["verdict"] == "candidate"]
                if not cands:
                    continue
                named = any(meta[x]["verdict"] == "eex_occupied" or x.startswith("SEED_")
                            for x in mm)
                if named:
                    ext += len(cands)
                elif len(mm) == 1:
                    sg += 1
                else:
                    nf += len(cands); nfam += 1
            rows.append({"E": E, "resolution": res, "clusters": len(sz),
                         "largest": sz.most_common(1)[0][1],
                         "singletons": sum(1 for v in sz.values() if v == 1),
                         "seed_control": stat,
                         "family_extension_candidates": ext,
                         "novel_families": nfam,
                         "candidates_in_novel_families": nf,
                         "candidate_singletons": sg})

    print("\n=== FAMILY EXTENSION -- the informative number ===")
    print("`candidate` is DEFINED as 'no named Eex family hits it', so candidates")
    print("not co-clustering with a named member is a restatement of the definition,")
    print("not a finding. The measurable quantity is how many the clustering RESCUES")
    print("into a known family -- i.e. how far iteration can extend what is known.")
    print("\n%-9s %-8s %-12s %-22s %-14s %s"
          % ("E", "res", "seed", "family_extension", "novel_families", "cand_singletons"))
    for r in rows:
        if not r["seed_control"].startswith("PASS"):
            continue
        print("%-9.0e %-8s %-12s %-22d %-14d %d"
              % (r["E"], r["resolution"], "PASS", r["family_extension_candidates"],
                 r["novel_families"], r["candidate_singletons"]))

    dest = os.path.join(PROJ, "data", "anchors", "slot_cluster_resolution_sweep.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
