#!/usr/bin/env python3
"""Round 0: cluster the decontaminated v1.1 slot occupants.

Input: 681 unique sequences surviving decontamination (job 26048996), plus the
four seed controls.

METHOD, each choice already forced by a control:
  phmmer, not MMseqs2       the k-mer prefilter returns ZERO cross-hits on 69-76 aa
                            proteins; the four seeds recover only 1 of 2 families
  Leiden CPM, not connected components
                            single linkage percolates -- 80.3% of sequences enter
                            one component at 1e-03, 41.2% at 1e-05
  E-value, not bitscore     bitscore scales with length and these span 30-526 aa

SEED CONTROL -- GATES EVERYTHING. A usable (threshold, resolution) must:
    RP4 trbK  <-> R751 trbK    SAME cluster   (true edge E = 6.0e-07)
    pKM101 eex <-> R388 cand   SAME cluster   (true edge E = 4.8e-09)
    TrbK_RP4  vs  Eex_IncN     DIFFERENT clusters
Nothing downstream is read from a configuration that fails this.

1e-07 SPLITS the RP4/R751 pair, whose real edge is 6.0e-07. "Stricter is safer"
is false here, and the sweep is run wide enough to show it rather than assert it.

REPRESENTATIVES: medoid, never `max(v, key=len)`. The longest member is the most
domain-rich and most promiscuous query -- the root cause of the retracted phmmer
recall, whose "recovered" VirB5 proteins had median length 816 against a family
median of 238. Clusters are also reported with their sizes so a 3.7% length-outlier
cluster cannot contribute a representative at equal weight with a 5,200-member one.
"""
import argparse, collections, csv, json, math, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

OUT = "/global/scratch/users/kh36969/exclusion_gene/slot_v11"
SEEDFA = "/global/scratch/users/kh36969/exclusion_gene/mmtest/seeds.faa"
SEED_FAM = {"RP4_trbK_TrbK_RP4family": "TrbK_RP4",
            "R751_trbK_TrbK_RP4family": "TrbK_RP4",
            "pKM101_eex_Eex_IncNfamily": "Eex_IncN",
            "R388_cand_Eex_IncNfamily": "Eex_IncN"}
THRESHOLDS = [1e-03, 1e-05, 1e-07]
RESOLUTIONS = [0.005, 0.01, 0.05, 0.1]


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="VirB5+1")
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer, igraph
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    names, aas = [], []
    fa = os.path.join(OUT, "slot_%s_clean.faa" % a.tag)
    nm = None
    for line in open(fa):
        if line[0] == ">": nm = line[1:].split()[0]
        else: names.append(nm); aas.append(line.rstrip("\n"))
    ncand = len(names)
    nm = None
    for line in open(SEEDFA):
        if line[0] == ">": nm = line[1:].strip()
        else: names.append(nm); aas.append(line.rstrip("\n"))
    print("candidates %d + seeds %d = %d sequences" % (ncand, len(names)-ncand, len(names)),
          flush=True)
    recs = json.load(open(os.path.join(OUT, "slot_%s_records.json" % a.tag)))["records"]
    seq2n = collections.Counter(recs.values())

    seqs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
            for i, s in enumerate(aas)]
    best = {}
    for top in pyhmmer.phmmer(seqs, seqs, cpus=a.cpus, E=10.0):
        i = int(dec(top.query.name))
        for h in top:
            j = int(dec(h.name))
            if i == j: continue
            k = (min(i, j), max(i, j))
            if k not in best or h.evalue < best[k]: best[k] = h.evalue
    print("all-vs-all pairs with E<=10: %d\n" % len(best), flush=True)

    sidx = {n: i for i, n in enumerate(names) if n in SEED_FAM}
    print("=== seed separation on THIS background ===")
    for x, y, lab in (("RP4_trbK_TrbK_RP4family", "R751_trbK_TrbK_RP4family", "within TrbK_RP4"),
                      ("pKM101_eex_Eex_IncNfamily", "R388_cand_Eex_IncNfamily", "within Eex_IncN"),
                      ("RP4_trbK_TrbK_RP4family", "pKM101_eex_Eex_IncNfamily", "ACROSS families")):
        i, j = sidx[x], sidx[y]
        e = best.get((min(i, j), max(i, j)))
        print("  %-18s E = %s" % (lab, ("%.1e" % e) if e else "no edge (E>10)"))

    print("\n=== sweep: seed control must pass before anything is read ===")
    print("  %-9s %-7s %7s %7s %9s  %s" % ("E<=", "res", "clust", "n>=2", "largest", "seed control"))
    results = {}
    for th in THRESHOLDS:
        el = [(u, v) for (u, v), e in best.items() if e <= th]
        g = igraph.Graph(n=len(names), edges=el)
        for res in RESOLUTIONS:
            part = g.community_leiden(objective_function="CPM", resolution=res, n_iterations=-1)
            mem = part.membership
            ok_t = mem[sidx["RP4_trbK_TrbK_RP4family"]] == mem[sidx["R751_trbK_TrbK_RP4family"]]
            ok_e = mem[sidx["pKM101_eex_Eex_IncNfamily"]] == mem[sidx["R388_cand_Eex_IncNfamily"]]
            ok_s = mem[sidx["RP4_trbK_TrbK_RP4family"]] != mem[sidx["pKM101_eex_Eex_IncNfamily"]]
            sizes = collections.Counter(mem)
            verdict = ("PASS" if (ok_t and ok_e and ok_s) else
                       "FAIL:" + ",".join(n for n, o in (("trbK_split", ok_t),
                                                         ("eex_split", ok_e),
                                                         ("families_merged", ok_s)) if not o))
            print("  %-9.0e %-7s %7d %7d %9d  %s"
                  % (th, res, len(sizes), sum(1 for c in sizes.values() if c >= 2),
                     max(sizes.values()), verdict))
            results[(th, res)] = (mem, verdict.startswith("PASS"))

    # ADOPTED: 1e-05 x 0.01.
    # The instruction specified resolution 0.05, which FAILS the seed control on
    # this background (job 26051504): it splits the Eex_IncN pair. The claimed
    # pass band "1e-03 to 1e-05 x 0.005-0.1" is not a rectangle here -- at 1e-05
    # only 0.005 and 0.01 pass. Threshold choice IS load-bearing on this pool.
    #
    # Why the band moved: seed separation is tighter on this background than the
    # documented values (within TrbK_RP4 2.3e-07 vs 6.0e-07; within Eex_IncN
    # 1.8e-09 vs 4.8e-09). The pool is 681 sequences, not the 4,410 the original
    # separation was measured on, and E-values scale with database size -- the
    # same property that invalidated the VirB5 E<=1e-5 anchor threshold (A13).
    #
    # E stays at 1e-05 because Round 1 must search at the clustering threshold;
    # moving to 1e-03 to gain resolution headroom would inflate the Round 1
    # denominator for reasons unrelated to biology.
    th, res = 1e-05, 0.01
    mem, ok = results[(th, res)]
    print("\n=== adopted configuration: E<=%.0e, resolution %s -> %s ==="
          % (th, res, "PASS" if ok else "FAIL"))
    if not ok:
        print("  SEED CONTROL FAILED at the adopted configuration. Nothing is read.")
        return
    cl = collections.defaultdict(list)
    for i, m in enumerate(mem): cl[m].append(i)
    cl = {k: v for k, v in cl.items() if any(i < ncand for i in v) or
          any(names[i] in SEED_FAM for i in v)}
    sim = collections.defaultdict(dict)
    for (u, v), e in best.items():
        if e <= th:
            s = -math.log10(max(e, 1e-300))
            sim[u][v] = s; sim[v][u] = s
    rows = []
    for c, mem_i in sorted(cl.items(), key=lambda x: -len(x[1])):
        cand = [i for i in mem_i if i < ncand]
        if not cand: continue
        nrec = sum(seq2n.get(aas[i], 0) for i in cand)
        med = max(cand, key=lambda i: sum(sim[i].get(j, 0.0) for j in cand))
        sd = [names[i] for i in mem_i if names[i] in SEED_FAM]
        rows.append({"cluster": c, "n_unique": len(cand), "n_records": nrec,
                     "medoid": names[med], "medoid_aa": len(aas[med]),
                     "median_aa": sorted(len(aas[i]) for i in cand)[len(cand)//2],
                     "seeds_in_cluster": ";".join(sd),
                     "branch": ("joins_seed_family" if sd else
                                "novel_family" if len(cand) >= 2 else "singleton"),
                     "medoid_seq": aas[med],
                     "member_seqs": ";".join(aas[i] for i in cand)})
    br = collections.Counter(r["branch"] for r in rows)
    print("  clusters containing candidates: %d" % len(rows))
    print("  branches: %s" % ", ".join("%s %d" % kv for kv in br.most_common()))
    print("\n  %-8s %7s %8s %9s %9s  %s" % ("cluster", "uniq", "records", "medoid_aa",
                                            "median_aa", "branch / seeds"))
    for r in rows[:15]:
        print("  %-8d %7d %8d %9d %9d  %s %s"
              % (r["cluster"], r["n_unique"], r["n_records"], r["medoid_aa"],
                 r["median_aa"], r["branch"], r["seeds_in_cluster"]))
    dest = os.path.join(PROJ, "data", "anchors", "round0_%s_v11.tsv" % a.tag.replace("+", "p"))
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d clusters)" % (dest, len(rows)))
    tot_u = sum(r["n_unique"] for r in rows); tot_r = sum(r["n_records"] for r in rows)
    print("  totals: %d unique / %d records across %d clusters" % (tot_u, tot_r, len(rows)))


if __name__ == "__main__":
    main()
