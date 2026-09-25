#!/usr/bin/env python3
"""Toward a high-confidence tra gene cluster catalogue.

Deliberately NOT called an operon catalogue: this is HMM hits and gene positions.
No promoters, terminators or co-transcription evidence has been gathered, and
calling it an operon would foreclose the moron line later (RP4's Rho-independent
terminator between trbJ and trbK; IncC sfx's own constitutive promoter).

Two things, both of which change the denominator of everything downstream.

(1) JOINT COMPLETENESS, which has never been computed. Only marginal failure rates
    exist. It cannot be estimated by assuming independence -- the failures are
    already known not to coincide (VirB6/VirB8/VirB9 intersection 301, union 460),
    so the true 7/7 count is well above 0.95^7 x 7436 ~ 5,200.

(2) SECOND-PASS RECALL by phmmer against confirmed members. The same pattern has
    now appeared four times: PF07996 scores 0/1221 candidates at GA while phmmer
    matches 162/228; IncP TrbJ scores 17.7 and 22.2 against a GA of 24.7; MMseqs2
    returned zero cross-hits on the four seeds; GA systematically drops divergent
    members. If VirB5 falls from 27.1% to single digits, the anchor dilemma
    dissolves and joint completeness rises with it.

    This is `decontaminate.py` run forwards: the same confirmed-member phmmer, used
    to RECOVER anchors rather than to remove them from a nomination set.
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
import decontaminate as DC

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CORE = ["VirD4", "VirB4", "VirB6", "VirB8", "VirB9", "VirB10", "VirB11"]
E_RECALL = 1e-05
MAX_REPS = 25


def main():
    from assertions import require_compute_node
    require_compute_node()          # A11: no heavy scans on a login node
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=20)
    a = ap.parse_args()
    import igraph
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner, idx, aas = [], [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                owner.append(f[0]); idx.append(int(f[1])); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    accs = sorted(set(owner))
    N = len(accs)
    print("cached proteins %d over %d accessions" % (len(seqs), N), flush=True)

    ga = {}
    for name, fam, ev, minaa in DC.ANCHOR_FAMILIES:
        p = os.path.join(HMM, fam + ".hmm")
        if not os.path.exists(p):
            continue
        with pyhmmer.plan7.HMMFile(p) as fh:
            model = next(iter(fh))
        kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
        hits = set()
        prots = set()
        for top in pyhmmer.hmmsearch([model], seqs, cpus=a.cpus, **kw):
            for h in top:
                i = int(dec(h.name))
                if minaa and len(aas[i]) <= minaa:
                    continue
                hits.add(owner[i]); prots.add(aas[i])
        ga[name] = (hits, sorted(prots))
        print("  GA %-8s %5d plasmids (%.1f%% failure), %d unique proteins"
              % (name, len(hits), 100.0 * (N - len(hits)) / N, len(prots)), flush=True)

    def completeness(hitmap, label):
        per = collections.Counter()
        for acc in accs:
            per[sum(1 for k in CORE if acc in hitmap[k])] += 1
        print("\n  JOINT COMPLETENESS over the %d core anchors (%s)" % (len(CORE), label))
        cum = 0
        for k in range(len(CORE), -1, -1):
            cum += per.get(k, 0)
            print("    %d/%d  %5d plasmids (%5.1f%%)   cumulative >=%d: %5d (%.1f%%)"
                  % (k, len(CORE), per.get(k, 0), 100.0 * per.get(k, 0) / N,
                     k, cum, 100.0 * cum / N))
        return per
    before = completeness({k: v[0] for k, v in ga.items()}, "GA only")

    print("\n=== SECOND-PASS RECALL: phmmer vs confirmed members ===", flush=True)
    rows = []
    recovered = {}
    for name, fam, ev, minaa in DC.ANCHOR_FAMILIES:
        if name not in ga:
            continue
        hits, prots = ga[name]
        if len(prots) < 2:
            recovered[name] = set(hits); continue
        # representatives spanning the family's own sequence diversity
        dig = lambda L: [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
                         for i, s in enumerate(L)]
        du = dig(prots)
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
        g = igraph.Graph(n=len(prots), edges=[k for k, w in edges.items() if w <= 1e-5])
        mem = g.community_leiden(objective_function="CPM", resolution=0.05,
                                 n_iterations=10).membership
        byc = collections.defaultdict(list)
        for i, m in enumerate(mem):
            byc[m].append(i)
        reps = []
        for v in sorted(byc.values(), key=len, reverse=True):
            reps.append(max(v, key=lambda i: len(prots[i])))
            if len(reps) >= MAX_REPS:
                break
        new = set()
        for top in pyhmmer.phmmer(dig([prots[i] for i in reps]), seqs,
                                  cpus=a.cpus, E=E_RECALL):
            for h in top:
                i = int(dec(h.name))
                if minaa and len(aas[i]) <= minaa:
                    continue
                if owner[i] not in hits:
                    new.add(owner[i])
        recovered[name] = hits | new
        rows.append({"anchor": name, "ga_plasmids": len(hits),
                     "ga_failure_pct": round(100.0 * (N - len(hits)) / N, 2),
                     "reps_used": len(reps), "recovered_new": len(new),
                     "phmmer_plasmids": len(hits) + len(new),
                     "phmmer_failure_pct": round(100.0 * (N - len(hits) - len(new)) / N, 2)})
        print("  %-8s GA %5d (%5.1f%% fail) + %4d recovered -> %5d (%5.1f%% fail)"
              % (name, len(hits), 100.0 * (N - len(hits)) / N, len(new),
                 len(hits) + len(new), 100.0 * (N - len(hits) - len(new)) / N), flush=True)
    after = completeness(recovered, "GA + phmmer recall")

    print("\n  SHIFT in cumulative >=6/7: %.1f%% -> %.1f%%"
          % (100.0 * sum(before.get(k, 0) for k in (6, 7)) / N,
             100.0 * sum(after.get(k, 0) for k in (6, 7)) / N))
    dest = os.path.join(PROJ, "data", "anchors", "anchor_second_pass_recall.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
