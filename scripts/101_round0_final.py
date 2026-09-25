#!/usr/bin/env python3
"""Round 0, final: one admission criterion (lipobox), everything else recorded.

Why a single criterion. Of the axes tried, only lipobox survived every robustness
check:

    family count                anti-correlated -- VirB10 yields 62 vs target 35
    within-cluster presence     dead -- only 8 mixed clusters
    Asp@+2                      collapsed -- Asp is 5.1% of NF033894's own members,
                                a family containing a confirmed Asp-carrying eex
    in-slot / total             null max 0.745 vs real 0.771; margin 0.026
    lipobox                     holds under strict / relaxed / structural
                                definitions, 9-20x, all four seeds 100% under each,
                                and the occupancy correction is invariant to it

Adding a criterion whose null overlaps the signal costs true positives and buys no
precision -- that is measurable here, not a matter of principle.

Decontamination is NOT a second criterion, it is data cleaning, and under a
lipobox-only rule it is mandatory: VirB7/TraN is itself a small lipoprotein whose
detection fails on 97% of plasmids, so every unrecognised VirB7 passes a lipobox
test by construction. The same structural trap has already been paid for three
times (VirB5 in the candidate pool, VirB5 in the recovered set, VirB11 in the
null), because in a colinear operon anchor+1 IS the next anchor.

Two costs, stated rather than hidden:

  * IncI2 is REJECTED. Its dominant 85 aa occupant has no hydrophobic stretch under
    any definition. It is excluded BY THE CRITERION, not absent from the data, and
    belongs in a separate case study rather than the main line.
  * Scope narrows to MPF_T. The lipobox calibration rests on four MPF_T seeds. F
    traS is curated inner-membrane with no recorded lipobox; R64 ExcA is reported
    in both membrane-bound and soluble forms. So the claim is "small lipoproteins
    in the VirB5-adjacent slot of MPF_T plasmids", NOT "exclusion genes are
    lipoproteins". Extending to MPF_F/MPF_I requires recalibration.
"""
import argparse
import collections
import csv
import os
import re
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
import decontaminate as DC

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "TraS", "DUF4467"]
STRICT = re.compile(r"[LVI][ASTVIG][GASN]C")
RELAX = re.compile(r"[LVIAFMT][ASTVIGC][GASNC]C")
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
      'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,
      'Y':-1.3,'V':4.2}
E_CLUSTER, RESOLUTION, FAM_LIPO = 1e-05, 0.05, 0.50


def hreg(s):
    for i, c in enumerate(s[:25]):
        if c == "C" and i >= 8 and sum(KD.get(x, 0) for x in s[i-8:i]) / 8.0 >= 1.0:
            return True
    return False
DEFS = {"strict": lambda s: bool(STRICT.search(s[:40])),
        "relaxed": lambda s: bool(RELAX.search(s[:40])),
        "structural": hreg}


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

    def nominate_target():
        out = {}
        for acc in B6:
            b6 = best(B6, acc); b5 = best(B5, acc); b4 = best(B4, acc)
            t = 1 if b6[1] == 1 else -1
            if b5 and b5[1] == b6[1]:
                j = bypos.get((acc, b5[0] + t))                # +1 from VirB5
            elif b4 and b4[1] == b6[1]:
                arch = "canonical" if b4[0] * t < b6[0] * t else "rearranged"
                j = bypos.get((acc, b6[0] + (-1 if arch == "canonical" else 2) * t))
            else:
                continue
            if j is not None and strand[j] == b6[1]:
                out[acc] = j
        return out

    def nominate_neg():
        out = {}
        for acc in B10:
            b = best(B10, acc)
            j = bypos.get((acc, b[0] + (1 if b[1] == 1 else -1)))
            if j is not None and strand[j] == b[1]:
                out[acc] = j
        return out

    print("\nbuilding confirmed-member sets for all anchor families ...", flush=True)
    members = DC.confirmed_members(seqs, aas, cpus=a.cpus)
    for k, v in sorted(members.items()):
        print("  %-8s %d unique confirmed" % (k, len(v)))

    eexhit = {}
    for fam in EEX:
        for acc, v in scan(fam).items():
            for _, _, _, i in v:
                eexhit[i] = fam

    results = []
    for label, nomfn in (("TARGET VirB5-adjacent", nominate_target),
                         ("NEGATIVE VirB10+1", nominate_neg)):
        nom = nomfn()
        byseq = collections.defaultdict(list)
        for acc, j in nom.items():
            byseq[aas[j]].append(j)
        uniq = list(byseq)
        print("\n=== %s ===" % label, flush=True)
        print("  nominated: %d plasmids, %d unique proteins" % (len(nom), len(uniq)))
        kept, dropped = DC.decontaminate(uniq, members, cpus=a.cpus)
        if len(kept) < 10:
            continue
        du = [pyhmmer.easel.TextSequence(name=str(k).encode(), sequence=s).digitize(alpha)
              for k, s in enumerate(kept)]
        edges = {}
        for top in pyhmmer.phmmer(du, du, cpus=a.cpus, E=1.0):
            q = int(dec(top.query.name))
            for h in top:
                t2 = int(dec(h.name))
                if t2 == q:
                    continue
                k2 = (q, t2) if q < t2 else (t2, q)
                if k2 not in edges or h.evalue < edges[k2]:
                    edges[k2] = h.evalue
        g = igraph.Graph(n=len(kept), edges=[k2 for k2, w in edges.items() if w <= E_CLUSTER])
        mem = g.community_leiden(objective_function="CPM", resolution=RESOLUTION,
                                 n_iterations=10).membership
        byc = collections.defaultdict(list)
        for i, m in enumerate(mem):
            byc[m].append(i)
        fams = [v for v in byc.values() if len(v) >= 2]
        named = {i for i, s in enumerate(kept) if any(j in eexhit for j in byseq[s])}
        row = {"set": label, "nominated_plasmids": len(nom), "unique_nominated": len(uniq),
               "removed_as_anchor": len(uniq) - len(kept), "unique_clean": len(kept),
               "families_ge2": len(fams)}
        print("  families (>=2 members): %d" % len(fams))
        for dname, fn in DEFS.items():
            rich = [v for v in fams
                    if sum(1 for i in v if fn(kept[i])) / len(v) >= FAM_LIPO]
            novel = [v for v in rich if not (set(v) & named)]
            row["lipobox_rich_%s" % dname] = len(rich)
            row["novel_%s" % dname] = len(novel)
            print("    %-11s lipobox-rich families %3d   of which no named Eex member %3d"
                  % (dname, len(rich), len(novel)))
        results.append(row)

    dest = os.path.join(PROJ, "data", "anchors", "round0_final.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(results)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
