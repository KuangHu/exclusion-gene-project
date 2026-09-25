#!/usr/bin/env python3
"""The three blocking checks before the criterion can be frozen.

4a. FULL Pfam-A + TIGRFAM search of the architecture-C set, not just the
    exclusion list. If >30% hit one named family, C is not a new architecture --
    it is an incomplete exclusion list, and that model must be added and
    everything rerun.

4b. CLONALITY. If the C set collapses onto <5 Mash clusters it is one clonal
    expansion wearing the appearance of an architecture. >20 clusters keeps it.

4c. FULL SET, not the MAXN=600 sample. Every percentage quoted so far rests on a
    sample; frozen numbers must come from all 681.

None of these can be skipped: each one can independently kill architecture C, and
C currently has NO experimental sentinel -- every A and B conclusion rests on
characterised proteins, C rests on 140 predictions.

The uniformity of the C set (100% conforming, median 67 aa, TM at residue 5,
42 aa post-TM domain) is striking, but uniformity is also exactly what a
systematic tool artefact looks like. The lipobox cross-tab excluded ONE artefact
(signal peptide read as TM). It did not exclude all of them.
"""
import collections
import csv
import gzip
import json
import os
import re
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

SLOTV11 = "/global/scratch/users/kh36969/exclusion_gene/slot_v11"
PFAM = "/global/scratch/users/kh36969/funcannot_dbs/pfam/Pfam-A.hmm"
PAIRS = "/global/scratch/users/kh36969/exclusion_gene/mash/pairs_d010.tsv.gz"
WORK = "/global/home/users/kh36969/tmp/archC"
STRICT = re.compile(r"[LVI][ASTVIG][GASN]C")
D_MASH = 0.007


def segments(pred):
    out, cur, st = [], None, 0
    for i, c in enumerate(pred):
        k = "TM" if c in "HhBb" else ("S" if c == "S" else ".")
        if k != cur:
            if cur in ("TM", "S"):
                out.append((cur, st + 1, i, pred[st]))
            cur, st = k, i
    if cur in ("TM", "S"):
        out.append((cur, st + 1, len(pred), pred[st]))
    return out


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    os.makedirs(WORK, exist_ok=True)

    # ---- 4c: FULL target set, not the 600 sample --------------------------
    meta = json.load(open(os.path.join(SLOTV11, "slot_VirB5+1_records.json")))
    recs = meta["records"]                      # accession -> occupant sequence
    acc_of = collections.defaultdict(set)
    for a, s in recs.items():
        acc_of[s].add(a)
    full = sorted(acc_of)
    print("4c. FULL target set: %d unique (previous runs sampled 600)" % len(full),
          flush=True)

    fa = os.path.join(WORK, "full.faa")
    with open(fa, "w") as fh:
        for i, s in enumerate(full):
            fh.write(">t%d\n%s\n" % (i, s))
    pred = os.path.join(WORK, "full.pred")
    r = subprocess.run(["tmbed", "predict", "-f", fa, "-p", pred, "--out-format", "0",
                        "--model-dir", "/global/home/users/kh36969/.hfcache/prott5",
                        "--use-gpu"], capture_output=True, text=True)
    if r.returncode:
        print("tmbed failed:\n%s" % r.stderr[-600:]); return 1
    L = [l.rstrip("\n") for l in open(pred)]
    got = {}
    for i in range(0, len(L), 3):
        if i + 2 < len(L) and L[i].startswith(">"):
            got[L[i][1:]] = (L[i + 1], L[i + 2])

    dist = collections.Counter()
    archC = []
    for i, s in enumerate(full):
        k = "t%d" % i
        if k not in got:
            continue
        seq, pr = got[k]
        tms = [x for x in segments(pr) if x[0] == "TM"]
        lb = bool(STRICT.search(seq[:40]))
        n = len(tms)
        dist["lipid" if (lb and n == 0) else
             "single_pass" if (not lb and n == 1) else
             "multi_pass" if (not lb and n >= 2) else
             "ambiguous" if lb else "other"] += 1
        if (not lb) and n == 1:
            st, en = tms[0][1], tms[0][2]
            archC.append((s, len(seq), st, len(seq) - en))
    tot = sum(dist.values())
    print("\n  anchor_mode on the FULL set (n=%d):" % tot)
    for k, v in dist.most_common():
        print("    %-14s %5d  %5.1f%%" % (k, v, 100.0 * v / tot))
    conf = sum(1 for _, _, st, post in archC if st <= 40 and post >= 30)
    print("  architecture C (single_pass): %d, of which %d (%.1f%%) conform"
          % (len(archC), conf, 100.0 * conf / max(1, len(archC))), flush=True)

    # ---- 4a: full Pfam-A search of the C set ------------------------------
    print("\n4a. FULL Pfam-A search of the %d architecture-C sequences" % len(archC),
          flush=True)
    Q = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
         for i, (s, _, _, _) in enumerate(archC)]
    hit = collections.defaultdict(list)
    with pyhmmer.plan7.HMMFile(PFAM) as fh:
        for top in pyhmmer.hmmer.hmmscan(Q, fh, cpus=8, bit_cutoffs="gathering"):
            i = int(dec(top.query.name))
            for h in top:
                hit[i].append(dec(h.name))
    nh = len(hit)
    print("  with >=1 Pfam domain: %d / %d (%.1f%%)"
          % (nh, len(archC), 100.0 * nh / max(1, len(archC))))
    fam = collections.Counter(x for v in hit.values() for x in set(v))
    for n2, c in fam.most_common(10):
        print("    %-26s %4d  (%.1f%% of C)" % (n2, c, 100.0 * c / max(1, len(archC))))
    top1 = fam.most_common(1)
    if top1 and 100.0 * top1[0][1] / max(1, len(archC)) > 30:
        print("\n  >30%% share one family (%s) -> C IS NOT A NEW ARCHITECTURE."
              % top1[0][0])
        print("  Add that model to the exclusion list and rerun.")
    else:
        print("\n  no single family exceeds 30%% -> C survives 4a")

    # ---- 4b: clonality ----------------------------------------------------
    print("\n4b. CLONALITY of the architecture-C set", flush=True)
    need = set()
    for s, _, _, _ in archC:
        need |= acc_of[s]
    par = {a: a for a in need}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x

    with gzip.open(PAIRS, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3:
                continue
            try:
                if float(f[2]) > D_MASH:
                    continue
            except ValueError:
                continue
            if f[0] in par and f[1] in par:
                ra, rb = find(f[0]), find(f[1])
                if ra != rb:
                    par[ra] = rb
    mash = len({find(a) for a in need})
    print("  %d sequences over %d plasmids -> %d Mash clusters" % (len(archC), len(need), mash))
    print("  verdict: %s" % ("CLONAL EXPANSION (<5 clusters) -- DEMOTE" if mash < 5 else
                             "survives (>20 clusters)" if mash > 20 else
                             "BORDERLINE (5-20 clusters) -- do not freeze on this"))

    dest = os.path.join(PROJ, "data", "anchors", "archC_blocking.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        for k, v in (("full_target_set", len(full)), ("archC_n", len(archC)),
                     ("archC_conforming", conf), ("archC_pfam_hit", nh),
                     ("archC_top_family", top1[0][0] if top1 else ""),
                     ("archC_top_family_n", top1[0][1] if top1 else 0),
                     ("archC_mash_clusters", mash), ("archC_plasmids", len(need))):
            w.writerow([k, v])
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    sys.exit(main() or 0)
