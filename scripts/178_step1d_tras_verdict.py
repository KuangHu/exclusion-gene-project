#!/usr/bin/env python3
"""STEP 1d -- the corrected MPF_F verdict. Threshold sweep, ALL TraS families,
full-set matched percentiles, and a permutation null.

Step 1c reported "partial, allele-dependent power" off four TraS families ranked
5 / 97 / 178 / 271. Four defects in that reading, each fixed here:

  1. TraS split into 41 families at 90% id and only the FOUR LARGEST were scored.
     Those four are the BEST case, not the typical one: the other 37 have fewer
     members and lower n_poly, so they will rank worse. The design's post-hoc
     merge needs a large share of them recovered, not one.

  2. 41 may be an artefact of the clustering threshold, not biology. A
     specificity protein under diversifying selection will shatter at 90%. The
     threshold is swept here and the verdict is recomputed at more than one
     value; a conclusion that flips between thresholds is a conclusion about
     mmseqs, not about conjugation.

     CAVEAT ON CALIBRATION, stated because it changes what this sweep can do:
     the design proposed calibrating against the known F / R100 / R64
     specificity groups. That is NOT available. Every plasmid in this TraS set
     was named by PF10624, which detects R100-type ONLY (1 of 14 traS-carrying
     seeds). F-type and R64-type TraS are absent from the set by construction.
     So the 41 families are substructure WITHIN one specificity group, and there
     is no 3-group target to calibrate onto. Threshold ROBUSTNESS is used
     instead: the verdict must survive at every threshold, or it is not a
     verdict.

  3. Matched percentiles were computed inside the top-2000 slice. Low-scoring
     matched families were excluded from their own comparison set, which inflates
     every matched percentile, and inflates them ASYMMETRICALLY -- worst for the
     families with many matched peers. Recomputed here over all families.

  4. A percentile says "better than other families", never "better than chance".
     Under strong phylogenetic structure those differ. Partner allele labels are
     therefore PERMUTED WITHIN each polymorphic cluster (cluster sizes and allele
     frequency spectra preserved), 1000x, giving an empirical p. That separates
     "real but swamped" from "no signal at all" -- the first is a sensitivity
     problem, the second is biology, and rank alone cannot tell them apart.

ALSO REPORTED: for each TraS family, the purity of its TraG allele association.
If the top-ranked family maps onto one clean TraG group while R100-type spreads
across many, the low rank of R100-type is a real finding -- R100 exclusion
specificity is broad -- and not a miss. That distinction decides whether the
method needs fixing or is behaving correctly.

PRE-REGISTERED VERDICT, written before the run:

    >= 50% of TraS families, WEIGHTED BY MEMBER COUNT, at
    matched percentile >= 95 AND permutation p < 0.01

  pass -> build the ICE backbone, proceed to G/B/C
  fail -> the method can VERIFY but not DISCOVER; skip the ICE step

Weighted-half rather than "at least one in the top 10" because on MPF_G there is
no merge prior. Here we know which families to merge -- they are all called TraS
and all sit at traG+1. On MPF_G, a slot held by twenty fragment families of which
the method surfaces two would present as two isolated mid-ranked genes with
nothing to indicate they are one slot. Reconstructing the slot, not hitting one
family, is the bar.
"""
import collections
import csv
import gzip
import math
import os
import random
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
PAIRS = "/global/scratch/users/kh36969/exclusion_gene/mash/pairs_d010.tsv.gz"
WORK = "/global/scratch/users/kh36969/exclusion_gene/recovery/step1d"
STEP1 = "/global/scratch/users/kh36969/exclusion_gene/recovery/step1"
REL = os.path.join(PROJ, "data", "release")
ANCH = os.path.join(PROJ, "data", "anchors")
HMM = os.path.join(PROJ, "data", "hmm")
D_MASH = 0.007
PARTNER_ID = 0.95
SWEEP = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95]   # TraS-only, cheap
SCORE_AT = [0.50, 0.70, 0.90]                  # full re-cluster + scoring
NPERM = 1000
PASS_MATCHED = 95.0
PASS_P = 0.01


def mi(pairs):
    n = len(pairs)
    if n < 2:
        return 0.0
    jx = collections.Counter(pairs)
    mx = collections.Counter(p[0] for p in pairs)
    my = collections.Counter(p[1] for p in pairs)
    out = 0.0
    for (a, b), c in jx.items():
        p = c / n
        out += p * math.log2(p / ((mx[a] / n) * (my[b] / n)))
    return max(out, 0.0)


def cluster(fasta, out, ident):
    r = subprocess.run(["mmseqs", "easy-cluster", fasta, out, out + "_tmp",
                        "--min-seq-id", str(ident), "-c", "0.8", "--cov-mode", "0",
                        "-v", "1"], capture_output=True, text=True)
    if r.returncode:
        raise SystemExit("mmseqs failed at %.2f:\n%s" % (ident, r.stderr[-700:]))
    m = {}
    for line in open(out + "_cluster.tsv"):
        rep, mem = line.rstrip("\n").split("\t")[:2]
        m[mem] = rep
    return m


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    rng = random.Random(0)
    os.makedirs(WORK, exist_ok=True)

    accs = {l.split("\t")[0] for i, l in
            enumerate(open(os.path.join(REL, "v1/catalogue_MPF_F_v1.tsv"))) if i}
    prot, owner, gidx = [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                if f[0] in accs:
                    owner.append(f[0]); gidx.append(int(f[1])); prot.append(line.rstrip("\n"))
    print("MPF_F: %d plasmids, %d proteins" % (len(accs), len(prot)), flush=True)
    pos = {}
    for i in range(len(prot)):
        pos[(owner[i], gidx[i])] = i

    # ---- backbone + partner, as in step 1c ---------------------------------
    par = {a: a for a in accs}

    def find(x):
        r = x
        while par[r] != r:
            r = par[r]
        while par[x] != r:
            par[x], x = r, par[x]
        return r

    with gzip.open(PAIRS, "rt") as fh:
        for line in fh:
            f = line.split("\t", 3)
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
    clu = {a: find(a) for a in accs}
    allele = {}
    for line in open(os.path.join(STEP1, "p_F_%d_cluster.tsv" % int(PARTNER_ID * 100))):
        rep, mem = line.rstrip("\n").split("\t")[:2]
        allele[mem] = rep
    per = collections.defaultdict(set)
    for a in accs:
        if a in allele:
            per[clu[a]].add(allele[a])
    poly = {k for k, v in per.items() if len(v) >= 2}
    keep = [a for a in accs if a in allele and clu[a] in poly]
    kset = set(keep)
    bycl = collections.defaultdict(list)
    for a in keep:
        bycl[clu[a]].append(a)
    N = len(keep)
    print("polymorphic clusters %d | plasmids %d\n" % (len(poly), N), flush=True)

    # ---- the named TraS slot proteins --------------------------------------
    sof = {r["accession"]: r for r in csv.DictReader(
        open(os.path.join(REL, "v1.1/slot_occupant_families.tsv")), delimiter="\t")}
    tras_idx = {}
    # coords -> gene index via cache headers
    want = {a for a, r in sof.items()
            if a in accs and r.get("slot_occupant_family") == "TraS" and r.get("slot_occupant_coords")}
    coord = {a: sof[a]["slot_occupant_coords"] for a in want}
    hdr = {}
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        for line in open(os.path.join(CACHE, fn)):
            if line[0] != ">":
                continue
            f = line[1:].rstrip("\n").split("|")
            if f[0] in want:
                hdr[(f[0], "%s..%s" % (f[2], f[3]))] = int(f[1])
    for a in want:
        gi = hdr.get((a, coord[a]))
        if gi is not None:
            tras_idx[a] = gi
    print("named (R100-type) TraS slot proteins located: %d of %d\n"
          % (len(tras_idx), len(want)), flush=True)

    # ---- PART A: does 41 survive a threshold sweep? -------------------------
    fa_t = os.path.join(WORK, "tras.faa")
    with open(fa_t, "w") as fh:
        for a, gi in sorted(tras_idx.items()):
            fh.write(">%s#%d\n%s\n" % (a, gi, prot[pos[(a, gi)]]))
    print("=== PART A: TraS-only clustering sweep ===")
    print("  (no 3-group target exists: the whole set is R100-type by PF10624)")
    print("  %-8s %10s %14s" % ("min_id", "families", "largest family"))
    for ident in SWEEP:
        m = cluster(fa_t, os.path.join(WORK, "t_%d" % int(ident * 100)), ident)
        cnt = collections.Counter(m.values())
        print("  %-8.2f %10d %14d" % (ident, len(cnt), cnt.most_common(1)[0][1]), flush=True)

    # ---- PART B: full clustering + scoring at several thresholds ------------
    fa_all = os.path.join(WORK, "all.faa")
    with open(fa_all, "w") as fh:
        for i in range(len(prot)):
            fh.write(">%s#%d\n%s\n" % (owner[i], gidx[i], prot[i]))

    import pyhmmer as ph
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    sub = [ph.easel.TextSequence(name=str(i).encode(), sequence=prot[i]).digitize(alpha)
           for i in range(len(prot))]
    with ph.plan7.HMMFile(os.path.join(HMM, "TraG_N.hmm")) as fh:
        hmm = next(iter(fh))
    trag_hits = set()
    for top in ph.hmmsearch([hmm], sub, cpus=16, bit_cutoffs="gathering"):
        for h in top:
            trag_hits.add(int(dec(h.name)))
    print("\nTraG proteins: %d\n" % len(trag_hits), flush=True)

    allrows, verdicts = [], []
    for ident in SCORE_AT:
        print("=" * 76)
        print("SCORING AT min_seq_id = %.2f" % ident)
        print("=" * 76, flush=True)
        fam = cluster(fa_all, os.path.join(WORK, "all_%d" % int(ident * 100)), ident)
        members = collections.defaultdict(set)
        for nm, rep in fam.items():
            members[rep].add(nm.split("#")[0])
        trag_fams = {fam["%s#%d" % (owner[i], gidx[i])] for i in trag_hits
                     if "%s#%d" % (owner[i], gidx[i]) in fam}

        tras_fams = collections.Counter()
        for a, gi in tras_idx.items():
            f = fam.get("%s#%d" % (a, gi))
            if f:
                tras_fams[f] += 1
        print("  gene families %d | TraS families %d | TraG families excluded %d"
              % (len(members), len(tras_fams), len(trag_fams)), flush=True)

        def cond_mi(present, lab=None):
            tot = 0.0
            for c, els in bycl.items():
                obs = [(1 if a in present else 0, (lab or allele)[a]) for a in els]
                if len(set(x for x, _ in obs)) < 2:
                    continue
                tot += (len(els) / N) * mi(obs)
            return tot

        scores = {}
        for f, mem in members.items():
            if f in trag_fams:
                continue
            p = mem & kset
            if not (2 <= len(p) <= N - 2):
                continue
            scores[f] = (cond_mi(p), len(p))
        print("  scored %d families\n" % len(scores), flush=True)

        # permutation null, shared across families
        perms = []
        for _ in range(NPERM):
            lab = {}
            for c, els in bycl.items():
                vals = [allele[a] for a in els]
                rng.shuffle(vals)
                for a, v in zip(els, vals):
                    lab[a] = v
            perms.append(lab)

        print("  %-24s %7s %8s %8s %9s %8s %9s %8s"
              % ("TraS family", "members", "n_poly", "cond_MI", "rank", "glob%",
                 "matched%", "perm_p"))
        tot_mem = sum(tras_fams.values())
        passed_mem = 0
        for f, nmem in tras_fams.most_common():
            if f not in scores:
                print("  %-24s %7d %8s  not scored (constant or partner family)"
                      % (f[:24], nmem, "-"))
                continue
            s, npoly = scores[f]
            rank = 1 + sum(1 for v, _ in scores.values() if v > s)
            gp = 100.0 * (1 - rank / len(scores))
            prev = npoly / N
            lo, hi = prev * 0.8, prev * 1.2
            matched = [v for v, n2 in scores.values() if lo <= n2 / N <= hi]
            mp = 100.0 * (1 - sum(1 for v in matched if v > s) / len(matched)) if matched else float("nan")
            memset = members[f] & kset
            null = [cond_mi(memset, lab) for lab in perms]
            pv = (1 + sum(1 for v in null if v >= s)) / (NPERM + 1)
            ok = (mp >= PASS_MATCHED) and (pv < PASS_P)
            if ok:
                passed_mem += nmem
            print("  %-24s %7d %8d %8.5f %9d %7.2f %8.2f %8.4f %s"
                  % (f[:24], nmem, npoly, s, rank, gp, mp, pv, "PASS" if ok else ""),
                  flush=True)
            # TraG purity
            tg = collections.Counter(allele[a] for a in memset)
            pur = 100.0 * tg.most_common(1)[0][1] / len(memset) if memset else 0.0
            allrows.append({"min_seq_id": ident, "family": f, "members": nmem,
                            "n_poly": npoly, "cond_mi": round(s, 6), "rank": rank,
                            "global_pct": round(gp, 2), "matched_pct": round(mp, 2),
                            "perm_p": round(pv, 5), "n_matched": len(matched),
                            "trag_purity_pct": round(pur, 1),
                            "n_trag_alleles": len(tg), "passes": int(ok)})
        frac = 100.0 * passed_mem / tot_mem if tot_mem else 0.0
        print("\n  member-weighted share passing (matched>=%.0f AND p<%.2f): %.1f%%"
              % (PASS_MATCHED, PASS_P, frac))
        verdicts.append((ident, frac))
        print()

    print("=" * 76)
    print("PRE-REGISTERED VERDICT: >=50%% of TraS families by member weight")
    print("=" * 76)
    for ident, frac in verdicts:
        print("  min_seq_id %.2f : %5.1f%%  %s" % (ident, frac, "PASS" if frac >= 50 else "FAIL"))
    allpass = all(f >= 50 for _, f in verdicts)
    anypass = any(f >= 50 for _, f in verdicts)
    print()
    if allpass:
        print("  PASS at every threshold. The method DISCOVERS, not merely verifies.")
        print("  Build the ICE backbone and proceed to MPF_G/B/C.")
    elif anypass:
        print("  THRESHOLD-DEPENDENT. The verdict flips across mmseqs settings, so it")
        print("  is a statement about clustering, not about conjugation. Do not")
        print("  proceed on it; the threshold must be resolved on external grounds.")
    else:
        print("  FAIL at every threshold. The method can VERIFY a known EEx but")
        print("  cannot RECONSTRUCT a slot from fragments. On MPF_G/B/C, where no")
        print("  merge prior exists, its output would not be interpretable.")
        print("  Skip the ICE backbone. Spend the effort on the structural layer.")

    dest = os.path.join(ANCH, "step1d_tras_verdict.tsv")
    if allrows:
        with open(dest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(allrows)
        print("\nwrote %s (%d rows)" % (dest, len(allrows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
