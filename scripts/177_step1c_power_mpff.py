#!/usr/bin/env python3
"""STEP 1c -- POWER TEST on MPF_F. Does 146 polymorphic clusters detect a known EEx?

Step 1b gave the honest sample size: 146 independent recombination contexts for
MPF_F, 89 for MPF_T, 36 for MPF_I. The tempting next move is to enlarge the sample
(Mash-sketch the ICEs). That is the wrong order. Enlarging the sample only helps if
the effect is detectable in principle at this scale -- and that costs minutes to
check, against hours to build an ICE backbone.

    If TraS lands at the 99.9th percentile of the null, 146 is enough.
    If TraS lands at the 85th, 1,147 ICEs will not rescue it, because the
    problem is the effect size under this statistic, not the sample size.

THE STATISTIC. Backbone is the confound: partner allele tracks lineage, so does
everything else on the plasmid, and plain MI ranks the whole backbone first. The
analysis is therefore CONDITIONED on backbone cluster and restricted to the
polymorphic ones, where backbone is held constant and only the partner varies:

    I(gene ; partner allele | backbone cluster)

summed over clusters and weighted by cluster size. Within a cluster, ancestry is
controlled by construction, so what is left is co-occurrence.

THE NULL IS EVERY OTHER GENE FAMILY on the same plasmids, scored identically. No
threshold is invented: TraS's rank in that distribution IS the readout.

  Conditional MI is positively biased for families at intermediate frequency, so
  a global quantile alone would flatter any gene that happens to be present in
  about half the set. TraS is therefore also ranked among FREQUENCY-MATCHED
  families (within +/-20% prevalence). Both numbers are printed; the matched one
  is the one to believe.

  The partner's OWN family is excluded. It has perfect dependence with its own
  allele partition by definition and would otherwise sit at rank 1 and make the
  scan look like it worked.

WHAT ANCHORS "TraS" HERE. PF10624 detects R100-type TraS only (1 of 14 traS-
carrying seeds; F/R64 TraS share no detectable similarity with it). So the
positive control is the FAMILY CONTAINING R100's TraS protein, located by its
catalogue coordinates -- NC_002134.1 gene_index 92 -- not by the model and not by
the gene name. The families occupying the TraG+1 slot are reported alongside it,
since that is where F-type and R64-type TraS would live if they are present.
"""
import collections
import csv
import gzip
import math
import os
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
PAIRS = "/global/scratch/users/kh36969/exclusion_gene/mash/pairs_d010.tsv.gz"
WORK = "/global/scratch/users/kh36969/exclusion_gene/recovery/step1c"
STEP1 = "/global/scratch/users/kh36969/exclusion_gene/recovery/step1"
REL = os.path.join(PROJ, "data", "release")
ANCH = os.path.join(PROJ, "data", "anchors")
HMM = os.path.join(PROJ, "data", "hmm")
D_MASH = 0.007
PARTNER_ID = 0.95          # swept in step 1; 0.95 gave the most contexts (146)
FAMILY_ID = 0.90           # tight layer -- the TraS lesson says do not over-merge
R100 = ("NC_002134.1", 92)   # traS, 159 aa, located in step 0
R100_TRAG = "NC_002134.1"


def mi(pairs):
    """Mutual information of a list of (x, y) observations, in bits."""
    n = len(pairs)
    if n < 2:
        return 0.0
    jx = collections.Counter(pairs)
    mx = collections.Counter(p[0] for p in pairs)
    my = collections.Counter(p[1] for p in pairs)
    out = 0.0
    for (a, b), c in jx.items():
        pxy = c / n
        out += pxy * math.log2(pxy / ((mx[a] / n) * (my[b] / n)))
    return max(out, 0.0)


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    os.makedirs(WORK, exist_ok=True)

    accs = {l.split("\t")[0] for i, l in
            enumerate(open(os.path.join(REL, "v1/catalogue_MPF_F_v1.tsv"))) if i}
    print("MPF_F plasmids: %d" % len(accs), flush=True)

    # ---- proteins -----------------------------------------------------------
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
                    owner.append(f[0]); gidx.append(int(f[1]))
                    prot.append(line.rstrip("\n"))
    print("proteins: %d\n" % len(prot), flush=True)

    # ---- backbone clusters --------------------------------------------------
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
    print("backbone clusters: %d" % len(set(clu.values())), flush=True)

    # ---- partner alleles (reuse step 1's clustering) ------------------------
    pfile = os.path.join(STEP1, "p_F_%d_cluster.tsv" % int(PARTNER_ID * 100))
    if not os.path.exists(pfile):
        print("missing %s -- run step 1 first" % pfile); return 1
    allele = {}
    for line in open(pfile):
        rep, mem = line.rstrip("\n").split("\t")[:2]
        allele[mem] = rep
    per = collections.defaultdict(set)
    for a in accs:
        if a in allele:
            per[clu[a]].add(allele[a])
    poly = {k for k, v in per.items() if len(v) >= 2}
    keep = [a for a in accs if a in allele and clu[a] in poly]
    print("polymorphic clusters: %d | plasmids inside them: %d\n"
          % (len(poly), len(keep)), flush=True)
    if len(poly) < 10:
        print("too few contexts to test. stopping."); return 1

    # ---- gene families ------------------------------------------------------
    # Families are built over ALL MPF_F proteins, not only those on polymorphic
    # clusters. The first attempt clustered the polymorphic subset and the
    # positive control vanished: R100's own backbone cluster is MONOMORPHIC
    # (8 plasmids, 1 TraG allele), so its TraS protein was never assigned a
    # family and the run stopped at its own assertion. Family membership must be
    # defined globally; only the SCORING is restricted to polymorphic clusters.
    kset = set(keep)
    fa = os.path.join(WORK, "mpff.faa")
    ids = []
    with open(fa, "w") as fh:
        for i, p in enumerate(prot):
            nm = "%s#%d" % (owner[i], gidx[i])
            ids.append(nm)
            fh.write(">%s\n%s\n" % (nm, p))
    print("clustering %d proteins at %.2f ..." % (len(ids), FAMILY_ID), flush=True)
    out = os.path.join(WORK, "fam")
    r = subprocess.run(["mmseqs", "easy-cluster", fa, out, out + "_tmp",
                        "--min-seq-id", str(FAMILY_ID), "-c", "0.8",
                        "--cov-mode", "0", "-v", "1"], capture_output=True, text=True)
    if r.returncode:
        print("mmseqs failed:\n%s" % r.stderr[-800:]); return 1
    fam = {}
    for line in open(out + "_cluster.tsv"):
        rep, mem = line.rstrip("\n").split("\t")[:2]
        fam[mem] = rep
    print("gene families: %d\n" % len(set(fam.values())), flush=True)

    # presence/absence per plasmid
    members = collections.defaultdict(set)
    for nm, rep in fam.items():
        members[rep].add(nm.split("#")[0])

    # ---- locate the positive control ---------------------------------------
    key = "%s#%d" % R100
    pos_fam = fam.get(key)
    if pos_fam is None:
        print("*** R100 traS (%s) is not in the clustered set." % key)
        print("    Either R100 is not inside a polymorphic cluster, or the")
        print("    coordinates moved. Without the positive control this test")
        print("    cannot be read. Stopping.")
        return 1
    partner_fam = fam.get("%s#%d" % (R100_TRAG, -1))   # not used; partner excluded below
    npos_in = len(members[pos_fam] & kset)
    print("positive control: R100 traS -> family %s" % pos_fam[:34])
    print("  family size overall            : %d plasmids" % len(members[pos_fam]))
    print("  inside polymorphic clusters    : %d" % npos_in)
    if npos_in < 20:
        print("\n*** Only %d control plasmids fall inside the %d polymorphic"
              % (npos_in, len(poly)))
        print("    clusters. The positive control cannot be scored. Stopping.")
        return 1

    # exclude the partner's own family: any family whose members ARE the partner
    partner_prots = set(allele)          # accessions, not proteins
    # identify by TraG_N hits instead
    import pyhmmer as ph
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    sub = [ph.easel.TextSequence(name=str(i).encode(), sequence=prot[i]).digitize(alpha)
           for i in range(len(prot))]
    submap = list(range(len(prot)))
    with ph.plan7.HMMFile(os.path.join(HMM, "TraG_N.hmm")) as fh:
        hmm = next(iter(fh))
    trag_fams = set()
    for top in ph.hmmsearch([hmm], sub, cpus=16, bit_cutoffs="gathering"):
        for h in top:
            i = submap[int(dec(h.name))]
            f = fam.get("%s#%d" % (owner[i], gidx[i]))
            if f:
                trag_fams.add(f)
    print("partner (TraG) families excluded: %d\n" % len(trag_fams), flush=True)

    # ---- conditional MI over polymorphic clusters ---------------------------
    bycl = collections.defaultdict(list)
    for a in keep:
        bycl[clu[a]].append(a)
    N = len(keep)

    def cond_mi(present):
        tot = 0.0
        for c, els in bycl.items():
            obs = [(1 if a in present else 0, allele[a]) for a in els]
            if len(set(x for x, _ in obs)) < 2:
                continue          # family constant in this cluster: no information
            tot += (len(els) / N) * mi(obs)
        return tot

    scores = []
    for f, mem in members.items():
        if f in trag_fams:
            continue
        p = mem & kset
        if not (2 <= len(p) <= N - 2):
            continue              # constant across the whole set
        scores.append((cond_mi(p), f, len(p)))
    scores.sort(reverse=True)
    print("scored %d families\n" % len(scores), flush=True)

    rank = next((i for i, (s, f, n) in enumerate(scores) if f == pos_fam), None)
    if rank is None:
        print("positive control was filtered out (constant or partner). stopping.")
        return 1
    sc, _, npos = scores[rank]
    pct = 100.0 * (1 - rank / len(scores))
    prev = npos / N

    # frequency-matched null
    lo, hi = prev * 0.8, prev * 1.2
    matched = [s for s, f, n in scores if lo <= n / N <= hi]
    mrank = sum(1 for s in matched if s > sc)
    mpct = 100.0 * (1 - mrank / len(matched)) if matched else float("nan")

    print("=" * 72)
    print("POWER TEST -- MPF_F, %d polymorphic clusters, %d plasmids" % (len(poly), N))
    print("=" * 72)
    print("  positive control (R100-type TraS family)")
    print("    conditional MI      %.5f bits" % sc)
    print("    prevalence          %d / %d  (%.1f%%)" % (npos, N, 100 * prev))
    print("    rank                %d of %d" % (rank + 1, len(scores)))
    print("    global percentile   %.2f" % pct)
    print("    frequency-matched   %.2f  (n=%d families at %.1f-%.1f%% prevalence)"
          % (mpct, len(matched), 100 * lo, 100 * hi))
    print("\n  top 12 families by conditional MI:")
    for i, (s, f, n) in enumerate(scores[:12]):
        tag = "  <== POSITIVE CONTROL" if f == pos_fam else ""
        print("    %2d  %.5f  n=%-5d %s%s" % (i + 1, s, n, f[:34], tag))

    print("\n=== VERDICT ===")
    if mpct >= 99.0:
        print("  TraS is at the %.2f frequency-matched percentile. %d contexts" % (mpct, len(poly)))
        print("  DETECT a known entry-exclusion gene. Enlarging the sample is")
        print("  justified and the method may proceed.")
    elif mpct >= 90.0:
        print("  TraS is at the %.2f percentile -- detectable but not separated." % mpct)
        print("  More contexts may help. Decide on effect size, not on hope.")
    else:
        print("  TraS is at the %.2f frequency-matched percentile." % mpct)
        print("  146 contexts DO NOT detect a known positive. Adding ICEs will not")
        print("  fix this: the limit is the effect size under this statistic, not N.")
        print("  Do NOT build the ICE backbone on the strength of this design.")

    dest = os.path.join(ANCH, "step1c_power_mpff.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["rank", "cond_mi_bits", "n_plasmids", "family", "is_positive_control"])
        for i, (s, f, n) in enumerate(scores[:2000]):
            w.writerow([i + 1, round(s, 6), n, f, int(f == pos_fam)])
    print("\nwrote %s" % dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
