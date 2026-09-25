#!/usr/bin/env python3
"""STEP 1b -- count the recombinants, which is what the module actually spends.

Step 1 reported effective N in the thousands and marked most classes "usable".
That verdict is not safe, and the reason is visible in its own sweep: effective N
barely moves as partner identity tightens. MPF_T goes from 317 TrbL alleles to
1,139 -- a 3.6x split -- and gains 42 units. If subdividing the partner 3.6x
subdivides the backbone hardly at all, then partner allele is nearly a FUNCTION of
backbone cluster, and the two columns of the contingency table are not independent
observations of anything.

The de-duplication unit only buys something where a Mash cluster contains MORE
THAN ONE partner allele. Those elements -- same backbone, different partner -- are
the recombinants, and they are the only observations that carry information about
covariation as opposed to shared ancestry. Everything else is the phylogeny.

So the number that decides this module is not effective N. It is:

    how many backbone clusters hold >= 2 partner alleles,
    and how many elements live in them

Reported per class per identity threshold, alongside the elements involved, so the
go/no-go can be read directly instead of inferred from a total that looks large
because it is dominated by singleton clusters.

Reads the Mash pair file ONCE for all classes (Step 1 read it per class, eight
times over) and reuses the mmseqs clusterings Step 1 already wrote.
"""
import collections
import csv
import gzip
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

PAIRS = "/global/scratch/users/kh36969/exclusion_gene/mash/pairs_d010.tsv.gz"
WORK = "/global/scratch/users/kh36969/exclusion_gene/recovery/step1"
REL = os.path.join(PROJ, "data", "release")
ANCH = os.path.join(PROJ, "data", "anchors")
D_MASH = 0.007
IDS = [0.50, 0.70, 0.80, 0.90, 0.95]
WITH_PARTNER = ["T", "F", "I", "FA", "B"]


def load_classes():
    pls = {}
    for c, p in (("T", "v1.1/catalogue_MPF_T_v1.1.tsv"),
                 ("F", "v1/catalogue_MPF_F_v1.tsv"),
                 ("I", "v1.1/catalogue_MPF_I_v1.tsv")):
        pls[c] = {l.split("\t")[0] for i, l in enumerate(open(os.path.join(REL, p))) if i}
    fa = list(csv.DictReader(open(os.path.join(REL, "v1.1/catalogue_MPF_FA_v1.tsv")),
                             delimiter="\t"))
    pls["FA"] = {r["accession"] for r in fa if r["fa_subclass"] == "FA"}
    rows = list(csv.DictReader(open(os.path.join(ANCH, "profile_counts_plsdb.tsv")),
                               delimiter="\t"))
    pls["B"] = {r["element"] for r in rows if int(r["n_MPF_B"]) >= 5}
    return pls


def main():
    require_compute_node()
    pls = load_classes()
    allacc = set().union(*pls.values())
    print("classes: %s" % ", ".join("MPF_%s %d" % (c, len(pls[c])) for c in WITH_PARTNER))
    print("union of members: %d\n" % len(allacc), flush=True)

    # ---- ONE pass over the pair file -------------------------------------
    par = {a: a for a in allacc}

    def find(x):
        r = x
        while par[r] != r:
            r = par[r]
        while par[x] != r:
            par[x], x = r, par[x]
        return r

    n = 0
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
            a, b = f[0], f[1]
            if a in par and b in par:
                ra, rb = find(a), find(b)
                if ra != rb:
                    par[ra] = rb
                    n += 1
    clu = {a: find(a) for a in allacc}
    print("one pass done: %d merging edges at d<=%.3f\n" % (n, D_MASH), flush=True)

    rows = []
    for c in WITH_PARTNER:
        accs = pls[c]
        cl = {a: clu[a] for a in accs}
        nb = len(set(cl.values()))
        print("=" * 76)
        print("MPF_%-4s  %d plasmids, %d backbone clusters" % (c, len(accs), nb))
        print("  %-7s %9s %11s %12s %11s  %s"
              % ("min_id", "alleles", "poly cl", "elements", "% of class", "max alleles/cl"))
        for ident in IDS:
            p = os.path.join(WORK, "p_%s_%d_cluster.tsv" % (c, int(ident * 100)))
            if not os.path.exists(p):
                print("  %-7.2f  (mmseqs output missing: %s)" % (ident, os.path.basename(p)))
                continue
            allele = {}
            for line in open(p):
                rep, mem = line.rstrip("\n").split("\t")[:2]
                allele[mem] = rep
            per = collections.defaultdict(set)
            members = collections.defaultdict(list)
            for a in accs:
                if a in allele:
                    per[cl[a]].add(allele[a])
                    members[cl[a]].append(a)
            poly = [k for k, v in per.items() if len(v) >= 2]
            nel = sum(len(members[k]) for k in poly)
            mx = max((len(v) for v in per.values()), default=0)
            print("  %-7.2f %9d %11d %12d %10.1f%%  %d"
                  % (ident, len(set(allele.values())), len(poly), nel,
                     100.0 * nel / len(accs), mx), flush=True)
            rows.append({"class": "MPF_" + c, "min_seq_id": ident,
                         "n_backbone_clusters": nb,
                         "n_polymorphic_clusters": len(poly),
                         "n_elements_in_polymorphic": nel,
                         "pct_of_class": round(100.0 * nel / len(accs), 1),
                         "max_alleles_per_cluster": mx})
        print()

    print("=" * 76)
    print("READING")
    print("=" * 76)
    print("""
  'poly cl' = backbone clusters holding >=2 partner alleles.
  'elements' = plasmids inside those clusters.

  These are the ONLY observations that separate covariation from shared ancestry.
  A class whose polymorphic-cluster count is small has no evidence about partner
  swapping no matter how many plasmids it contains, because every other element
  is a repeat of its cluster's single (backbone, partner) combination.

  This does NOT say the biology is absent -- it says PLSDB does not sample it.
  A negative here is a statement about the database, not about conjugation.
""")
    dest = os.path.join(ANCH, "step1b_recombinants.tsv")
    if rows:
        with open(dest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(rows)
        print("wrote %s (%d rows)" % (dest, len(rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
