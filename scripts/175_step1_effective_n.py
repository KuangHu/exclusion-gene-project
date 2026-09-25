#!/usr/bin/env python3
"""STEP 1 -- de-duplicate without destroying the recombinants, and report effective N.

PLSDB holds hundreds of near-identical clinical plasmids. Left in, they fill one
cell of the contingency table on their own and the mutual information reports
sequencing bias as biology.

But the naive fix breaks the experiment. Taking one representative per Mash
cluster discards exactly the informative elements: the signal lives in plasmids
with the SAME backbone and a DIFFERENT partner allele, and whole-plasmid distance
merges precisely those.

    de-duplication unit = (backbone cluster  x  partner allele)

Same backbone, different partner -> both kept. Same backbone, same partner ->
one. This is the only unit under which a recombinant survives de-duplication.

BACKBONE, by dataset:
  plasmids  Mash single-linkage at d <= 0.007, from pairs_d010.tsv.gz. That is
            DEDUPLICATION ONLY, not independence -- the graph percolates at
            d = 0.016-0.018, so no biologically meaningful threshold exists.
  ICEs      no Mash sketch exists for these. VirB4 allele at 95% id is used as a
            backbone proxy: it is the one gene mandatory in all eight class
            definitions. It is a PROXY and is labelled as one; it will
            under-merge relative to true Mash clustering.

THE PARTNER THRESHOLD IS NOT CHOSEN HERE. Step 2 calibrates it against the known
TraG specificity groups (F / R100 / R64). Step 1 therefore SWEEPS it and prints
effective N at every value, so the go/no-go decision can be read off without
depending on a number that has not been calibrated yet.

WHAT THIS STEP IS FOR: one number per class. If a class has only tens of
independent observations, no amount of downstream statistics rescues it and
MPF_G/B/C are not worth running. That decision is meant to be made here, before
the expensive work, not after.

Classes with no established partner (G, C, FATA) get backbone counts only. Naming
a partner for them by guesswork would smuggle in the answer that Step 6 exists to
discover without supervision.
"""
import collections
import csv
import gzip
import os
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
PAIRS = "/global/scratch/users/kh36969/exclusion_gene/mash/pairs_d010.tsv.gz"
CJ = os.path.join(PROJ, "data", "CONJScan", "profiles")
HMM = os.path.join(PROJ, "data", "hmm")
REL = os.path.join(PROJ, "data", "release")
ANCH = os.path.join(PROJ, "data", "anchors")
WORK = "/global/scratch/users/kh36969/exclusion_gene/recovery/step1"
D_MASH = 0.007
IDS = [0.50, 0.70, 0.80, 0.90, 0.95]

# class -> (partner profile, where it lives, why this gene)
PARTNER = {
    "T":    ("TrbL",          HMM, "VirB6 homologue; TrbK's neighbour in RP4"),
    "F":    ("TraG_N",        HMM, "TraG; specificity groups F/R100/R64 are known"),
    "I":    ("T4SS_I_traY",   CJ,  "TraY; donor TraY is ExcA's target (Sakuma 2013)"),
    "FA":   ("FA_orf15",      CJ,  "marks ConG (script 136)"),
    "B":    ("T4SS_B_traJ",   CJ,  "TraJ_B = VirB6 equivalent (Guglielmini)"),
    # G, C, FATA: partner not established -- deliberately absent
}


def load_classes():
    """accession sets per class, plasmids and ICEs, from what is already frozen."""
    pls, ice = {}, {}
    cat = [("T", "v1.1/catalogue_MPF_T_v1.1.tsv"), ("F", "v1/catalogue_MPF_F_v1.tsv"),
           ("I", "v1.1/catalogue_MPF_I_v1.tsv")]
    for c, p in cat:
        pls[c] = {l.split("\t")[0] for i, l in enumerate(open(os.path.join(REL, p))) if i}
    fa = list(csv.DictReader(open(os.path.join(REL, "v1.1/catalogue_MPF_FA_v1.tsv")),
                             delimiter="\t"))
    pls["FA"] = {r["accession"] for r in fa if r["fa_subclass"] == "FA"}
    pls["FATA"] = {r["accession"] for r in fa if r["fa_subclass"] == "FATA"}

    # new classes: thresholds READ OFF the step-2 tables, per class, not one number
    NEW = {"G": 3, "B": 5, "C": 2, "FATA": 2}
    for ds, dest in (("plsdb", pls), ("ice", ice)):
        p = os.path.join(ANCH, "profile_counts_%s.tsv" % ds)
        if not os.path.exists(p):
            continue
        rows = list(csv.DictReader(open(p), delimiter="\t"))
        for c, n in NEW.items():
            s = {r["element"] for r in rows if int(r["n_MPF_" + c]) >= n}
            if ds == "plsdb":
                dest[c] = dest.get(c, set()) | s
            else:
                dest[c] = s
        for c in ("T", "F", "I", "FA"):
            s = {r["element"] for r in rows if int(r["n_MPF_" + c]) >= 2}
            if ds == "ice":
                dest[c] = s
    return pls, ice, NEW


def mash_clusters(accs):
    """Single linkage at d <= 0.007. Returns accession -> cluster id."""
    par = {a: a for a in accs}

    def find(x):
        r = x
        while par[r] != r:
            r = par[r]
        while par[x] != r:
            par[x], x = r, par[x]
        return r

    kept = 0
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
            a, b = f[0], f[1]
            if a in par and b in par:
                ra, rb = find(a), find(b)
                if ra != rb:
                    par[ra] = rb
                    kept += 1
    return {a: find(a) for a in accs}, kept


def best_hit(profile_dir, profile, seqs, owner):
    """element -> its best-scoring protein for this profile."""
    import pyhmmer
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    path = os.path.join(profile_dir, profile + ".hmm")
    if not os.path.exists(path):
        return None
    with pyhmmer.plan7.HMMFile(path) as fh:
        hmm = next(iter(fh))
    best = {}
    for top in pyhmmer.hmmsearch([hmm], seqs, cpus=8, bit_cutoffs="gathering"):
        for h in top:
            i = int(dec(h.name))
            e = owner[i]
            if e not in best or h.score > best[e][0]:
                best[e] = (h.score, i)
    return best


def mmseqs_clusters(fasta, out, ident):
    tmp = out + "_tmp"
    r = subprocess.run(["mmseqs", "easy-cluster", fasta, out, tmp,
                        "--min-seq-id", str(ident), "-c", "0.8", "--cov-mode", "0",
                        "-v", "1"], capture_output=True, text=True)
    if r.returncode:
        raise SystemExit("mmseqs failed (%s):\n%s" % (ident, r.stderr[-800:]))
    m = {}
    for line in open(out + "_cluster.tsv"):
        rep, mem = line.rstrip("\n").split("\t")[:2]
        m[mem] = rep
    return m


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    os.makedirs(WORK, exist_ok=True)
    pls, ice, NEW = load_classes()

    print("=== class membership ===")
    print("  %-6s %9s %9s   %s" % ("class", "plasmids", "ICEs", "partner"))
    for c in ("T", "F", "I", "FA", "FATA", "G", "B", "C"):
        pr = PARTNER.get(c)
        print("  MPF_%-3s %9d %9d   %s"
              % (c, len(pls.get(c, ())), len(ice.get(c, ())),
                 pr[0] if pr else "-- not established --"))
    print("\n  new-class thresholds used: %s" % NEW)

    # ---- proteins for every element we care about, one pass ----------------
    want = set().union(*pls.values()) if pls else set()
    seqs, owner = [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                if nm.split("|")[0] in want:
                    owner.append(nm.split("|")[0])
                    seqs.append(pyhmmer.easel.TextSequence(
                        name=str(len(seqs)).encode(),
                        sequence=line.rstrip("\n")).digitize(alpha))
    print("\nloaded %d proteins over %d plasmids\n" % (len(seqs), len(want)), flush=True)

    rows = []
    for c in ("T", "F", "I", "FA", "FATA", "G", "B", "C"):
        accs = pls.get(c, set())
        if not accs:
            continue
        print("=" * 72)
        print("MPF_%s  --  %d plasmids" % (c, len(accs)))
        clu, edges = mash_clusters(accs)
        n_mash = len(set(clu.values()))
        print("  Mash d<=%.3f : %d clusters (%d merging edges)  collapse %.1fx"
              % (D_MASH, n_mash, edges, len(accs) / n_mash), flush=True)

        pr = PARTNER.get(c)
        if not pr:
            print("  partner not established -- backbone count only.")
            rows.append({"class": "MPF_" + c, "dataset": "plsdb", "n_records": len(accs),
                         "n_backbone": n_mash, "partner": "", "min_seq_id": "",
                         "n_effective": "", "note": "partner not established"})
            continue

        prof, pdir, why = pr
        best = best_hit(pdir, prof, seqs, owner)
        got = {e: v for e, v in (best or {}).items() if e in accs}
        print("  partner %s (%s): found on %d/%d (%.1f%%)"
              % (prof, why, len(got), len(accs), 100.0 * len(got) / len(accs)), flush=True)
        if not got:
            print("  no partner sequences -- cannot form the de-duplication unit.")
            continue

        fa = os.path.join(WORK, "partner_%s.faa" % c)
        with open(fa, "w") as fh:
            for e, (sc, i) in sorted(got.items()):
                fh.write(">%s\n%s\n" % (e, seqs[i].textize().sequence))

        print("  %-9s %10s %12s %12s  %s"
              % ("min_id", "alleles", "effective N", "vs backbone", "reading"))
        for ident in IDS:
            m = mmseqs_clusters(fa, os.path.join(WORK, "p_%s_%d" % (c, int(ident * 100))),
                                ident)
            alle = len(set(m.values()))
            units = {(clu[e], m.get(e, "?")) for e in got}
            neff = len(units)
            print("  %-9.2f %10d %12d %12s  %s"
                  % (ident, alle, neff, "+%d" % (neff - n_mash),
                     "usable" if neff >= 100 else
                     "THIN -- statistics will not carry" if neff >= 30 else
                     "TOO FEW -- do not run the scan"), flush=True)
            rows.append({"class": "MPF_" + c, "dataset": "plsdb",
                         "n_records": len(accs), "n_backbone": n_mash,
                         "partner": prof, "min_seq_id": ident,
                         "n_partner_alleles": alle, "n_effective": neff,
                         "note": ""})
        print()

    dest = os.path.join(ANCH, "step1_effective_n.tsv")
    if rows:
        keys = sorted({k for r in rows for k in r})
        with open(dest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", lineterminator="\n",
                               restval="")
            w.writeheader(); w.writerows(rows)
        print("wrote %s (%d rows)" % (dest, len(rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
