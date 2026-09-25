#!/usr/bin/env python3
"""Is the MPF_FA residual actually Gram-positive, or is it junk?

Job 25991904 rejected `T4SS_virb4` alone (79.5% cross-class) and produced the
conjunction residual `virb4 AND NOT (T or F or I)` = 4,747 plasmids.

It reported per-profile overlaps -- FATA_cd419a 53.4%, FA_orf14 23.6%, FA_orf13
12.4% ... -- but those OVERLAP EACH OTHER and cannot be added. The question they
do not answer:

    how many of the 4,747 carry ANY FA/FATA profile, and how many carry NONE?

A residual whose members mostly carry no Gram-positive profile at all is not an
MPF_FA container; it is whatever else has a VirB4-family ATPase.

THE CONTROL, without which the union is uninterpretable: the same union rate
measured on MPF_T and MPF_F plasmids. The FA/FATA profiles are only evidence of
Gram-positive identity if they are RARE on known Gram-negative systems. If MPF_T
plasmids carry them at a similar rate, the profiles are not class-specific and
the residual's rate means nothing.

Cheap because it is targeted: one full-cache scan for virb4, then the 34 FA/FATA
profiles over only the residual + two control samples, not the whole 7.7M.
"""
import collections, csv, os, random, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
REL = os.path.join(PROJ, "data", "release")
CATS = [("MPF_T", os.path.join(REL, "v1.1", "catalogue_MPF_T_v1.1.tsv")),
        ("MPF_F", os.path.join(REL, "v1", "catalogue_MPF_F_v1.tsv")),
        ("MPF_I", os.path.join(REL, "v1.1", "catalogue_MPF_I_v1.tsv"))]
NCTRL = 2000


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    admitted = {}
    for name, path in CATS:
        admitted[name] = {l.split("\t")[0] for i, l in enumerate(open(path)) if i}
    allprev = set().union(*admitted.values())

    seqs, owner = [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                owner.append(nm.split("|")[0])
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    print("cache: %d proteins over %d plasmids" % (len(seqs), len(set(owner))), flush=True)

    with pyhmmer.plan7.HMMFile(os.path.join(CJ, "T4SS_virb4.hmm")) as fh:
        m = next(iter(fh))
    vb4 = set()
    for top in pyhmmer.hmmsearch([m], seqs, cpus=16, bit_cutoffs="gathering"):
        for h in top: vb4.add(owner[int(dec(h.name))])
    res = vb4 - allprev
    print("virb4 %d ; residual %d" % (len(vb4), len(res)), flush=True)

    rng = random.Random(0)
    ct = set(rng.sample(sorted(admitted["MPF_T"]), min(NCTRL, len(admitted["MPF_T"]))))
    cf = set(rng.sample(sorted(admitted["MPF_F"]), min(NCTRL, len(admitted["MPF_F"]))))
    keep = res | ct | cf
    idx = [i for i in range(len(seqs)) if owner[i] in keep]
    sub = [seqs[i] for i in idx]
    sowner = [owner[i] for i in idx]
    for k, s in enumerate(sub): s.name = str(k).encode()
    print("targeted subset: %d proteins over %d plasmids\n" % (len(sub), len(keep)), flush=True)

    profs = sorted(f[:-4] for f in os.listdir(CJ)
                   if (f.startswith("FA_") or f.startswith("FATA_")) and f.endswith(".hmm"))
    nprof = collections.Counter()
    for p in profs:
        with pyhmmer.plan7.HMMFile(os.path.join(CJ, p + ".hmm")) as fh:
            mm = next(iter(fh))
        hits = set()
        for top in pyhmmer.hmmsearch([mm], sub, cpus=16, bit_cutoffs="gathering"):
            for h in top: hits.add(sowner[int(dec(h.name))])
        for acc in hits: nprof[acc] += 1
    print("=== how many Gram-positive profiles does each plasmid carry? ===")
    print("  %-22s %7s %9s %9s %9s %9s"
          % ("set", "n", ">=1", ">=2", ">=3", "NONE"))
    out = []
    for lab, s in (("residual (virb4 & !TFI)", res),
                   ("MPF_T control", ct), ("MPF_F control", cf)):
        n = len(s)
        g1 = sum(1 for a in s if nprof[a] >= 1)
        g2 = sum(1 for a in s if nprof[a] >= 2)
        g3 = sum(1 for a in s if nprof[a] >= 3)
        print("  %-22s %7d %8.1f%% %8.1f%% %8.1f%% %8.1f%%"
              % (lab, n, 100.0*g1/n, 100.0*g2/n, 100.0*g3/n, 100.0*(n-g1)/n))
        out.append({"set": lab, "n": n, "ge1_pct": round(100.0*g1/n, 1),
                    "ge2_pct": round(100.0*g2/n, 1), "ge3_pct": round(100.0*g3/n, 1),
                    "none_pct": round(100.0*(n-g1)/n, 1)})
    r, t, f = out[0], out[1], out[2]
    print("\n=== reading ===")
    print("  residual >=1 Gram+ profile : %.1f%%" % r["ge1_pct"])
    print("  MPF_T background           : %.1f%%" % t["ge1_pct"])
    print("  MPF_F background           : %.1f%%" % f["ge1_pct"])
    bg = max(t["ge1_pct"], f["ge1_pct"])
    if bg > 0:
        print("  enrichment over the worse background: %.1fx" % (r["ge1_pct"]/bg))
    print("\n  If the residual's rate is not well above BOTH backgrounds, the")
    print("  FA/FATA profiles are not Gram-positive-specific and the residual is")
    print("  NOT established as an MPF_FA container.")
    print("  %.1f%% of the residual carries NO Gram-positive profile at all." % r["none_pct"])
    dest = os.path.join(PROJ, "data", "anchors", "mpffa_residual_identity.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(out)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
