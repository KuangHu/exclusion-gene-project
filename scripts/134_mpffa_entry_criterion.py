#!/usr/bin/env python3
"""MPF_FA §2: entry criterion, measured against the three existing catalogues.

`T4SS_virb4` fires on all five Gram-positive seeds (91.0-279.8, job 25974613) and
is the only profile that does. But it is MANDATORY IN EVERY CONJScan plasmid class
definition -- FA, FATA, T, F, I alike -- so by construction it will also admit
MPF_T and MPF_F plasmids. Seed recall is already known (5/5) and is NOT the
question here.

THE QUESTION IS CROSS-CLASS CONTAMINATION.

"Already admitted" is read from the CATALOGUES, not from a rescan -- those files
ARE the admitted sets:
    MPF_T  data/release/v1.1/catalogue_MPF_T_v1.1.tsv   7,436
    MPF_F  data/release/v1/catalogue_MPF_F_v1.tsv      11,120
    MPF_I  data/release/v1.1/catalogue_MPF_I_v1.tsv     3,732

Decision rule stated BEFORE the measurement:
    contamination <= 30%   -> `T4SS_virb4` alone is usable as entry
    contamination >  30%   -> it is not; test the CONJUNCTION
                              `T4SS_virb4 AND NOT (T or F or I)`
                              and the FA/FATA accessory profiles as alternates

The known failure mode this guards against: MPF_F rejected four candidate entry
families precisely because they swallowed 6,000-9,000 plasmids of other classes.

The opposite trap is also known: generic `T4SS_virb4` scores ZERO on R64 while the
class-specific `T4SS_I_traU` scores 1,460. A generic profile can fail on a class
entirely. Here it does fire on all five Gram+ seeds, so that failure mode is
excluded FOR THIS CLASS -- but only for this class.

MOB-suite purity is NOT used as a gate: MPF_I showed it to be class-dependent.
"""
import collections, csv, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
REL = os.path.join(PROJ, "data", "release")
CATS = [("MPF_T", os.path.join(REL, "v1.1", "catalogue_MPF_T_v1.1.tsv")),
        ("MPF_F", os.path.join(REL, "v1", "catalogue_MPF_F_v1.tsv")),
        ("MPF_I", os.path.join(REL, "v1.1", "catalogue_MPF_I_v1.tsv"))]
LIMIT = 30.0


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    admitted, allprev = {}, set()
    for name, path in CATS:
        s = {l.split("\t")[0] for i, l in enumerate(open(path)) if i}
        admitted[name] = s; allprev |= s
        print("  %-6s already admitted: %6d" % (name, len(s)))
    print("  union of the three          : %6d\n" % len(allprev), flush=True)

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
    NDB = len(set(owner))
    print("cache: %d proteins over %d plasmids\n" % (len(seqs), NDB), flush=True)

    _cache = {}

    def scan(p):
        """Cached: every profile is searched over the 7.7M proteins at most once."""
        if p in _cache:
            return _cache[p]
        with pyhmmer.plan7.HMMFile(p) as fh: m = next(iter(fh))
        out = set()
        for top in pyhmmer.hmmsearch([m], seqs, cpus=16, bit_cutoffs="gathering"):
            for h in top: out.add(owner[int(dec(h.name))])
        _cache[p] = out
        return out

    # profiles that fired on >=1 seed in job 25974613, plus the mandatory set
    sr = collections.defaultdict(set)
    for r in csv.DictReader(open(os.path.join(PROJ, "data", "anchors",
                                              "mpffa_seed_readout.tsv")), delimiter="\t"):
        sr[r["profile"]].add(r["seed"])
    fa = sorted(f[:-4] for f in os.listdir(CJ) if f.startswith("FA_"))
    fata = sorted(f[:-4] for f in os.listdir(CJ) if f.startswith("FATA_"))
    cands = ["T4SS_virb4", "T4SS_tcpA", "T4SS_t4cp1", "T4SS_t4cp2"] + fa + fata

    print("=== candidate entry criteria: admission and CROSS-CLASS contamination ===")
    print("  contamination = share of admitted plasmids ALREADY in T, F or I")
    print("  %-18s %8s %7s %9s %9s %9s %10s %6s"
          % ("profile", "admitted", "%db", "in MPF_T", "in MPF_F", "in MPF_I",
             "contam%", "seeds"))
    rows = []
    for p in cands:
        path = os.path.join(CJ, p + ".hmm")
        if not os.path.exists(path): continue
        s = scan(path)
        if not s: continue
        t, f, i = (len(s & admitted[k]) for k in ("MPF_T", "MPF_F", "MPF_I"))
        c = 100.0 * len(s & allprev) / len(s)
        ns = len(sr.get(p, ()))
        rows.append({"profile": p, "admitted": len(s), "pct_db": round(100.0*len(s)/NDB, 2),
                     "in_MPF_T": t, "in_MPF_F": f, "in_MPF_I": i,
                     "contamination_pct": round(c, 1), "n_seeds": ns,
                     "residual_novel": len(s - allprev)})
        if p in ("T4SS_virb4", "T4SS_tcpA") or ns or len(s) > 200:
            print("  %-18s %8d %6.2f%% %9d %9d %9d %9.1f%% %6d"
                  % (p, len(s), 100.0*len(s)/NDB, t, f, i, c, ns), flush=True)
        if p == "T4SS_virb4":
            vb4 = s

    print("\n=== the decision, against the rule stated before the run ===")
    c = 100.0 * len(vb4 & allprev) / len(vb4)
    print("  T4SS_virb4 admits %d plasmids (%.2f%% of the database)" % (len(vb4), 100.0*len(vb4)/NDB))
    print("  of those, %d are already in T/F/I -> contamination %.1f%%" % (len(vb4 & allprev), c))
    print("  threshold was %.0f%% -> %s" % (LIMIT, "USABLE ALONE" if c <= LIMIT else
                                            "NOT usable alone; test the conjunction"))
    res = vb4 - allprev
    print("\n=== conjunction: T4SS_virb4 AND NOT (T or F or I) ===")
    print("  residual: %d plasmids (%.2f%% of the database)" % (len(res), 100.0*len(res)/NDB))
    print("  this is the MPF_FA candidate container if the conjunction is adopted")
    print("\n  overlap of the residual with the FA/FATA accessory profiles:")
    for p in fa + fata:
        path = os.path.join(CJ, p + ".hmm")
        if not os.path.exists(path): continue
        s = scan(path)
        n = len(s & res)
        if n:
            print("    %-18s %6d of the residual (%.1f%%)  seeds %d"
                  % (p, n, 100.0*n/len(res), len(sr.get(p, ()))))
    dest = os.path.join(PROJ, "data", "anchors", "mpffa_entry_criterion.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(sorted(rows, key=lambda r: -r["admitted"]))
    print("\nwrote %s (%d rows)" % (dest, len(rows)))
    print("\n  NOTE on seeds: 5/5 was measured on the SEED FILES (job 25974613),")
    print("  not in PLSDB. pCF10 is absent from PLSDB and ICEBs1 is chromosomal,")
    print("  so a PLSDB-based seed check would read 0/5 for reasons that have")
    print("  nothing to do with the criterion (the A12 failure mode).")


if __name__ == "__main__":
    main()
