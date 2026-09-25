#!/usr/bin/env python3
"""Re-gate the TM counter against REAL proteins, not synthetic polymers.

The first gate (job: 144) passed at 50 of 60 parameter settings. That is not a
good counter, it is a bad gate: the negatives were poly-E/K, charged-random and
proline-rich synthetic sequences, which ANY hydrophobicity method separates.

The failure mode that matters is a GLOBULAR protein with a hydrophobic core being
called transmembrane. So the negatives here are real, abundant, known-soluble
proteins pulled from the cache by their own HMMs:

    VirB11 / T2SSE      cytoplasmic ATPase
    VirB4 / CagE_TrbE   cytoplasmic ATPase
    VirD4 / T4SS-DNA_transf  coupling ATPase (largely cytoplasmic)

and the positives are real polytopic membrane proteins:

    VirB6 / TrbL        polytopic inner-membrane
    plus the four seed exclusion proteins (TrbK, Eex) -- known membrane-associated

A setting is adopted ONLY if it separates these. If the bands overlap, TM count is
reported UNAVAILABLE rather than as a weak criterion.
"""
import collections, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SEEDFA = "/global/scratch/users/kh36969/exclusion_gene/mmtest/seeds.faa"
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
SOLUBLE = [("VirB11", "anchor_T2SSE"), ("VirB4", "anchor_CagE_TrbE_VirB"),
           ("VirD4", "anchor_T4SS-DNA_transf")]
MEMBRANE = [("VirB6", "TrbL")]
NSAMP = 300


def tm_count(seq, win, thresh, min_len):
    n = len(seq)
    if n < win: return 0
    vals = [sum(KD.get(c, 0.0) for c in seq[i:i+win]) / win for i in range(n - win + 1)]
    need = max(1, min_len - win + 1)
    segs = run = 0
    for v in vals:
        if v >= thresh: run += 1
        else:
            if run >= need: segs += 1
            run = 0
    if run >= need: segs += 1
    return segs


def main():
    require_compute_node()
    import pyhmmer, random
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    seqs, aas = [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        for line in open(os.path.join(CACHE, fn)):
            if line[0] != ">":
                p = line.rstrip("\n"); aas.append(p)
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(aas)-1).encode(), sequence=p).digitize(alpha))
    print("cache %d proteins" % len(aas), flush=True)

    def grab(fam):
        with pyhmmer.plan7.HMMFile(os.path.join(HMM, fam + ".hmm")) as fh:
            m = next(iter(fh))
        out = []
        for top in pyhmmer.hmmsearch([m], seqs, cpus=16, bit_cutoffs="gathering"):
            for h in top: out.append(aas[int(dec(h.name))])
        return out

    rng = random.Random(0)
    groups = {}
    for lab, fam in SOLUBLE + MEMBRANE:
        g = grab(fam)
        groups[lab] = rng.sample(g, min(NSAMP, len(g)))
        print("  %-8s %6d proteins, sampled %d, median %d aa"
              % (lab, len(g), len(groups[lab]),
                 sorted(len(x) for x in groups[lab])[len(groups[lab])//2]), flush=True)
    seed = []
    nm = None
    for line in open(SEEDFA):
        if line[0] == ">": nm = line[1:].strip()
        else: seed.append(line.rstrip("\n"))
    groups["seed_eex"] = seed
    print("  %-8s %6d (the four named exclusion controls)" % ("seed_eex", len(seed)))

    sol = [s for lab, _ in SOLUBLE for s in groups[lab]]
    mem = groups["VirB6"]
    print("\n=== gate: real soluble ATPases vs real polytopic membrane ===")
    print("  %-5s %-7s %-7s %9s %9s %9s %9s  %s"
          % ("win", "thresh", "minlen", "sol med", "sol %>=1", "mem med", "mem %>=1", "verdict"))
    ok_set = []
    for win in (15, 17, 19, 21):
        for thresh in (1.2, 1.6, 2.0, 2.4):
            for min_len in (15, 18, 21):
                sc = [tm_count(s, win, thresh, min_len) for s in sol]
                mc = [tm_count(s, win, thresh, min_len) for s in mem]
                sp = 100.0 * sum(1 for x in sc if x >= 1) / len(sc)
                mp = 100.0 * sum(1 for x in mc if x >= 1) / len(mc)
                good = sp <= 20.0 and mp >= 80.0
                if good: ok_set.append((win, thresh, min_len, sp, mp))
                print("  %-5d %-7.1f %-7d %9d %8.1f%% %9d %8.1f%%  %s"
                      % (win, thresh, min_len, sorted(sc)[len(sc)//2], sp,
                         sorted(mc)[len(mc)//2], mp, "PASS" if good else ""))
    if not ok_set:
        print("\n  NO SETTING separates real soluble ATPases from real membrane proteins")
        print("  at (soluble <=20%% with >=1 TM) and (membrane >=80%% with >=1 TM).")
        print("  TM COUNT IS UNAVAILABLE. It is not reported as a weak criterion:")
        print("  a predictor failing its own controls yields numbers that look like")
        print("  data and are not. EMBOSS tmap was rejected for the same reason.")
        return 1
    win, thresh, min_len, sp, mp = max(ok_set, key=lambda x: x[4] - x[3])
    print("\n  ADOPTED win=%d thresh=%.1f min_len=%d (soluble %.1f%%, membrane %.1f%%)"
          % (win, thresh, min_len, sp, mp))
    for lab in ("seed_eex",):
        c = [tm_count(s, win, thresh, min_len) for s in groups[lab]]
        print("  %-10s TM counts %s" % (lab, c))
    with open(os.path.join(PROJ, "scripts", "lib", "tm_params.txt"), "w") as fh:
        fh.write("%d\t%.1f\t%d\n" % (win, thresh, min_len))
    return 0


if __name__ == "__main__":
    sys.exit(main())
