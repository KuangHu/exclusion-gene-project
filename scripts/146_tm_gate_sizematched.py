#!/usr/bin/env python3
"""TM gate, third attempt -- with the controls the TASK requires.

Attempt 1 (job 144) used synthetic polymers as negatives: passed 50/60 settings,
a bad gate not a good counter.

Attempt 2 (job 26074559) used real proteins but the WRONG POSITIVE: VirB6, a
polytopic inner-membrane transporter with median 351 aa and many TM segments. The
proteins we must actually detect are SMALL membrane-associated exclusion proteins
-- ExcA ~220 aa, TraS ~163 aa, TrbK 69 aa -- carrying one or two TM segments.
Separating polytopic transporters from large ATPases says nothing about whether a
163 aa TraS is distinguishable from a 163 aa soluble protein.

This gate:
  POSITIVES  the named exclusion families themselves, pulled by their own HMMs:
             PF10624 TraS, NF033891 ExcA, TIGR04359 TrbK, NF033894 Eex_IncN,
             plus the four seeds
  NEGATIVES  SIZE-MATCHED proteins drawn from the same plasmids -- for each
             positive, a random protein within +-15% length from the same
             accession. Size matching is essential: TM-segment count scales with
             length, so an unmatched negative set would separate on length alone
             and the counter would look good for the wrong reason.

Sweeps the COUNT THRESHOLD (>=1, >=2, >=3) as well as the window parameters; the
>=1 cut in attempt 2 was an arbitrary choice and the medians suggested it was the
wrong one.

Adopted only if some (win, thresh, min_len, k) gives positives >=80% with >=k TM
and size-matched negatives <=20%.
"""
import collections, os, random, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SEEDFA = "/global/scratch/users/kh36969/exclusion_gene/mmtest/seeds.faa"
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
EEX = ["PF10624", "NF033891", "TIGR04359", "NF033894"]


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
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    seqs, aas, owner = [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                p = line.rstrip("\n"); aas.append(p); owner.append(nm.split("|")[0])
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(aas)-1).encode(), sequence=p).digitize(alpha))
    print("cache %d proteins" % len(aas), flush=True)
    byacc = collections.defaultdict(list)
    for i, a in enumerate(owner): byacc[a].append(i)

    pos_idx = []
    for fam in EEX:
        p = os.path.join(HMM, fam + ".hmm")
        if not os.path.exists(p): continue
        with pyhmmer.plan7.HMMFile(p) as fh: m = next(iter(fh))
        got = set()
        for top in pyhmmer.hmmsearch([m], seqs, cpus=16, bit_cutoffs="gathering"):
            for h in top: got.add(int(dec(h.name)))
        pos_idx += sorted(got)
        print("  %-10s %d proteins" % (fam, len(got)), flush=True)
    pos_idx = sorted(set(pos_idx))
    pos = [aas[i] for i in pos_idx]
    nm = None
    for line in open(SEEDFA):
        if line[0] == ">": nm = line[1:].strip()
        else: pos.append(line.rstrip("\n"))
    print("  positives total %d, median %d aa"
          % (len(pos), sorted(len(x) for x in pos)[len(pos)//2]), flush=True)

    rng = random.Random(0)
    neg, miss = [], 0
    for i in pos_idx:
        L = len(aas[i]); acc = owner[i]
        cand = [j for j in byacc[acc] if j != i and 0.85*L <= len(aas[j]) <= 1.15*L]
        if cand: neg.append(aas[rng.choice(cand)])
        else: miss += 1
    print("  size-matched negatives %d (no match for %d), median %d aa"
          % (len(neg), miss, sorted(len(x) for x in neg)[len(neg)//2]), flush=True)

    print("\n=== gate: named exclusion proteins vs SIZE-MATCHED same-plasmid proteins ===")
    print("  %-5s %-7s %-7s %-3s %9s %9s  %s" % ("win","thresh","minlen","k","pos %>=k","neg %>=k","verdict"))
    ok = []
    for win in (15, 17, 19, 21):
        for thresh in (1.2, 1.6, 2.0, 2.4):
            for min_len in (15, 18, 21):
                pc = [tm_count(s, win, thresh, min_len) for s in pos]
                nc = [tm_count(s, win, thresh, min_len) for s in neg]
                for k in (1, 2, 3):
                    pp = 100.0*sum(1 for x in pc if x >= k)/len(pc)
                    np_ = 100.0*sum(1 for x in nc if x >= k)/len(nc)
                    good = pp >= 80.0 and np_ <= 20.0
                    if good: ok.append((win, thresh, min_len, k, pp, np_))
                    if good or (pp - np_) > 40:
                        print("  %-5d %-7.1f %-7d %-3d %8.1f%% %8.1f%%  %s"
                              % (win, thresh, min_len, k, pp, np_, "PASS" if good else ""))
    if not ok:
        print("\n  NO (window, threshold, min_len, k) separates named exclusion proteins")
        print("  from SIZE-MATCHED proteins on the same plasmids.")
        print("  TM COUNT IS UNAVAILABLE as a criterion. Reported as unavailable,")
        print("  not as a weak signal.")
        return 1
    b = max(ok, key=lambda x: x[4] - x[5])
    print("\n  ADOPTED win=%d thresh=%.1f min_len=%d k>=%d  (pos %.1f%%, neg %.1f%%)" % b)
    with open(os.path.join(PROJ, "scripts", "lib", "tm_params.txt"), "w") as fh:
        fh.write("%d\t%.1f\t%d\t%d\n" % (b[0], b[1], b[2], b[3]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
