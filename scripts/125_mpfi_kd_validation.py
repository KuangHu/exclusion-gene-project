#!/usr/bin/env python3
"""Is hydrophobicity the MPF_I admission criterion? Validate it as lipobox was.

MPF_I's own ExcA carries NO lipobox on any seed (R64 220 aa, R621a 204, pEK204
220) but is strongly hydrophobic -- max Kyte-Doolittle 1.90, 2.93, 1.90, at or
above the MPF_T exclusion controls at 1.04-1.97. So the candidates' 0% lipobox is
consistent with the class rather than disqualifying, and hydrophobicity is the
criterion to test.

That explanation was checked BEFORE the candidates were judged, not after, which
is the only ordering that is not post-hoc rationalisation.

BACKGROUND CONTROL. MPF_T's lipobox validation compared slot ORFs against ORFs in
size-matched intergenic gaps on the same plasmids. That control does not transfer:
the MPF_I candidates are CALLED GENES, not six-frame ORFs, so the matched control
is other called proteins of the same length band (200-240 aa) on the same
plasmids. Comparing against random ORFs or against all proteins would only show
that membrane proteins are hydrophobic -- not that the slot occupant is special.

Three groups, one measure (max KD over a 19-residue window):
    ExcA          the named exclusion protein, positive reference
    slot occupant TraY+1 on ExcA-negative plasmids, the candidates
    background    200-240 aa proteins elsewhere on the SAME plasmids
"""
import argparse
import collections
import csv
import os
import statistics
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
      'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
      'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
BAND = (200, 240)


def maxkd(s, w=19):
    if len(s) < w:
        return sum(KD.get(c, 0) for c in s) / max(1, len(s))
    return max(sum(KD.get(c, 0) for c in s[i:i + w]) / w for i in range(len(s) - w + 1))


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    csv.field_size_limit(10 ** 8)

    surv = {r["accession"]: r for r in csv.DictReader(
        open(os.path.join(PROJ, "data", "anchors", "mpfi_architecture_survey.tsv")),
        delimiter="\t")}
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
                owner.append(f[0]); idx.append(int(f[1])); strand.append(int(f[4]))
                aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    bypos = {(owner[i], idx[i]): i for i in range(len(owner))}
    byacc = collections.defaultdict(list)
    for i in range(len(owner)):
        byacc[owner[i]].append(i)
    print("cache: %d proteins" % len(seqs), flush=True)

    def scan(path):
        with pyhmmer.plan7.HMMFile(path) as fh:
            m = next(iter(fh))
        out = collections.defaultdict(list)
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                out[owner[i]].append((idx[i], strand[i], h.score, i))
        return out
    tray = scan(os.path.join(CJ, "T4SS_I_traY.hmm"))
    exca = scan(os.path.join(HMM, "NF033891.hmm"))
    print("TraY %d, ExcA %d" % (len(tray), len(exca)), flush=True)
    best = lambda d, acc: max(d[acc], key=lambda x: x[2]) if acc in d else None

    exc_seqs, cand_seqs, bg_seqs = [], [], []
    cand_accs = []
    for acc in surv:
        ty = best(tray, acc)
        if not ty:
            continue
        ex = best(exca, acc)
        if ex:
            exc_seqs.append(aas[ex[3]])
            continue
        t = 1 if ty[1] == 1 else -1
        j = bypos.get((acc, ty[0] + t))
        if j is None or strand[j] != ty[1]:
            continue
        cand_seqs.append(aas[j]); cand_accs.append(acc)
        # background: same-length called proteins elsewhere on THIS plasmid
        for k in byacc[acc]:
            if k == j or k == ty[3]:
                continue
            if BAND[0] <= len(aas[k]) <= BAND[1]:
                bg_seqs.append(aas[k])

    print("\n=== max Kyte-Doolittle (19-residue window) ===")
    print("%-34s %7s %8s %8s %8s %9s" % ("group", "n", "median", "Q1", "Q3", ">=1.90"))
    ref = None
    for lab, g in (("ExcA (named, positive reference)", exc_seqs),
                   ("TraY+1 occupant, ExcA-negative", cand_seqs),
                   ("background 200-240 aa, same plasmids", bg_seqs)):
        if not g:
            continue
        v = sorted(maxkd(s) for s in g)
        hi = sum(1 for x in v if x >= 1.90)
        print("%-34s %7d %8.2f %8.2f %8.2f %8.1f%%"
              % (lab, len(v), statistics.median(v), v[len(v)//4], v[3*len(v)//4],
                 100.0 * hi / len(v)))
        if lab.startswith("ExcA"):
            ref = statistics.median(v)
    print("\n  seed ExcA reference: R64 1.90, R621a 2.93, pEK204 1.90")

    uq = collections.Counter(cand_seqs)
    print("\n=== the unique candidate sequences ===")
    print("  %-6s %5s %8s  %s" % ("n", "aa", "maxKD", "N-terminus"))
    for s, n in uq.most_common(12):
        print("  %-6d %5d %8.2f  %s" % (n, len(s), maxkd(s), s[:36]))
    print("\n  candidates: %d records, %d unique" % (len(cand_seqs), len(uq)))

    ca = sorted(maxkd(s) for s in cand_seqs)
    bg = sorted(maxkd(s) for s in bg_seqs)
    if ca and bg:
        import bisect
        pct = 100.0 * bisect.bisect_left(bg, statistics.median(ca)) / len(bg)
        print("  the candidates' median KD sits at the %.1fth percentile of the "
              "same-length background" % pct)
        print("  (a criterion that separates should put it far above 50)")
    dest = os.path.join(PROJ, "data", "anchors", "mpfi_kd_validation.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["group", "aa_len", "max_kd"])
        for lab, g in (("excA", exc_seqs), ("candidate", cand_seqs), ("background", bg_seqs)):
            for s in g:
                w.writerow([lab, len(s), round(maxkd(s), 3)])
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
