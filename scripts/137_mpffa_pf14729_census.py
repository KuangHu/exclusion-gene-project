#!/usr/bin/env python3
"""MPF_FA: class-scale PF14729 census. Does the slot exist beyond ICEBs1?

The seed check FAILED to generalise (job 26035629): conG+3 is confirmed on ICEBs1
alone. On pAD1 the +2 occupant is a type II toxin-antitoxin system (antitoxin MazE
at +1); on pLS20 and pAM373 the occupant is ANTI-ORIENTED to the anchor and so is
almost certainly not co-transcribed; on pCF10 the anchor is absent entirely.

Two clean outcomes, decided in advance:

  DUF4467 recurs at a CONSISTENT offset  -> the slot exists; the seed failure is a
                                            sampling problem (all five seeds in one
                                            stratum), as with MPF_I's order#1
  scattered or absent                    -> the slot does not exist, and the
                                            catalogue must carry NO SLOT COLUMN
                                            rather than a warned one

THE TEST IS OFFSET CONSISTENCY, NOT HIT COUNT. 200 hits spread over -5..+8 is
indistinguishable from random adjacency.

Decision rule borrowed from MPF_T, stated before the run:
    a dominant offset holding >=30% of hits, coherent  -> architecture (IncI2 +2)
    top offset ~5-11% with no consistent alternative   -> noise (MPF_I +2)

Offsets are in the TRANSCRIPTION FRAME of the anchor: (idx_hit - idx_anchor) * s,
s = +1 if the anchor is on the plus strand else -1. Same-strand fraction is
reported separately -- an anti-oriented neighbour is not a slot occupant.

Also writes the admitted accession list, which §2 did not persist.
"""
import collections, csv, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
REL = os.path.join(PROJ, "data", "release")
CATS = [("MPF_T", os.path.join(REL, "v1.1", "catalogue_MPF_T_v1.1.tsv")),
        ("MPF_F", os.path.join(REL, "v1", "catalogue_MPF_F_v1.tsv")),
        ("MPF_I", os.path.join(REL, "v1.1", "catalogue_MPF_I_v1.tsv"))]
MIN_PROF = 2


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    allprev = set()
    for _, p in CATS:
        allprev |= {l.split("\t")[0] for i, l in enumerate(open(p)) if i}

    seqs, owner, idx, st, en, strand, aas = [], [], [], [], [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                owner.append(f[0]); idx.append(int(f[1])); st.append(int(f[2]))
                en.append(int(f[3])); strand.append(int(f[4])); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(), sequence=line.rstrip("\n")).digitize(alpha))
    print("cache: %d proteins over %d plasmids" % (len(seqs), len(set(owner))), flush=True)

    def search(path, pool, pool_owner):
        with pyhmmer.plan7.HMMFile(path) as fh: m = next(iter(fh))
        hits = []
        for top in pyhmmer.hmmsearch([m], pool, cpus=16, bit_cutoffs="gathering"):
            for h in top: hits.append((int(dec(h.name)), round(h.score, 1)))
        return hits

    vb4 = {owner[i] for i, _ in search(os.path.join(CJ, "T4SS_virb4.hmm"), seqs, owner)}
    res = vb4 - allprev
    print("virb4 %d ; residual %d" % (len(vb4), len(res)), flush=True)

    sel = [i for i in range(len(seqs)) if owner[i] in res]
    sub = [seqs[i] for i in sel]
    for k, s in enumerate(sub): s.name = str(k).encode()
    sown = [owner[i] for i in sel]
    profs = sorted(f[:-4] for f in os.listdir(CJ)
                   if (f.startswith("FA_") or f.startswith("FATA_")) and f.endswith(".hmm"))
    nprof = collections.Counter()
    anch = {}     # acc -> {profile: (global_i, score)}
    for p in profs:
        got = set()
        for k, sc in search(os.path.join(CJ, p + ".hmm"), sub, sown):
            g = sel[k]; a = owner[g]
            got.add(a)
            if p in ("FA_orf14", "FA_orf15"):
                d = anch.setdefault(a, {})
                if p not in d or sc > d[p][1]: d[p] = (g, sc)
        for a in got: nprof[a] += 1
    adm = sorted(a for a in res if nprof[a] >= MIN_PROF)
    print("admitted (>=%d FA/FATA profiles): %d\n" % (MIN_PROF, len(adm)), flush=True)
    with open(os.path.join(PROJ, "data", "anchors", "mpffa_admitted.tsv"), "w") as fh:
        fh.write("accession\tn_fa_fata_profiles\n")
        for a in adm: fh.write("%s\t%d\n" % (a, nprof[a]))

    admset = set(adm)
    sel2 = [i for i in range(len(seqs)) if owner[i] in admset]
    sub2 = [seqs[i] for i in sel2]
    for k, s in enumerate(sub2): s.name = str(k).encode()
    pf = search(os.path.join(HMM, "PF14729.hmm"), sub2, None)
    hits = collections.defaultdict(list)
    for k, sc in pf:
        g = sel2[k]; hits[owner[g]].append((g, sc))
    print("=== PF14729 (DUF4467) census over the %d admitted plasmids ===" % len(adm))
    print("  plasmids with >=1 DUF4467 hit: %d (%.1f%%)"
          % (len(hits), 100.0 * len(hits) / max(1, len(adm))))
    tot = sum(len(v) for v in hits.values())
    print("  total DUF4467 proteins: %d" % tot)
    print("\n  class-scale comparison (named family recoverable across the class):")
    print("    MPF_F TraS  2,725 | MPF_I ExcA 2,213 | MPF_T TrbK 224 | MPF_FA %d" % len(hits))
    if len(hits) < 100:
        print("    DUF4467 is effectively ABSENT at class scale. Offset shares below")
        print("    are reported for the record but CANNOT establish a slot.")
    if not hits:
        print("\n  NO DUF4467 ANYWHERE IN THE CLASS -> the slot does not exist beyond")
        print("  ICEBs1. The catalogue must carry NO slot column.")
        return

    rows = []
    for aname in ("FA_orf14", "FA_orf15"):
        off = collections.Counter(); same = 0; n = 0
        for acc, hl in hits.items():
            a = anch.get(acc, {}).get(aname)
            if not a: continue
            ai, asc = a
            s = 1 if strand[ai] == 1 else -1
            for g, sc in hl:
                d = (idx[g] - idx[ai]) * s
                off[d] += 1; n += 1
                if strand[g] == strand[ai]: same += 1
                rows.append({"accession": acc, "anchor": aname, "anchor_bit": asc,
                             "offset": d, "same_strand": int(strand[g] == strand[ai]),
                             "duf_bit": sc, "duf_aa": len(aas[g])})
        print("\n  --- offsets relative to %s (n=%d hit-anchor pairs, %d plasmids) ---"
              % (aname, n, sum(1 for a in hits if aname in anch.get(a, {}))))
        if not n:
            print("    anchor absent on every DUF4467-positive plasmid"); continue
        print("    same strand as anchor: %d/%d = %.1f%%" % (same, n, 100.0*same/n))
        print("    %8s %7s %8s" % ("offset", "n", "share"))
        for d, c in sorted(off.items(), key=lambda x: -x[1])[:10]:
            bar = "#" * int(40.0 * c / n)
            print("    %+8d %7d %7.1f%%  %s" % (d, c, 100.0*c/n, bar))
        top, tc = max(off.items(), key=lambda x: x[1])
        sh = 100.0 * tc / n
        # MINIMUM-n GUARD. The >=30% rule was calibrated on MPF_T/MPF_I where n
        # was in the thousands. A share computed on a handful of observations has
        # no power to discriminate: 4 of 11 is 36% and is also chance. Without
        # this guard the script reported "slot supported" from 11 observations.
        MIN_N = 100
        if n < MIN_N:
            print("\n    dominant offset %+d at %.1f%% (n=%d) -> NOT READ." % (top, sh, n))
            print("    n < %d: the >=30%% rule was calibrated on thousands of" % MIN_N)
            print("    observations and cannot discriminate at this size. A share")
            print("    is not evidence when the denominator is a handful.")
        else:
            print("\n    dominant offset %+d at %.1f%% (n=%d) -> %s" % (top, sh, n,
                  "CONSISTENT (>=30%): slot supported" if sh >= 30 else
                  "SCATTERED (<30%): indistinguishable from random adjacency"))
    if rows:
        dest = os.path.join(PROJ, "data", "anchors", "mpffa_pf14729_census.tsv")
        with open(dest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(rows)
        print("\nwrote %s (%d rows)" % (dest, len(rows)))


if __name__ == "__main__":
    main()
