#!/usr/bin/env python3
"""Remove unrecognised anchors from a positionally nominated set.

The same failure has now appeared three times:

    candidate pool   unrecognised VirB5   ->  186 unique removed
    recovered set    unrecognised VirB5   ->  190 unique removed
    VirB10 null      unrecognised VirB11  ->  under test

The cause is structural, not incidental. In a colinear virB operon the gene at
anchor+1 IS the next anchor -- VirB8|VirB9 are immediately adjacent on 6,954 of
7,002 plasmids. So a positional nomination set is built mostly out of anchors that
the HMM failed to detect, not out of random genes. Each anchor family carries its
own detection failure (VirB5 27.1%, VirB2 20.9%, VirB1 19.5%, VirB11 4.7%, ...),
and every one of those failures deposits its protein into the nomination set at a
perfectly fixed position.

Two consequences:

  * A nominated protein is guilty until cleared, and not only of being the ONE
    anchor you expected. Checking against VirB5 alone would have left VirB11
    contamination untouched.
  * The candidate set and the null set must be decontaminated IDENTICALLY, or
    their comparison is not a comparison. The in-slot calibration's 0.026 margin
    between null max (0.745) and the real families (0.771) was measured with the
    candidate side cleaned and the null side not.

phmmer against confirmed members is used rather than the family HMM because the
HMMs are precisely what missed these proteins: PF07996 scored 0/1221 candidates at
GA while phmmer matched 162 of the 228 in the VirB5 length band.
"""
import collections
import os

HMM_DIR = "/global/scratch/users/kh36969/exclusion_gene/hmm"

# every anchor family, with the threshold at which its confirmed members are taken
ANCHOR_FAMILIES = [
    ("VirB1", "anchor_SLT", None, 0),
    ("VirB2", "anchor_TrbC", None, 0),
    ("VirB3", "anchor_VirB3", None, 0),
    ("VirB4", "anchor_CagE_TrbE_VirB", None, 0),
    ("VirB4b", "anchor_VirB4_TrbE_N", None, 0),
    # VirB5: GA, NOT E<=1e-5. The E-value rule is WITHDRAWN (A13) -- it failed
    # both controls it was written for, was stricter than the GA it replaced, and
    # is database-size dependent. docs/virb5_evalue_threshold_invalid.md
    # GA alone gives 5,545 plasmids; the ADOPTED rule is the disjunction in
    # ANCHOR_DISJUNCTION below, which gives 6,881. A consumer that scans only
    # ANCHOR_FAMILIES gets the GA arm and UNDERCOUNTS VirB5 by 1,336 plasmids.
    ("VirB5", "anchor_T4SS", None, 150),
    ("VirB6", "TrbL", None, 0),
    ("VirB7", "anchor_P_T4SS_TraN", None, 0),
    ("VirB8", "anchor_VirB8", None, 0),
    ("VirB9", "anchor_CagX", None, 0),
    ("VirB10", "anchor_TrbI", None, 0),
    ("VirB11", "anchor_T2SSE", None, 0),
    ("VirD4", "anchor_T4SS-DNA_transf", None, 0),
    ("TraN", "anchor_F_T4SS_TraN", None, 0),
]
# Anchors whose detection is a DISJUNCTION: the Pfam arm above OR a CONJScan
# profile, with the SAME conjunct on both arms. Adopted only after passing a
# length-band and a cross-class control -- see data/anchors/anchor_disjunction_v1_1.tsv.
# VirB2, VirB3 and VirB6 were TESTED and REJECTED on purity; they are absent here
# deliberately and must not be added without re-running those controls.
ANCHOR_DISJUNCTION = {
    "VirB5": ("T4SS_T_virB5", 150),   # band 91.3%, cross-class 6.6%; +1,336 plasmids
    "VirB1": ("T4SS_T_virB1", 0),     # band 95.4%, cross-class 8.5%
}


def assert_no_evalue_anchor(families=None):
    """A13, applied to code: an anchor family may not carry an E-value cutoff."""
    bad = [(n, e) for n, f, e, m in (families or ANCHOR_FAMILIES) if e is not None]
    if bad:
        raise AssertionError(
            "A13: E-value cutoffs in ANCHOR_FAMILIES: %s\n"
            "  E = P * N makes these database-size dependent. Use GA or a bitscore.\n"
            "  See docs/virb5_evalue_threshold_invalid.md" % bad)
    return True


DEFAULT_E = 1e-3


def confirmed_members(seqs, aas, cpus=16):
    """family -> sorted unique protein sequences detected by that family's HMM."""
    import pyhmmer
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    out = {}
    for name, fam, ev, minaa in ANCHOR_FAMILIES:
        p = os.path.join(HMM_DIR, fam + ".hmm")
        if not os.path.exists(p):
            continue
        with pyhmmer.plan7.HMMFile(p) as fh:
            model = next(iter(fh))
        kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
        hits = set()
        for top in pyhmmer.hmmsearch([model], seqs, cpus=cpus, **kw):
            for h in top:
                i = int(dec(h.name))
                if minaa and len(aas[i]) <= minaa:
                    continue
                hits.add(aas[i])
        if hits:
            out[name] = sorted(hits)
    return out


def decontaminate(nominated, members, cpus=16, evalue=DEFAULT_E, verbose=True):
    """Drop nominated sequences homologous to ANY confirmed anchor member.

    `nominated` is a list of unique protein sequences. Returns (kept, dropped_by_family).
    """
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    dig = lambda L, t: [pyhmmer.easel.TextSequence(
        name=("%s%d" % (t, i)).encode(), sequence=s).digitize(alpha)
        for i, s in enumerate(L)]
    q = dig(nominated, "q")
    dropped = collections.defaultdict(set)
    for fam, ref in members.items():
        best = {}
        for top in pyhmmer.phmmer(q, dig(ref, "r"), cpus=cpus, E=1000.0):
            i = int(dec(top.query.name)[1:])
            for h in top:
                if i not in best or h.evalue < best[i]:
                    best[i] = h.evalue
        for i, e in best.items():
            if e <= evalue:
                dropped[fam].add(i)
    bad = set().union(*dropped.values()) if dropped else set()
    kept = [s for i, s in enumerate(nominated) if i not in bad]
    if verbose:
        print("  decontamination: %d nominated -> %d kept (%d removed, %.1f%%)"
              % (len(nominated), len(kept), len(bad),
                 100.0 * len(bad) / max(1, len(nominated))))
        for fam, v in sorted(dropped.items(), key=lambda kv: -len(kv[1])):
            if v:
                print("      %-8s %d" % (fam, len(v)))
    return kept, {k: sorted(v) for k, v in dropped.items()}
