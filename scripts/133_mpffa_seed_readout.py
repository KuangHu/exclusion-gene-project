#!/usr/bin/env python3
"""MPF_FA §1.1: five Gram-positive seeds against every FA and FATA profile.

CONJScan has NO `T4SS_FA_*` profile set. Gram-positive splits into two classes
with different prefixes (docs/mpf_fa_precheck.md):

    FA      7 profiles   FA_orf13 orf14 orf15 orf17a orf17b orf19 orf23
    FATA   26 profiles   prg* (pCF10), trs* (pGO1), cd* (C. difficile),
                         gbs* (S. agalactiae)

Every class-specific profile in BOTH sets is `accessory`; only virb4, a coupling
protein and a relaxase are mandatory. So this readout answers two questions that
must not be merged:

  Q1  which profile set does each seed belong to -- FA, FATA, both, neither?
  Q2  which individual profiles fire, and on what?

**ICEBs1's class is a MEASUREMENT, not an inference from B. subtilis.** pCF10 is
the FATA reference strain, so pCF10 hitting FATA is the POSITIVE CONTROL: if the
FATA profiles do not fire on pCF10, the scan is broken and nothing else is read.

Mandatory-gene profiles are scanned too (virb4, t4cp1/2, tcpA, MOB*) -- they are
what an entry criterion can actually be built on, since the class-specific ones
are optional by CONJScan's own definition.

Exclusion families are scanned LAST and POST-HOC. PF14729 on ICEBs1 yddJ is the
named positive control for the slot; it never enters the entry criterion.
"""
import collections, csv, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
GB = os.path.join(PROJ, "data", "seed", "genbank")
SEEDS = [("ICEBs1", "ICEBs1__CP171645.1"), ("pLS20", "pLS20__AB615352.1"),
         ("pCF10", "pCF10__AY855841.2"), ("pAD1", "pAD1__CP046109.1"),
         ("pAM373", "pAM373__AE002565.1")]
MAND = ["T4SS_virb4", "T4SS_tcpA", "T4SS_t4cp1", "T4SS_t4cp2"] + \
       ["T4SS_MOB" + x for x in ("B", "C", "F", "H", "P1", "P2", "P3", "Q", "T", "V")]
EEX = ["PF14729", "PF10624", "NF033891", "NF033894", "NF041429", "TIGR04359"]
CONTROL = ("pCF10", "FATA")   # pCF10 IS the FATA reference -- must fire


def main():
    require_compute_node()
    import pyhmmer
    from Bio import SeqIO
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner, tag, aas = [], [], [], []
    for name, fn in SEEDS:
        rec = next(SeqIO.parse(os.path.join(GB, fn + ".gb"), "genbank"))
        k = 0
        for f in rec.features:
            if f.type != "CDS":
                continue
            p = f.qualifiers.get("translation", [None])[0]
            if not p:
                continue
            g = (f.qualifiers.get("gene") or f.qualifiers.get("locus_tag")
                 or f.qualifiers.get("product") or ["?"])[0]
            owner.append(name); tag.append((k, g, int(f.location.start) + 1,
                                            int(f.location.end),
                                            1 if f.location.strand == 1 else -1))
            aas.append(p); k += 1
            seqs.append(pyhmmer.easel.TextSequence(
                name=str(len(seqs)).encode(), sequence=p).digitize(alpha))
        print("  %-8s %s  %d CDS with translation" % (name, fn, k), flush=True)
    print("total %d proteins\n" % len(seqs), flush=True)

    fa = sorted(f[:-4] for f in os.listdir(CJ) if f.startswith("FA_"))
    fata = sorted(f[:-4] for f in os.listdir(CJ) if f.startswith("FATA_"))
    print("profiles: FA %d, FATA %d, mandatory %d\n" % (len(fa), len(fata), len(MAND)))

    def scan(path):
        with pyhmmer.plan7.HMMFile(path) as fh:
            m = next(iter(fh))
        out = collections.defaultdict(list)
        for top in pyhmmer.hmmsearch([m], seqs, cpus=4, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                out[owner[i]].append((round(h.score, 1), tag[i], i))
        return out

    # A named control profile that is MISSING must not be skipped silently.
    # PF14729/conJ is the stated positive control for the MPF_FA slot; the first
    # run of this script reported nothing for it because PF14729.hmm did not
    # exist on disk and the loop skipped it. Absence of a model read as absence
    # of a hit is the same class of error as reading "we did not look" as a
    # negative result.
    miss = [p for p in EEX if not os.path.exists(os.path.join(HMM, p + ".hmm"))]
    if miss:
        raise SystemExit(
            "MISSING exclusion-family models: %s\n"
            "  These are NAMED CONTROLS. A run without them cannot report on the\n"
            "  slot's positive control. Fetch from Pfam-A (hmmfetch) before rerunning."
            % ", ".join(miss))

    res, rows = {}, []
    for grp, names, base in (("FA", fa, CJ), ("FATA", fata, CJ),
                             ("mandatory", MAND, CJ), ("eex", EEX, HMM)):
        for p in names:
            path = os.path.join(base, p + ".hmm")
            if not os.path.exists(path):
                if grp != "eex":
                    print("  note: %s/%s.hmm absent" % (grp, p))
                continue
            res[(grp, p)] = scan(path)

    # --- positive control FIRST ------------------------------------------
    seed, grp = CONTROL
    fired = [p for (g, p), r in res.items() if g == grp and seed in r]
    print("=== POSITIVE CONTROL: %s against the %s set ===" % (seed, grp))
    print("  %d of %d %s profiles fire on %s" % (len(fired), len(fata), grp, seed))
    if not fired:
        raise SystemExit("CONTROL FAILED: %s is the %s reference strain and no %s "
                         "profile fires. The scan is broken; no other result is read."
                         % (seed, grp, grp))
    print("  CONTROL PASSES -- proceeding\n")

    # --- Q1: class assignment -------------------------------------------
    print("=== Q1: which profile set does each seed belong to? ===")
    print("  %-8s %10s %10s   %s" % ("seed", "FA hits", "FATA hits", "assignment"))
    assign = {}
    for name, _ in SEEDS:
        nfa = sum(1 for (g, p), r in res.items() if g == "FA" and name in r)
        nft = sum(1 for (g, p), r in res.items() if g == "FATA" and name in r)
        # One hit on each side is NOT dual membership. Require a margin and a
        # minimum, or the seed is UNASSIGNED. Calling 1-vs-1 "both" overstates it.
        if max(nfa, nft) < 2:
            a = "UNASSIGNED"
        elif nfa >= 2 * max(1, nft):
            a = "FA"
        elif nft >= 2 * max(1, nfa):
            a = "FATA"
        else:
            a = "AMBIGUOUS"
        assign[name] = (nfa, nft, a)
        print("  %-8s %10d %10d   %s" % (name, nfa, nft, a))
    print("\n  UNASSIGNED = fewer than 2 hits on the stronger side. Every")
    print("  class-specific profile is accessory, so a genuine MPF_FA element may")
    print("  hit almost none -- that is not evidence for either set.")
    print("  AMBIGUOUS = both sides fire without a 2x margin.\n")

    # --- Q2: family x seed matrix ---------------------------------------
    print("=== Q2: family x seed matrix (bitscore at GA; '.' = no hit) ===")
    hdr = "  %-18s %-6s" % ("profile", "set") + "".join("%9s" % n for n, _ in SEEDS)
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for (g, p) in sorted(res, key=lambda k: (["mandatory", "FA", "FATA", "eex"].index(k[0]), k[1])):
        r = res[(g, p)]
        if not r:
            continue
        line = "  %-18s %-6s" % (p, g)
        for n, _ in SEEDS:
            line += "%9s" % (max(x[0] for x in r[n]) if n in r else ".")
        print(line)
        for n, _ in SEEDS:
            for sc, t, i in sorted(r.get(n, []), key=lambda x: -x[0])[:1]:
                rows.append({"profile": p, "set": g, "seed": n, "bitscore": sc,
                             "gene": t[1], "cds_index": t[0], "start": t[2],
                             "end": t[3], "strand": "+" if t[4] == 1 else "-",
                             "aa_len": len(aas[i])})
    dest = os.path.join(PROJ, "data", "anchors", "mpffa_seed_readout.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))

    print("\n=== ICEBs1 class assignment (the question this step exists to answer) ===")
    nfa, nft, a = assign["ICEBs1"]
    print("  FA %d, FATA %d -> %s" % (nfa, nft, a))
    print("  Measured, not inferred from B. subtilis.")
    print("\n=== named positive control for the slot ===")
    any_eex = False
    for p in EEX:
        r = res.get(("eex", p), {})
        for n in r:
            for sc, t, i in sorted(r[n], key=lambda x: -x[0])[:1]:
                any_eex = True
                print("  %-10s %-8s %-14s %6.1f  %d aa  %d..%d"
                      % (p, n, t[1], sc, len(aas[i]), t[2], t[3]))
    if not any_eex:
        print("  NONE of the six exclusion families fires on ANY of the five seeds.")
        print("  All six models were present, so this is a measured negative --")
        print("  the MPF_FA slot has NO named-family positive control.")
    pf = res.get(("eex", "PF14729"), {})
    print("\n  PF14729 (DUF4467) on ICEBs1 conJ -- the STATED control:")
    if "ICEBs1" in pf:
        for sc, t, i in sorted(pf["ICEBs1"], key=lambda x: -x[0])[:1]:
            print("    FIRES: %s %.1f bits, %d aa at %d..%d"
                  % (t[1], sc, len(aas[i]), t[2], t[3]))
    else:
        print("    DOES NOT FIRE. The RefSeq product string on conJ says")
        print("    'DUF4467 domain-containing protein', but that is RefSeq's own")
        print("    annotation, not our measurement. The control is NOT established.")


if __name__ == "__main__":
    main()
