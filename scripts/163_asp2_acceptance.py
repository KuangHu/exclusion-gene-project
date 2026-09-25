#!/usr/bin/env python3
"""A-class acceptance: does Asp@+2 hold across the 842 named sequences?

The panel gave 2 testable positives (TrbK, Eex_pKM101 -- both Asp) against 1
negative (TraT -- Gly). 2-vs-1 is not a rule, only a consistent direction. The
rule's credibility comes from being taken from independent literature rather than
fitted to those three, but its RECALL is untested and has two known gaps:

  SPECIES RANGE. The +2 Asp Lol-avoidance rule is established in E. coli /
  Enterobacteriaceae. Other lineages sort differently, and PLSDB spans many
  non-enterobacterial hosts. A hard Asp@+2 gate would lose true positives
  SYSTEMATICALLY BY HOST TAXON -- a biased loss, which would then bias any
  host-distribution conclusion drawn from the families.

  GRAM-POSITIVES. No outer membrane, so nothing to sort against and the rule does
  not apply at all. On Gram-positive ICEs -- the top target after ICEBs1 -- the
  A-class repair simply does not exist.

So this measures the pass rate per family AND per host taxon. TrbK and Eex_IncN
are known true positives: if they pass at 90%+ the cross-species loss is
tolerable and Asp@+2 can gate; at ~60% the loss is quantified and the rule must
be demoted to a scoring field, not a gate.

Mature position 1 is the Cys; position 2 is the residue IMMEDIATELY after it.
This project has already made that off-by-one once (read position 3 instead of
mature position 2), so it is stated explicitly here.
"""
import collections
import csv
import os
import re
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
PLSDB = "/global/scratch/users/kh36969/plsdb"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
FAM = {"TIGR04359": "TrbK", "PF20084": "TrbK_Pfam", "NF033894": "Eex_IncN",
       "NF041429": "EexR", "NF033891": "ExcA", "PF10624": "TraS",
       "PF14729": "DUF4467", "PF05818": "TraT", "PrgA_Sea1": "PrgA_Sea1"}
STRICT = re.compile(r"[LVI][ASTVIG][GASN]C")
RELAXED = re.compile(r"[LVIMFWY][ASTVIGN][GASNDQ]C")
ENTERO = {"Escherichia", "Klebsiella", "Salmonella", "Enterobacter", "Citrobacter",
          "Shigella", "Serratia", "Proteus", "Yersinia", "Cronobacter", "Raoultella",
          "Providencia", "Morganella", "Pantoea", "Edwardsiella", "Kluyvera"}


def plus2(s):
    """Return (cys_pos, mature+2 residue) or (None, None)."""
    m = STRICT.search(s[:40]) or RELAXED.search(s[:40])
    if not m:
        return None, None
    c = m.end() - 1                      # 0-based Cys
    return c + 1, (s[c + 1] if c + 1 < len(s) else None)


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    host = {}
    tp = os.path.join(PLSDB, "taxonomy.csv")
    if os.path.exists(tp):
        csv.field_size_limit(10 ** 7)
        for r in csv.DictReader(open(tp)):
            g = r.get("genus") or r.get("GENUS") or ""
            a = r.get("NUCCORE_ACC") or r.get("accession") or ""
            if a:
                host[a] = g.split()[0] if g else ""

    rows = [r for r in csv.DictReader(open(SOF), delimiter="\t") if r["slot_occupant_coords"]]
    want = collections.defaultdict(set)
    for r in rows:
        want[r["accession"]].add(r["slot_occupant_coords"])
    seq_of, acc_of = {}, collections.defaultdict(set)
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                if acc in want:
                    k = "%s..%s" % (f[2], f[3])
                    if k in want[acc]:
                        p = line.rstrip("\n")
                        seq_of[(acc, k)] = p
                        acc_of[p].add(acc)
    for r in rows:
        pass
    occ = sorted({s for s in seq_of.values()})
    print("slot occupants: %d unique" % len(occ), flush=True)

    dig = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
           for i, s in enumerate(occ)]
    byfam = collections.defaultdict(set)
    for acc, lab in FAM.items():
        p = os.path.join(HMM, acc + ".hmm")
        if not os.path.exists(p):
            continue
        with pyhmmer.plan7.HMMFile(p) as fh:
            m = next(iter(fh))
        try:
            for top in pyhmmer.hmmsearch([m], dig, cpus=8, bit_cutoffs="gathering"):
                for h in top:
                    byfam[lab].add(int(dec(h.name)))
        except pyhmmer.errors.MissingCutoffs:
            for top in pyhmmer.hmmsearch([m], dig, cpus=8, E=1e-5):
                for h in top:
                    byfam[lab].add(int(dec(h.name)))

    print("\n=== Asp@+2 pass rate per named family ===")
    print("  %-11s %6s %9s %9s %10s  %s" %
          ("family", "n", "lipobox", "Asp@+2", "of lipobox", "verdict"))
    out = []
    for lab in sorted(byfam, key=lambda x: -len(byfam[x])):
        idx = sorted(byfam[lab])
        if not idx:
            continue
        S = [occ[i] for i in idx]
        lb = [x for x in S if plus2(x)[0]]
        asp = [x for x in lb if plus2(x)[1] == "D"]
        n = len(S)
        pct = 100.0 * len(asp) / len(lb) if lb else float("nan")
        v = ("gate OK (>=90%)" if pct >= 90 else
             "DEMOTE to score (<90%)" if pct == pct else "no lipobox")
        print("  %-11s %6d %8.1f%% %9d %9.1f%%  %s"
              % (lab, n, 100.0 * len(lb) / n, len(asp), pct, v))
        out.append({"family": lab, "n": n, "lipobox_pct": round(100.0 * len(lb) / n, 1),
                    "n_asp": len(asp),
                    "asp_of_lipobox_pct": round(pct, 1) if pct == pct else ""})

    print("\n=== residue at mature +2, across all lipobox-positive named sequences ===")
    allr = collections.Counter()
    for lab in byfam:
        for i in byfam[lab]:
            c, r2 = plus2(occ[i])
            if c and r2:
                allr[r2] += 1
    tot = sum(allr.values())
    for r2, c in allr.most_common(8):
        print("  %s  %5d  %5.1f%%%s" % (r2, c, 100.0 * c / tot,
                                        "   <- Lol avoidance (inner membrane)" if r2 == "D" else ""))

    print("\n=== by host taxon: is the loss biased? ===")
    ent = collections.Counter(); oth = collections.Counter()
    for lab in byfam:
        for i in byfam[lab]:
            c, r2 = plus2(occ[i])
            if not (c and r2):
                continue
            gen = ""
            for a in acc_of[occ[i]]:
                gen = host.get(a, "") or gen
            (ent if gen in ENTERO else oth)[r2 == "D"] += 1
    for lab, d in (("Enterobacteriaceae", ent), ("other / unknown", oth)):
        t = sum(d.values())
        if t:
            print("  %-20s n=%5d  Asp@+2 %5.1f%%" % (lab, t, 100.0 * d[True] / t))
    print("\n  If the two rates differ materially, a hard Asp@+2 gate loses true")
    print("  positives BY HOST TAXON, which would bias any host-distribution claim.")

    dest = os.path.join(PROJ, "data", "anchors", "asp2_acceptance.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
