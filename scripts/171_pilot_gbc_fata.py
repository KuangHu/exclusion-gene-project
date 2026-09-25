#!/usr/bin/env python3
"""PILOT -- can the four remaining MPF classes be detected at all?

MPF_T, F, I and the FA/FATA pool are catalogued. The remaining classes are G, B,
C and FATA-as-its-own-class. Before any database-scale run, one question has to
be answered on known material: DO THE PROFILES FIRE ON THE PROTOTYPE, AND ONLY
THERE?

This is a feasibility probe, not a measurement. It produces no counts over PLSDB
and nothing here may be quoted as a result. The database scan is separately
blocked anyway -- 10 of 34 cds_cache_full shards are on the dead Lustre OST and a
scan on the surviving 24 would silently under-count every class.

WHAT MAKES THIS READABLE: all eight classes are scanned against every seed, so
each seed carries its own negative controls. CONJScan's mandatory genes are
IDENTICAL across all eight classes (virb4/traU + t4cp1/2 + one of ten MOB
relaxases), so class identity lives ENTIRELY in the accessory profiles. That is
the thing being tested.

  PERFECT-SCORE BASELINE: pCF10 is the published MPF_FATA reference and is
  already a seed in this project (script 133 used it as the FATA control). If
  FATA profiles do not top the table on pCF10, the scan is broken and NO other
  row may be read. Same stop condition as script 133.

Seeds, each the prototype named in the Guglielmini 2014 classification:

    MPF_G     AJ627386   ICEHin1056, H. influenzae (Tfc1-Tfc24)
    MPF_B     AF289050   CTnDOT transfer region, B. thetaiotaomicron
    MPF_C     BA000020   pCC7120alpha, Nostoc sp. PCC 7120
    MPF_FATA  pCF10      AY855841.2 -- the baseline

AF289050 is a transfer REGION, not a whole element, so its absolute profile count
is not comparable to the others. Only the within-seed ranking is.

Proteins come from the GenBank CDS translations, not from the ORF caller, so this
has no dependency on the blocked cache.
"""
import collections
import csv
import os
import sys
import xml.etree.ElementTree as ET

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CJ = os.path.join(PROJ, "data", "CONJScan", "profiles")
DEFS = os.path.join(PROJ, "data", "CONJScan", "definitions", "Plasmids")
PILOT = os.path.join(PROJ, "data", "seed", "pilot_mpf")
GB = os.path.join(PROJ, "data", "seed", "genbank")

SEEDS = [("ICEHin1056", os.path.join(PILOT, "AJ627386.gb"), "G"),
         ("CTnDOT",     os.path.join(PILOT, "AF289050.gb"), "B"),
         ("pCC7120a",   os.path.join(PILOT, "BA000020.gb"), "C"),
         ("pCF10",      os.path.join(GB, "pCF10__AY855841.2.gb"), "FATA")]
BASELINE = ("pCF10", "FATA")
CLASSES = ["T", "F", "I", "FA", "FATA", "G", "B", "C"]
# shared by every class definition -- cannot discriminate, so excluded from the
# accessory sets and reported separately
GENERIC = ("T4SS_virb4", "T4SS_I_traU", "T4SS_t4cp1", "T4SS_t4cp2", "T4SS_tcpA")


def accessory(cls):
    """Class-specific profiles: everything in the definition minus the generics."""
    root = ET.parse(os.path.join(DEFS, "T4SS_type%s.xml" % cls)).getroot()
    names = {g.get("name") for g in root.iter("gene")}
    return sorted(n for n in names
                  if n not in GENERIC and not n.startswith("T4SS_MOB"))


def main():
    require_compute_node()
    import pyhmmer
    from Bio import SeqIO
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    acc_sets = {c: accessory(c) for c in CLASSES}
    for c in CLASSES:
        print("  MPF_%-5s %2d accessory profiles" % (c, len(acc_sets[c])))
    # a profile shared by two classes would make the ranking meaningless
    owner = collections.Counter(p for c in CLASSES for p in acc_sets[c])
    shared = {p: n for p, n in owner.items() if n > 1}
    print("\n  profiles claimed by more than one class: %d %s"
          % (len(shared), sorted(shared) if shared else ""))

    seqs, tag = [], []
    for name, path, _ in SEEDS:
        if not os.path.exists(path):
            raise SystemExit("seed missing: %s" % path)
        n = 0
        for rec in SeqIO.parse(path, "genbank"):
            for f in rec.features:
                if f.type != "CDS":
                    continue
                p = f.qualifiers.get("translation", [None])[0]
                if not p:
                    continue
                seqs.append(pyhmmer.easel.TextSequence(
                    name=("%d" % len(tag)).encode(), sequence=p).digitize(alpha))
                tag.append(name)
                n += 1
        print("  %-12s %4d CDS with translations" % (name, n))
        if n == 0:
            raise SystemExit("%s yielded no proteins -- refusing to score it" % name)

    # ---- scan every profile of every class against every seed ---------------
    fired = collections.defaultdict(set)      # (seed, class) -> profiles
    missing = []
    allp = sorted({p for c in CLASSES for p in acc_sets[c]} | set(GENERIC))
    for prof in allp:
        path = os.path.join(CJ, prof + ".hmm")
        if not os.path.exists(path):
            missing.append(prof)
            continue
        with pyhmmer.plan7.HMMFile(path) as fh:
            hmm = next(iter(fh))
        for top in pyhmmer.hmmsearch([hmm], seqs, cpus=4, bit_cutoffs="gathering"):
            for h in top:
                s = tag[int(dec(h.name))]
                for c in CLASSES:
                    if prof in acc_sets[c]:
                        fired[(s, c)].add(prof)
                if prof in GENERIC:
                    fired[(s, "generic")].add(prof)
    if missing:
        raise SystemExit("profiles missing from disk: %s\n"
                         "A silent skip here would report a class as absent."
                         % missing)

    # ---- the baseline, before anything else is read -------------------------
    bs, bc = BASELINE
    top = max(CLASSES, key=lambda c: len(fired[(bs, c)]))
    print("\n=== BASELINE: %s must top out on MPF_%s ===" % (bs, bc))
    for c in CLASSES:
        print("    MPF_%-5s %2d" % (c, len(fired[(bs, c)])))
    if top != bc or not fired[(bs, bc)]:
        print("\n  BASELINE FAILED (top = MPF_%s). The scan is broken." % top)
        print("  STOPPING -- no other row may be read.")
        return 1
    print("  BASELINE PASSED: top = MPF_%s with %d profiles.\n"
          % (bc, len(fired[(bs, bc)])))

    # ---- the table ----------------------------------------------------------
    print("=== accessory profiles firing, seed x class ===")
    print("  %-12s %s  %s" % ("seed", " ".join("%5s" % ("MPF_" + c) for c in CLASSES), "generic"))
    rows = []
    for name, _, own in SEEDS:
        cells = " ".join("%5d" % len(fired[(name, c)]) for c in CLASSES)
        print("  %-12s %s  %5d" % (name, cells, len(fired[(name, "generic")])))
        rows.append({"seed": name, "expected_class": own,
                     "generic": len(fired[(name, "generic")]),
                     **{("MPF_" + c): len(fired[(name, c)]) for c in CLASSES},
                     "own_profiles": ";".join(sorted(fired[(name, own)]))})

    print("\n=== VERDICT: does each seed top out on its OWN class? ===")
    ok = True
    for name, _, own in SEEDS:
        best = max(CLASSES, key=lambda c: len(fired[(name, c)]))
        n_own, n_best = len(fired[(name, own)]), len(fired[(name, best)])
        good = (best == own and n_own > 0)
        ok &= good
        print("  %-12s expected MPF_%-5s got MPF_%-5s (%d vs %d)  %s"
              % (name, own, best, n_own, n_best, "PASS" if good else "FAIL"))
    print("\n  %s" % ("All four prototypes are separated by their accessory profiles."
                      " Class detection is FEASIBLE; the entry criterion still has to"
                      " be built against cross-class backgrounds." if ok else
                      "At least one prototype is NOT separated. Do not build an entry"
                      " criterion on these profiles until this is understood."))
    for name, _, own in SEEDS:
        if fired[(name, own)]:
            print("\n  %s / MPF_%s fired: %s"
                  % (name, own, ", ".join(sorted(fired[(name, own)]))))

    dest = os.path.join(PROJ, "data", "anchors", "pilot_gbc_fata.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
