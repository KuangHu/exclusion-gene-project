#!/usr/bin/env python3
"""Corrected slot-occupant classification: decontamination is part of the DEFINITION.

TWO DEFECTS FIXED.

1. ANCHOR CONTAMINATION. `slot_occupant_families.tsv` (script 148) classified every
   slot occupant using ONLY the six exclusion families and never decontaminated
   against the anchor families. The MPF_T slot iteration DID decontaminate --
   removing 57.9% of nominations, dominated by VirB6 -- so the same project ran two
   paths with different standards. The Pfam scan then found 562 anchor proteins
   (24.5%) sitting in the candidate pool: TrbL 324, VirB8 122, TraG_N 50, SLT 22,
   TrwB_AAD_bind 20, TraG-D_C 20, T4SS-DNA_transf 19, T2SSE 17.

   Decontamination is now part of what "slot occupant" MEANS, not a step some
   analyses happen to run.

2. THE EXCLUSION FAMILY SET WAS INCOMPLETE, in two different ways:

   a. MODEL CHOICE. TrbK was scanned only as TIGR04359. Pfam's own TrbK model
      (PF20084) finds 34 members TIGR04359 misses. Both are now used.

   b. AN ENTIRE CATEGORY WAS ABSENT. The six families covered ENTRY exclusion only.
      SURFACE exclusion was explicitly scoped out -- `data/elements/R27.yaml` records
      eexB as "279 aa, outer membrane, SURFACE exclusion -- out of scope under the
      entry-only cut". But F carries BOTH traS (entry) and traT (surface), IncC both
      eexC and sfx, R27 both eexA and eexB. Half the biology was unreachable by
      construction.

      Added: PF05818 (TraT) and a PrgA/Sea1 profile built from pCF10 prgA + pAD1
      sea1 (both 891 aa, product "surface exclusion protein SEA1/PrgA").

      STILL MISSING and recorded as such: IncC sfx (data/elements/IncC.yaml has
      `sequence: status: NOT_RECOVERED`) and pLS20 ses (absent from the GenBank
      record). Sfx also would not be found by this slot frame even with a model --
      it abuts traN, not VirB6.

Output is a mutually exclusive stratification, so no number double-counts.
"""
import collections, csv, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import decontaminate as DC

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
ENTRY = {"TIGR04359": "TrbK", "NF033894": "Eex_IncN", "NF041429": "EexR",
         "NF033891": "ExcA", "PF10624": "TraS", "PF14729": "DUF4467",
         "PF20084": "TrbK_Pfam"}
SURFACE = {"PF05818": "TraT", "PrgA_Sea1": "PrgA_Sea1"}


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    DC.assert_no_evalue_anchor()

    rows = [r for r in csv.DictReader(open(SOF), delimiter="\t")
            if r["slot_status"] in ("candidate", "eex_occupied") and r["slot_occupant_coords"]]
    want = collections.defaultdict(set)
    for r in rows: want[r["accession"]].add(r["slot_occupant_coords"])
    seq_of = {}
    allseq, allaa = [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]; p = line.rstrip("\n")
                if acc in want:
                    k = "%s..%s" % (f[2], f[3])
                    if k in want[acc]: seq_of[(acc, k)] = p
                # Confirmed anchor members are drawn ONLY from the slot-bearing
                # plasmids, not all 72,438. Scanning 14 HMMs over the full 7.7M
                # cache needs ~5.6 h at 8 cpus and exceeds the wall (job 26104797
                # was killed for this). The restricted set is also the more
                # appropriate reference population: it is the same plasmids the
                # occupants come from.
                if acc in want:
                    allaa.append(p)
                    allseq.append(pyhmmer.easel.TextSequence(
                        name=str(len(allaa)-1).encode(), sequence=p).digitize(alpha))
    occ = collections.defaultdict(set)
    for r in rows:
        s = seq_of.get((r["accession"], r["slot_occupant_coords"]))
        if s: occ[s].add(r["accession"])
    uniq = sorted(occ)
    nrec = sum(len(v) for v in occ.values())
    print("slot occupants: %d records -> %d unique" % (nrec, len(uniq)), flush=True)

    dig = lambda L: [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
                     for i, s in enumerate(L)]
    U = dig(uniq)

    def scan(path, ev=None):
        if not os.path.exists(path): return set()
        with pyhmmer.plan7.HMMFile(path) as fh: m = next(iter(fh))
        kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
        got = set()
        try:
            for top in pyhmmer.hmmsearch([m], U, cpus=8, **kw):
                for h in top: got.add(int(dec(h.name)))
        except pyhmmer.errors.MissingCutoffs:
            for top in pyhmmer.hmmsearch([m], U, cpus=8, E=1e-5):
                for h in top: got.add(int(dec(h.name)))
        return got

    ent = collections.defaultdict(set)
    for acc, nmf in ENTRY.items(): ent[nmf] |= scan(os.path.join(HMM, acc + ".hmm"))
    surf = collections.defaultdict(set)
    for acc, nmf in SURFACE.items(): surf[nmf] |= scan(os.path.join(HMM, acc + ".hmm"), 1e-5)
    E = set().union(*ent.values()) if ent else set()
    S = set().union(*surf.values()) if surf else set()
    print("  entry-exclusion families : %s"
          % ", ".join("%s %d" % (k, len(v)) for k, v in sorted(ent.items()) if v), flush=True)
    print("  SURFACE-exclusion families: %s"
          % ", ".join("%s %d" % (k, len(v)) for k, v in sorted(surf.items()) if v), flush=True)

    print("\n  reference proteins for anchor families: %d (slot-bearing plasmids only)"
          % len(allaa), flush=True)
    print("  decontaminating against the %d anchor families ..." % len(DC.ANCHOR_FAMILIES),
          flush=True)
    members = DC.confirmed_members(allseq, allaa, cpus=8)
    kept, dropped = DC.decontaminate(uniq, members, cpus=8, verbose=False)
    anchor = {i for i, s in enumerate(uniq) if s not in set(kept)}
    print("  anchor-family occupants removed: %d of %d unique (%.1f%%)"
          % (len(anchor), len(uniq), 100.0*len(anchor)/len(uniq)))
    for f, s in sorted(dropped.items(), key=lambda x: -len(x[1]))[:8]:
        print("      %-12s %d" % (f, len(s)))

    # mutually exclusive, in priority order
    strat = collections.OrderedDict()
    strat["entry exclusion (named)"] = E - S
    strat["surface exclusion (named)"] = S - E
    strat["both entry+surface hit"] = E & S
    assigned = E | S
    strat["ANCHOR family (was miscounted as candidate)"] = anchor - assigned
    assigned |= anchor
    strat["unassigned candidate"] = set(range(len(uniq))) - assigned
    print("\n=== MUTUALLY EXCLUSIVE STRATIFICATION (unique sequences) ===")
    print("  %-46s %7s %8s %9s" % ("class", "unique", "pct", "records"))
    for k, v in strat.items():
        rec = sum(len(occ[uniq[i]]) for i in v)
        print("  %-46s %7d %7.1f%% %9d" % (k, len(v), 100.0*len(v)/len(uniq), rec))
    print("  %-46s %7d %7.1f%%" % ("TOTAL", sum(len(v) for v in strat.values()), 100.0))

    dest = os.path.join(PROJ, "data", "anchors", "slot_strata_corrected.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["class", "n_unique", "n_records"])
        for k, v in strat.items():
            w.writerow([k, len(v), sum(len(occ[uniq[i]]) for i in v)])
    print("\nwrote %s" % dest)
    print("\n  NOT COVERED, recorded: IncC sfx (sequence NOT_RECOVERED) and pLS20 ses")
    print("  (absent from the GenBank record). sfx abuts traN, so this slot frame")
    print("  would miss it even with a model.")


if __name__ == "__main__":
    main()
