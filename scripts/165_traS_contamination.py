#!/usr/bin/env python3
"""Is the TraS family contaminated with TraT?

The anomaly: the TraS family is 85.8% lipobox-positive, but the type member --
F's own TraS (P6 in the panel) -- is lipobox-NEGATIVE by all three definitions
and TMbed gives it 4 TM. F TraS is a textbook polytopic inner-membrane protein
and should not carry a lipobox at all.

A family that is 85.8% lipobox-positive whose type member carries none is
internally contradictory.

Most likely cause: traS and traT are ADJACENT genes on F (88,274-88,795 and
88,817-89,551 -- 22 bp apart). If the naming step pulled TraT-type sequences into
the TraS family, then the TraS count is inflated AND the conclusion that plasmids
are TraS-dominant may be partly TraT.

This is the same assignment-confusion class already caught six times in this
project (P1 regex hitting Mycoplasma phage P1, PGAP dropping R100 traS's /gene,
the pED208 sfx collision, the partial SXT record, EexR annotated only as
'hypothetical', anchors leaking into the candidate pool).

Also checks the second anomaly: TIGR04359 (TrbK) is 93.9% lipobox while PF20084
(TrbK_Pfam) is 11.4%. Two models for one family with an 8-fold difference means
they are recognising different things; merging them would dilute the definition.
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
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
STRICT = re.compile(r"[LVI][ASTVIG][GASN]C")


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    rows = [r for r in csv.DictReader(open(SOF), delimiter="\t") if r["slot_occupant_coords"]]
    want = collections.defaultdict(set)
    for r in rows:
        want[r["accession"]].add(r["slot_occupant_coords"])
    seq_of = {}
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
                        seq_of[(acc, k)] = line.rstrip("\n")
    occ = sorted({s for s in seq_of.values()})
    dig = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
           for i, s in enumerate(occ)]

    def scan(acc, ev=None):
        p = os.path.join(HMM, acc + ".hmm")
        if not os.path.exists(p):
            return {}
        with pyhmmer.plan7.HMMFile(p) as fh:
            m = next(iter(fh))
        best = {}
        try:
            it = pyhmmer.hmmsearch([m], dig, cpus=8, bit_cutoffs="gathering")
            for top in it:
                for h in top:
                    i = int(dec(h.name))
                    if i not in best or h.score > best[i]:
                        best[i] = round(h.score, 1)
        except pyhmmer.errors.MissingCutoffs:
            for top in pyhmmer.hmmsearch([m], dig, cpus=8, E=1e-5):
                for h in top:
                    i = int(dec(h.name))
                    if i not in best or h.score > best[i]:
                        best[i] = round(h.score, 1)
        return best

    traS = scan("PF10624")
    traT = scan("PF05818")
    trbk_t = scan("TIGR04359")
    trbk_p = scan("PF20084")

    print("=== 1. IS THE TraS FAMILY CONTAMINATED WITH TraT? ===")
    print("  TraS (PF10624) members            %5d" % len(traS))
    print("  TraT (PF05818) members            %5d" % len(traT))
    both = set(traS) & set(traT)
    print("  hit by BOTH models                %5d  (%.1f%% of TraS)"
          % (len(both), 100.0 * len(both) / max(1, len(traS))))
    lp = [i for i in traS if STRICT.search(occ[i][:40])]
    print("\n  TraS members that are lipobox-positive: %d (%.1f%%)"
          % (len(lp), 100.0 * len(lp) / max(1, len(traS))))
    lp_traT = [i for i in lp if i in traT]
    print("  of those, ALSO hit by PF05818 (TraT)  : %d (%.1f%%)"
          % (len(lp_traT), 100.0 * len(lp_traT) / max(1, len(lp))))
    if both:
        print("\n  score comparison on the dual hits (higher score = better assignment):")
        wins_s = sum(1 for i in both if traS[i] >= traT[i])
        print("    TraS scores higher: %d | TraT scores higher: %d"
              % (wins_s, len(both) - wins_s))
        ex = sorted(both, key=lambda i: traT[i] - traS[i], reverse=True)[:5]
        print("    %-6s %8s %8s %6s  %s" % ("idx", "TraS", "TraT", "aa", "lipobox"))
        for i in ex:
            print("    %-6d %8.1f %8.1f %6d  %s"
                  % (i, traS[i], traT[i], len(occ[i]),
                     "yes" if STRICT.search(occ[i][:40]) else "no"))
    if len(lp_traT) < 0.5 * len(lp):
        print("\n  -> The lipobox-positive TraS members are NOT mostly TraT.")
        print("     The contamination hypothesis does NOT explain the anomaly.")
    else:
        print("\n  -> CONTAMINATION CONFIRMED: most lipobox-positive TraS members")
        print("     are also TraT. The TraS count is inflated.")

    print("\n=== 2. TIGR04359 vs PF20084 -- same family or different things? ===")
    ov = set(trbk_t) & set(trbk_p)
    print("  TIGR04359 (TrbK)     %4d members, lipobox %.1f%%"
          % (len(trbk_t), 100.0 * sum(1 for i in trbk_t if STRICT.search(occ[i][:40]))
             / max(1, len(trbk_t))))
    print("  PF20084  (TrbK_Pfam) %4d members, lipobox %.1f%%"
          % (len(trbk_p), 100.0 * sum(1 for i in trbk_p if STRICT.search(occ[i][:40]))
             / max(1, len(trbk_p))))
    print("  overlap              %4d  (%.1f%% of TIGR, %.1f%% of Pfam)"
          % (len(ov), 100.0 * len(ov) / max(1, len(trbk_t)),
             100.0 * len(ov) / max(1, len(trbk_p))))
    if len(ov) < 0.3 * min(len(trbk_t), len(trbk_p)):
        print("  -> LOW overlap: the two models recognise DIFFERENT protein sets.")
        print("     Keep them as separate subfamilies; do not merge.")
    else:
        print("  -> Substantial overlap; the lipobox difference is within one family.")

    dest = os.path.join(PROJ, "data", "anchors", "traS_contamination.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        for k, v in (("traS_members", len(traS)), ("traT_members", len(traT)),
                     ("both_models", len(both)),
                     ("traS_lipobox_pos", len(lp)), ("traS_lipobox_pos_also_traT", len(lp_traT)),
                     ("trbk_tigr", len(trbk_t)), ("trbk_pfam", len(trbk_p)),
                     ("trbk_overlap", len(ov))):
            w.writerow([k, v])
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
