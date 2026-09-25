#!/usr/bin/env python3
"""Is MPF_F's slot-occupant length signal real, or an artefact of the negatives?

Measured: MPF_F target occupants median 177 aa vs negative slots 126 and 132.
BUT the negative slots retain only 12-19% of their nominations after
decontamination, so their survivors are UNUSUAL CASES -- small leftover ORFs at
tightly packed anchor junctions -- not a representative background.

The right control is the one used for MPF_I hydrophobicity: proteins from THE SAME
PLASMIDS at other positions, which is a representative background rather than a
selected remainder.

Reported separately for named (TraS) and unnamed occupants, since the mixed
median (173) conflates two distributions with disjoint IQRs (named 159-165,
unnamed 172-189).
"""
import collections, csv, json, os, statistics, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CAT = os.path.join(PROJ, "data", "release", "v1", "catalogue_MPF_F_v1.tsv")


def band(v):
    v = sorted(v); n = len(v)
    return n, v[n//4], v[n//2], v[3*n//4]


def main():
    require_compute_node()
    cat = [r for r in csv.DictReader(open(CAT), delimiter="\t")
           if r["slot_ready__traG"] == "1" and r["slot_aa_len"]]
    want = {r["accession"] for r in cat}
    slot_len = {r["accession"]: int(r["slot_aa_len"]) for r in cat}
    named = {r["accession"] for r in cat if r["slot_eex_family"]}
    slot_coord = {}
    for r in cat:
        c = r.get("slot_coords") or ""
        slot_coord[r["accession"]] = c

    bg = collections.defaultdict(list)
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                if acc in want: bg[acc].append(len(line.rstrip("\n")))
    print("plasmids %d ; background proteins %d"
          % (len(bg), sum(len(v) for v in bg.values())), flush=True)

    allbg = [L for acc in bg for L in bg[acc]]
    n, q1, med, q3 = band(allbg)
    print("\n=== ALL proteins on the SAME %d plasmids (representative background) ===" % len(bg))
    print("  n=%d  Q1 %d  median %d  Q3 %d" % (n, q1, med, q3))

    for lab, sel in (("slot occupant, NAMED (TraS)", [slot_len[a] for a in named]),
                     ("slot occupant, UNNAMED", [slot_len[a] for a in want - named]),
                     ("slot occupant, ALL", [slot_len[a] for a in want])):
        n, q1, med, q3 = band(sel)
        print("  %-28s n=%5d  Q1 %4d  median %4d  Q3 %4d" % (lab, n, q1, med, q3))

    print("\n=== enrichment: how often is a slot occupant longer than a random")
    print("    protein drawn from its OWN plasmid? (0.50 = no signal) ===")
    import random
    rng = random.Random(0)
    for lab, accs in (("NAMED (TraS)", named), ("UNNAMED", want - named)):
        wins = tot = 0
        for a in accs:
            pool = bg.get(a)
            if not pool or len(pool) < 2: continue
            for _ in range(20):
                if slot_len[a] > rng.choice(pool): wins += 1
                tot += 1
        print("  %-14s P(slot longer than random same-plasmid protein) = %.3f  (n=%d draws)"
              % (lab, wins / max(1, tot), tot))
    print("\n  The negative-slot comparison (177 vs 126/132) is NOT this measure.")
    print("  If P is near 0.50 the slot occupant is unremarkable for its plasmid and")
    print("  the apparent length signal came from the negatives being a selected")
    print("  remainder (12-19%% survival), not from the slot being special.")


if __name__ == "__main__":
    main()
