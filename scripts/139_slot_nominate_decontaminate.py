#!/usr/bin/env python3
"""MPF_T slot iteration, step 1: nominate and decontaminate, on the v1.1 set.

Same code as the original Round 0 path, different input. The nomination set is
the v1.1 `slot_ready__virb5_virb6 = 1` plasmids (6,548), NOT v1's 5,272. The
1,336 `conjscan_only` plasmids the v1.1 disjunction added are divergent IncP-type
systems PF07996 cannot see -- the most likely place for a new family, and never
clustered before.

SLOT DEFINITION. `docs/slot_anchor_is_virb5.md` establishes the rule as
**VirB5 +1 in ALL THREE architectures**, taken in the VirB5/VirB6 strand frame:

    "VirB5 -> eex adjacency is preserved in all three; VirB6 moves from 3' of eex
     to 5' of VirB5. A VirB6-based rule cannot be repaired by a constant offset --
     it needs -1 in canonical and +2 in both inverted groups. The VirB5 rule is
     +1 everywhere."

The instruction for this run said "+1 for canonical/IncI2, -1 for pEC4115". That
-1 was measured relative to VIRB6, and the same file records a correction for
precisely this confusion: an earlier version took the transcription frame from
VirB4's strand, which made one inverted group appear to put eex on the opposite
side of VirB5. So +1 is used as primary AND the -1 occupant is recorded in its own
column, so nothing is lost if the pEC4115 -1 is real in the VirB5 frame.

DECONTAMINATION (§1). phmmer the nominated occupants against the confirmed members
of ALL FOURTEEN anchor families, not just the expected one. In a colinear operon
anchor+1 IS the next anchor, so a positional nomination set is contaminated with
missed anchors by construction -- this has already cost 186 unique sequences
(candidate pool), 190 (recovered set) and a null max inflated to 0.745 (VirB10).

phmmer against confirmed MEMBERS, not the family HMM: the HMM is the thing that
missed these proteins (PF07996 scores 0/1221 of the candidate pool at GA; phmmer
against confirmed members matched 162/228).

DIRECTION ASYMMETRY: only the decontamination direction is used. A false positive
here loses a candidate (safe); a false positive in the recall direction inflates
completeness (unsafe).

KNOWN RESIDUAL CONTAMINATION: VirB7-class proteins (97% undetectable, 48 aa ORFs
the caller often does not emit) cannot be removed and stay in the pool.

The `--anchor` flag exists so the negative slots of §4(a) traverse the IDENTICAL
path including decontamination, not a parallel implementation.
"""
import argparse, collections, csv, json, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import decontaminate as DC

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CAT = os.path.join(PROJ, "data", "release", "v1.1", "catalogue_MPF_T_v1.1.tsv")
OUT = "/global/scratch/users/kh36969/exclusion_gene/slot_v11"


def parse_cell(c):
    """'ga:idx:start..end:strand' -> (idx, start, end, strand)"""
    if not c or c == "none":
        return None
    _, gi, coords, sd = c.split(":")
    s, e = coords.split("..")
    return int(gi), int(s), int(e), 1 if sd == "+" else -1


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchor", default="VirB5")
    ap.add_argument("--offset", type=int, default=1)
    ap.add_argument("--cpus", type=int, default=16)
    ap.add_argument("--tag", default=None)
    a = ap.parse_args()
    tag = a.tag or ("%s%+d" % (a.anchor, a.offset))
    os.makedirs(OUT, exist_ok=True)
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()

    DC.assert_no_evalue_anchor()          # A13 in code

    cat = [r for r in csv.DictReader(open(CAT), delimiter="\t")
           if r["slot_ready__virb5_virb6"] == "1"]
    want = {r["accession"] for r in cat}
    print("nomination set: %d plasmids (v1.1 slot_ready=1)" % len(cat), flush=True)

    seqs, owner, idx, aas = [], [], [], []
    byacc = collections.defaultdict(dict)
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                if acc not in want: continue
                p = line.rstrip("\n")
                byacc[acc][int(f[1])] = len(seqs)
                owner.append(acc); idx.append(int(f[1])); aas.append(p)
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(), sequence=p).digitize(alpha))
    print("cache subset: %d proteins over %d plasmids\n" % (len(seqs), len(byacc)), flush=True)

    nom, minus1, stat = {}, {}, collections.Counter()
    for r in cat:
        acc = r["accession"]
        p = parse_cell(r.get(a.anchor, ""))
        if not p:
            stat["anchor_absent"] += 1; continue
        gi, _, _, sd = p
        for off, store in ((a.offset, nom), (-a.offset, minus1)):
            j = byacc[acc].get(gi + off * sd)
            if j is None:
                if store is nom: stat["no_cds_at_offset"] += 1
                continue
            store[acc] = aas[j]
        if acc in nom: stat["nominated"] += 1
    print("=== nomination (%s, +%d in the anchor's strand frame) ===" % (a.anchor, a.offset))
    for k, v in stat.most_common(): print("  %-20s %d" % (k, v))
    uniq = sorted(set(nom.values()))
    print("  records %d | unique sequences %d" % (len(nom), len(uniq)), flush=True)
    print("  (-%d occupant recorded for %d plasmids, not used downstream)"
          % (a.offset, len(minus1)))

    print("\n=== decontamination vs all %d anchor families ===" % len(DC.ANCHOR_FAMILIES),
          flush=True)
    members = DC.confirmed_members(seqs, aas, cpus=a.cpus)
    print("  confirmed members per family: %s"
          % ", ".join("%s %d" % (k, len(v)) for k, v in sorted(members.items())), flush=True)
    kept, dropped = DC.decontaminate(uniq, members, cpus=a.cpus)
    print("\n  unique in  %d" % len(uniq))
    print("  dropped    %d" % (len(uniq) - len(kept)))
    for fam, s in sorted(dropped.items(), key=lambda x: -len(x[1])):
        print("    %-10s %d" % (fam, len(s)))
    print("  kept       %d unique sequences" % len(kept))
    keptset = set(kept)
    recs = {k: v for k, v in nom.items() if v in keptset}
    print("\n  after decontamination: records %d | unique %d" % (len(recs), len(kept)))
    print("  KNOWN RESIDUAL: VirB7-class (97%% undetectable) cannot be removed.")

    fa = os.path.join(OUT, "slot_%s_clean.faa" % tag)
    with open(fa, "w") as fh:
        for i, s in enumerate(kept):
            fh.write(">cand%d len=%d\n%s\n" % (i, len(s), s))
    with open(os.path.join(OUT, "slot_%s_records.json" % tag), "w") as fh:
        json.dump({"anchor": a.anchor, "offset": a.offset,
                   "nomination_set": len(cat), "records": recs,
                   "minus_offset_occupant": minus1,
                   "dropped_by_family": {k: sorted(v) for k, v in dropped.items()}}, fh)
    print("\nwrote %s (%d sequences)" % (fa, len(kept)))
    L = sorted(len(s) for s in kept)
    if L:
        print("  length: median %d  Q1 %d  Q3 %d  min %d  max %d"
              % (L[len(L)//2], L[len(L)//4], L[3*len(L)//4], L[0], L[-1]))


if __name__ == "__main__":
    main()
