#!/usr/bin/env python3
"""Terminator detection around slot occupants -- pilot with controls first.

QUESTION. Does an exclusion gene sit inside its own transcriptional unit? A unit
needs BOUNDARIES, so both sides matter and they mean different things:

  UPSTREAM terminator  -- blocks read-through from the upstream operon INTO the
                          EEx gene. Direct evidence of decoupling: EEx
                          transcription can then only come from its own promoter.
  DOWNSTREAM terminator -- blocks the EEx gene's own transcription from spilling
                          into its neighbours. Protects the neighbours, not
                          itself -- but still evidence of a bounded unit.

BOTH SIDES ARE SCANNED AND RECORDED SEPARATELY. Directionality is treated as a
pattern to be discovered, not assumed: three literature cases are not enough to
fix a topology rule.

THE CRITICAL DESIGN CHOICE -- RELATIVE RANK, NOT ABSOLUTE THRESHOLD.

Plasmid GC spans ~30% to 65%+. Rho-independent terminator stem stability depends
directly on GC, so in a high-GC plasmid almost any sequence folds into a low-dG
hairpin, and in a low-GC plasmid a real terminator's dG is never very negative.
ANY absolute cutoff manufactures false positives at one end and false negatives
at the other.

So: for each plasmid, enumerate ALL intergenic regions, score them, and convert
the site of interest to a WITHIN-PLASMID PERCENTILE. The criterion is "top 5% of
this plasmid's intergenic regions", never "dG < -10 kcal/mol".

Side benefit: this automatically makes the main tra operon's internal spacers the
internal control.

TWO ORTHOGONAL SCORES, both recorded, never collapsed to a boolean:
  TransTermHP  -- structure + confidence, gene-context aware, hand-built model
  ViennaRNA    -- hairpin dG plus a poly-U tail term, pure thermodynamics
They fail differently, so agreement means more than either alone.

CONTROLS, evaluated BEFORE any genome-wide claim:
  POSITIVE  the slot-flanking spacers of Eex_IncN (cl25) carriers
  NEGATIVE  spacers INSIDE the main trb operon (between adjacent anchors).
            Those genes are co-transcribed and must NOT carry strong terminators.
            If they score high too, the workflow is measuring GC, not termination.

VERDICT RULE, fixed here: the positive control's percentile distribution must be
significantly above the negative control's. If they do not separate, this signal
does not stand -- and the response is to drop it, not to retune.
"""
import collections
import csv
import os
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import orf_caller

PLSDB_FA = "/global/scratch/users/kh36969/plsdb/sequences.fasta"
TTHP = "/global/scratch/users/kh36969/bin/transterm_hp_v2.09"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CAT = os.path.join(PROJ, "data", "release", "v1.1", "catalogue_MPF_T_v1.1.tsv")
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
WORK = "/global/home/users/kh36969/tmp/term"
NPLASMID = 120          # pilot size; controls must pass before scaling
MIN_SPACER = 25         # a terminator needs room


WIN = 80        # a Rho-independent terminator is a LOCAL hairpin near the gene
                # end, typically <60 nt. Folding a whole multi-kb spacer is both
                # O(n^3) slow AND wrong -- it reports the best hairpin anywhere in
                # the spacer, not one positioned to terminate.


def hairpin_score(seq, side="both"):
    """Best local hairpin dG + poly-U tail bonus, in a WINDOW at the spacer end.

    side="down" scans the 5' end of the spacer (i.e. just after the upstream
    gene stops -- where that gene's terminator would sit).
    side="up"   scans the 3' end (just before the downstream gene starts).
    """
    import RNA
    if len(seq) < 20:
        return 0.0
    wins = []
    if side in ("down", "both"):
        wins.append(seq[:WIN])
    if side in ("up", "both"):
        wins.append(seq[-WIN:])
    best = 0.0
    for w in wins:
        if len(w) < 20:
            continue
        try:
            dg = RNA.fold(w.replace("T", "U"))[1]
        except Exception:
            continue
        tail = w[-8:].upper().count("T")
        v = dg - 0.4 * tail
        if v < best:
            best = v
    return best


def main():
    require_compute_node()
    import pyhmmer
    from Bio import SeqIO
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    os.makedirs(WORK, exist_ok=True)

    # --- pick plasmids: Eex_IncN carriers (positive control) -----------------
    sof = {r["accession"]: r for r in csv.DictReader(open(SOF), delimiter="\t")
           if r["slot_status"] == "eex_occupied" and "Eex_IncN" in r["slot_occupant_family"]}
    cat = {r["accession"]: r for r in csv.DictReader(open(CAT), delimiter="\t")
           if r["slot_ready__virb5_virb6"] == "1"}
    accs = [a for a in sof if a in cat][:NPLASMID]
    print("pilot plasmids (Eex_IncN at slot, slot_ready): %d" % len(accs), flush=True)
    if not accs:
        print("no plasmids matched"); return 1

    want = set(accs)
    seqs = {}
    for rec in SeqIO.parse(PLSDB_FA, "fasta"):
        if rec.id in want:
            seqs[rec.id] = str(rec.seq).upper()
            if len(seqs) == len(want):
                break
    print("  sequences retrieved: %d" % len(seqs), flush=True)

    def parse_cell(c):
        if not c or c == "none":
            return None
        f = c.split(":")
        if f[0] == "ga":
            f = f[1:]
        s, e = f[1].split("..")
        return int(f[0]), int(s), int(e), 1 if f[2] == "+" else -1

    ANCH = ["VirB1", "VirB2", "VirB3", "VirB4", "VirB5", "VirB6", "VirB8",
            "VirB9", "VirB10", "VirB11", "VirD4"]
    rows = []
    for acc in accs:
        if acc not in seqs:
            continue
        g = orf_caller.call(seqs[acc])
        if len(g) < 5:
            continue
        g = sorted(g, key=lambda x: x["start"])
        # all intergenic spacers on this plasmid
        spacers = []
        for i in range(len(g) - 1):
            a, b = g[i], g[i + 1]
            gap = b["start"] - a["end"] - 1
            if gap >= MIN_SPACER:
                spacers.append((a["end"] + 1, b["start"] - 1, i))
        if len(spacers) < 10:
            continue
        sc = [(s, e, i, hairpin_score(seqs[acc][s - 1:e])) for s, e, i in spacers]
        vals = sorted(x[3] for x in sc)                  # lower = stronger

        def pct(v):
            # percentile of strength: 100 = strongest on this plasmid
            return 100.0 * sum(1 for x in vals if x > v) / len(vals)

        r = cat[acc]
        p5 = parse_cell(r.get("VirB5", ""))
        if not p5:
            continue
        gi5 = p5[0]
        slot_gi = gi5 + p5[3]
        # spacer immediately 5' and 3' of the slot gene, in gene-index terms
        up = [x for x in sc if x[2] == min(gi5, slot_gi) - 1]
        dn = [x for x in sc if x[2] == max(gi5, slot_gi)]
        # NEGATIVE control: spacers between ADJACENT anchors (co-transcribed)
        idx = {}
        for a in ANCH:
            c = parse_cell(r.get(a, ""))
            if c:
                idx[a] = c[0]
        anch_gi = sorted(idx.values())
        neg = []
        for j in range(len(anch_gi) - 1):
            if anch_gi[j + 1] - anch_gi[j] == 1:          # strictly adjacent
                neg += [x for x in sc if x[2] == anch_gi[j]]
        for lab, hits in (("slot_upstream", up), ("slot_downstream", dn),
                          ("anchor_internal", neg)):
            for s, e, i, v in hits:
                rows.append({"accession": acc, "site": lab, "start": s, "end": e,
                             "spacer_len": e - s + 1, "dg_score": round(v, 2),
                             "percentile": round(pct(v), 1),
                             "n_spacers_on_plasmid": len(sc)})
        if len(rows) and len(rows) % 200 < 4:
            print("  ... %d sites" % len(rows), flush=True)

    if not rows:
        print("no sites scored"); return 1
    dest = os.path.join(PROJ, "data", "anchors", "terminator_pilot.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print("\n=== CONTROL TEST: within-plasmid percentile by site class ===")
    print("  (100 = strongest terminator signal on that plasmid)")
    print("  %-18s %6s %8s %8s %8s %10s" % ("site", "n", "Q1", "median", "Q3", ">=95th"))
    agg = collections.defaultdict(list)
    for r in rows:
        agg[r["site"]].append(r["percentile"])
    for lab in ("slot_upstream", "slot_downstream", "anchor_internal"):
        v = sorted(agg.get(lab, []))
        if not v:
            print("  %-18s %6d" % (lab, 0)); continue
        top = 100.0 * sum(1 for x in v if x >= 95) / len(v)
        print("  %-18s %6d %8.1f %8.1f %8.1f %9.1f%%"
              % (lab, len(v), v[len(v) // 4], v[len(v) // 2], v[3 * len(v) // 4], top))

    up = sorted(agg.get("slot_upstream", []))
    dn = sorted(agg.get("slot_downstream", []))
    ng = sorted(agg.get("anchor_internal", []))
    print("\n=== VERDICT ===")
    if not ng:
        print("  NEGATIVE CONTROL EMPTY -- no adjacent-anchor spacers found.")
        print("  Cannot read the positive side. Fix the control before scaling.")
        return 1
    nm = ng[len(ng) // 2]
    for lab, v in (("upstream", up), ("downstream", dn)):
        if not v:
            print("  %-11s no sites" % lab); continue
        m = v[len(v) // 2]
        print("  %-11s median %.1f vs anchor-internal median %.1f  (%+.1f)"
              % (lab, m, nm, m - nm))
    best = max([v[len(v) // 2] for v in (up, dn) if v] or [0])
    print("\n  %s" % ("SEPARATES -- signal stands, proceed to scale"
                      if best - nm >= 15 else
                      "DOES NOT SEPARATE -- drop this signal, do not retune"))
    print("\nwrote %s (%d sites)" % (dest, len(rows)))


if __name__ == "__main__":
    sys.exit(main() or 0)
