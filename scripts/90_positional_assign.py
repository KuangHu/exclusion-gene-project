#!/usr/bin/env python3
"""Positional CDS assignment for VirB skeleton positions — steps 1 and 5 only.

WHAT THIS IS AND IS NOT
-----------------------
INTERPOLATION between two GA-confirmed anchors, never extrapolation:

    allowed:   [VirB4 HMM✓] [?] [?] [VirB6 HMM✓]    both sides pinned
    forbidden: [?] [?] [VirB6 HMM✓]                 one side only

An unassigned gene between two confirmed anchors is *constrained*; a gene counted
outward from a single anchor is *assumed*. Only the first is done here.

Per the design, only step 1 (length priors) and step 5 (leave-one-out validation)
are implemented. Step 5 is the go/no-go: it measures the error rate of positional
assignment on cases where the HMM answer is known, by hiding that answer. Steps 3
and 4 (actually assigning gaps) are NOT written, and must not be, until step 5's
per-position accuracy is known.

RED LINES (enforced in code, not comments)
  1. positional assignments must never feed slot_status — slot occupancy is
     decided by HMM or length only. Nothing here writes slot_status.
  2. E4_positional never upgrades. There is no promotion path in this module.
  3. exclusion-gene families (TIGR04359, NF033894, NF041429, NF033891, PF10624,
     PF14729) are absent from the anchors AND from the length priors. Their
     presence would make the priors partly a function of the target.

VirB7 IS EXCLUDED. It has no Pfam family at all (pKM101 traN, 48 aa), and its
length band would overlap the 69-76 aa slot occupants this project is looking
for. Assigning that position would be actively harmful.
"""
import argparse
import csv
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import orf_caller
from assertions import AssertionFailure

HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
STEP1A = "/global/scratch/users/kh36969/exclusion_gene/step1a"
PLSDB = "/global/scratch/users/kh36969/plsdb/sequences.fasta"
EXTRA = "/global/scratch/users/kh36969/exclusion_gene/mash/extra_controls.fa"

# VirB position -> (model file, GA-only). VirB7 deliberately absent.
SKELETON = [
    ("VirB1",  "anchor_SLT"),
    ("VirB2",  "anchor_TrbC"),
    ("VirB3",  "anchor_VirB3"),
    ("VirB4",  "anchor_CagE_TrbE_VirB"),
    ("VirB5",  "anchor_T4SS"),
    ("VirB6",  "TrbL"),
    ("VirB8",  "anchor_VirB8"),
    ("VirB9",  "anchor_CagX"),
    ("VirB10", "anchor_TrbI"),
    ("VirB11", "anchor_T2SSE"),
]
FORBIDDEN = {"TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"}


def _check_red_line_3(models):
    """Red line 3, enforced: no exclusion family may reach the priors."""
    for label, hmm, _ in models:
        acc = (hmm.accession or b"").decode() if isinstance(hmm.accession, bytes) \
            else (hmm.accession or "")
        if acc.split(".")[0] in FORBIDDEN:
            raise AssertionFailure(
                "RED LINE 3 VIOLATION: exclusion family %s reached the skeleton "
                "anchors. Length priors built from it would be a function of the "
                "target." % acc)


def load_models():
    import pyhmmer
    out = []
    for label, fn in SKELETON:
        path = os.path.join(HMM, fn + ".hmm")
        with pyhmmer.plan7.HMMFile(path) as fh:
            for hmm in fh:
                ga = None
                try:
                    ga = hmm.cutoffs.gathering[0]
                except Exception:
                    pass
                out.append((label, hmm, ga))
    _check_red_line_3(out)
    return out


def virb4_positive():
    accs = set()
    for fn in sorted(os.listdir(STEP1A)):
        if not fn.startswith("plasmids_"):
            continue
        with open(os.path.join(STEP1A, fn)) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5 or f[0] == "accession":
                    continue
                if "VirB4" in f[4]:
                    accs.add(f[0])
    return accs


def scan_plasmid(seq, models, alpha, pyhmmer):
    """Return calls plus, per gene index, the best GA-level VirB label."""
    calls = orf_caller.call(seq)
    if not calls:
        return calls, {}
    digs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=c["aa"]).digitize(alpha)
            for i, c in enumerate(calls)]
    best = {}
    for label, hmm, ga in models:
        for top in pyhmmer.hmmsearch([hmm], digs, cpus=1, bit_cutoffs="gathering"):
            for h in top:
                nm = h.name
                i = int(nm.decode() if isinstance(nm, bytes) else nm)
                if i not in best or h.score > best[i][1]:
                    best[i] = (label, h.score)
    return calls, best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="0 = all VirB4+ plasmids")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    import pyhmmer
    from Bio import SeqIO
    alpha = pyhmmer.easel.Alphabet.amino()
    models = load_models()
    PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = a.out or os.path.join(PROJ, "data", "positional")
    os.makedirs(out_dir, exist_ok=True)

    want = virb4_positive()
    sys.stderr.write("VirB4+ plasmids available: %d\n" % len(want))
    if a.limit:
        want = set(sorted(want)[:a.limit])
        sys.stderr.write("limited to %d\n" % len(want))

    # ---------------- STEP 1: length priors, GA hits only ----------------
    lengths = {lab: [] for lab, _ in SKELETON}
    skeletons = []          # (accession, [(idx, label|None, aa_len, strand)])
    n = 0
    for path in (PLSDB, EXTRA):
        if not os.path.exists(path):
            continue
        for rec in SeqIO.parse(path, "fasta"):
            if rec.id not in want:
                continue
            calls, best = scan_plasmid(str(rec.seq).upper(), models, alpha, pyhmmer)
            if not calls:
                continue
            row = []
            for i, c in enumerate(calls):
                lab = best.get(i, (None,))[0]
                if lab:
                    lengths[lab].append(c["aa_len"])
                row.append((i, lab, c["aa_len"], c["strand"]))
            skeletons.append((rec.id, row))
            n += 1
            if n % 250 == 0:
                sys.stderr.write("  %d plasmids\n" % n)
    sys.stderr.write("scanned %d plasmids\n" % n)

    prior = {}
    print("\n" + "=" * 78)
    print("STEP 1 -- length priors from GA-level HMM hits only")
    print("=" * 78)
    print("  %-8s %7s %8s %8s %8s %8s %8s   accept band (IQR x1.5)"
          % ("VirB", "n", "min", "Q1", "median", "Q3", "max"))
    for lab, _ in SKELETON:
        v = sorted(lengths[lab])
        if len(v) < 20:
            print("  %-8s %7d   too few instances for a prior" % (lab, len(v)))
            continue
        q1 = v[len(v) // 4]; q3 = v[3 * len(v) // 4]; iqr = q3 - q1
        lo = max(1, int(q1 - 1.5 * iqr)); hi = int(q3 + 1.5 * iqr)
        prior[lab] = (lo, hi)
        print("  %-8s %7d %8d %8d %8d %8d %8d   %d-%d"
              % (lab, len(v), v[0], q1, statistics.median(v), q3, v[-1], lo, hi))
    print("\n  VirB7: EXCLUDED -- no Pfam family, and its length band would overlap")
    print("         the 69-76 aa slot occupants this project is searching for.")

    with open(os.path.join(out_dir, "length_priors.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["virb", "n", "min", "q1", "median", "q3", "max",
                    "accept_lo", "accept_hi"])
        for lab, _ in SKELETON:
            v = sorted(lengths[lab])
            if lab not in prior:
                continue
            w.writerow([lab, len(v), v[0], v[len(v) // 4],
                        int(statistics.median(v)), v[3 * len(v) // 4], v[-1],
                        prior[lab][0], prior[lab][1]])

    # ---------------- STEP 5: leave-one-out validation ----------------
    # Hide one gene's HMM label, then try to recover it from position alone,
    # using ONLY the four gap conditions. This is the go/no-go.
    print("\n" + "=" * 78)
    print("STEP 5 -- leave-one-out: can position alone recover a known label?")
    print("=" * 78)
    order = [lab for lab, _ in SKELETON]
    opos = {lab: i for i, lab in enumerate(order)}
    res = {lab: {"n": 0, "attempted": 0, "correct": 0, "len_ok": 0,
                 "rej_no_flank": 0, "rej_strand": 0, "rej_order": 0,
                 "rej_count": 0, "rej_multi": 0} for lab in order}
    for acc, row in skeletons:
        labelled = [(i, lab, aa, st) for i, lab, aa, st in row if lab]
        for k, (idx, true_lab, aa, st) in enumerate(labelled):
            res[true_lab]["n"] += 1
            # neighbours in the labelled sequence, with the target hidden
            left = labelled[k - 1] if k > 0 else None
            right = labelled[k + 1] if k + 1 < len(labelled) else None
            if left is None or right is None:
                res[true_lab]["rej_no_flank"] += 1
                continue                      # not two-sided: refuse (no extrapolation)
            # condition 4: same strand as both flanks
            if not (left[3] == st == right[3]):
                res[true_lab]["rej_strand"] += 1
                continue
            # condition 2: flanks adjacent-or-near in canonical virB order
            li, ri = opos.get(left[1]), opos.get(right[1])
            if li is None or ri is None or not (ri > li):
                res[true_lab]["rej_order"] += 1
                continue
            missing = order[li + 1:ri]
            # condition 3: gap gene count must equal the number of missing positions
            n_between = right[0] - left[0] - 1
            if len(missing) != 1:
                res[true_lab]["rej_multi"] += 1
                continue
            if n_between != len(missing):
                # THE IMPORTANT REJECTION: extra genes sit in the gap. In MPF_T
                # this is exactly the slot case -- VirB4|VirB5|<occupant>|VirB6
                # gives 2 genes for 1 missing position, so the count check
                # refuses. Correct and conservative, but it means positional
                # assignment can never reach the slot in an OCCUPIED plasmid.
                res[true_lab]["rej_count"] += 1
                continue
            res[true_lab]["attempted"] += 1
            guess = missing[0]
            if guess == true_lab:
                res[true_lab]["correct"] += 1
                lo_hi = prior.get(guess)
                if lo_hi and lo_hi[0] <= aa <= lo_hi[1]:
                    res[true_lab]["len_ok"] += 1

    print("  %-8s %8s %10s %9s %9s   %s"
          % ("VirB", "known", "attempted", "correct", "accuracy", "verdict"))
    print("  " + "-" * 74)
    rows = []
    for lab in order:
        r = res[lab]
        acc_pct = 100.0 * r["correct"] / r["attempted"] if r["attempted"] else None
        if r["attempted"] < 20:
            v = "too few interpolatable cases"
        elif acc_pct >= 95:
            v = "USABLE (>=95%)"
        elif acc_pct >= 85:
            v = "marginal"
        else:
            v = "NOT USABLE"
        print("  %-8s %8d %10d %9d %9s   %s"
              % (lab, r["n"], r["attempted"], r["correct"],
                 ("%.1f%%" % acc_pct) if acc_pct is not None else "-", v))
        rows.append({"virb": lab, "known": r["n"],
                     "rej_no_flank": r["rej_no_flank"], "rej_strand": r["rej_strand"],
                     "rej_order": r["rej_order"], "rej_multi": r["rej_multi"],
                     "rej_count": r["rej_count"], "attempted": r["attempted"],
                     "correct": r["correct"],
                     "accuracy_pct": round(acc_pct, 1) if acc_pct is not None else "",
                     "length_prior_ok": r["len_ok"], "verdict": v})
    with open(os.path.join(out_dir, "loo_validation.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    print("\n  REJECTION BREAKDOWN -- why a known gene was NOT interpolatable")
    print("  %-8s %8s %9s %8s %8s %8s %8s"
          % ("VirB", "known", "no_flank", "strand", "order", "multi", "count"))
    for lab in order:
        r = res[lab]
        print("  %-8s %8d %9d %8d %8d %8d %8d"
              % (lab, r["n"], r["rej_no_flank"], r["rej_strand"], r["rej_order"],
                 r["rej_multi"], r["rej_count"]))
    print("\n  'count' = extra genes in the gap. For VirB5 this is the SLOT case:")
    print("  VirB4|VirB5|<occupant>|VirB6 gives 2 genes for 1 missing position.")

    print("\n  GO/NO-GO: steps 3 and 4 (actual gap assignment) are NOT implemented")
    print("  and must not be until the per-position accuracy above is judged.")
    print("  VirB5 matters most: its accuracy on IncN/IncW (where the HMM works)")
    print("  is the only basis for trusting a positional call on IncP (where")
    print("  PF07996 fails qualitatively -- scores 18/22 vs GA 24.7).")
    print("\nwrote %s" % out_dir)


if __name__ == "__main__":
    main()
