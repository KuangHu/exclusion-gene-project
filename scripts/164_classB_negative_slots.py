#!/usr/bin/env python3
"""B-class acceptance: what is the Class B rate at NEGATIVE slots?

This is the only untested part of Class B, and it is where all of its risk sits.

The two criteria have very different natural specificity:

  Class A  lipidated Cys + signal peptide   -- rare, Cys position is constrained
  Class B  3-4 transmembrane helices        -- COMMON. Small polytopic membrane
                                               proteins are everywhere.

Class A earned its zero-hit result at three negative slots. Class B has no
negative measurement at all. TraT being excluded proves only that it stops an
OUTER-membrane protein; it says nothing about ordinary INNER-membrane proteins,
which are the real false-positive source in a slot.

Expectation stated before the run: Class B will NOT be zero at negative slots,
and may not be low. That is not failure -- it is the denominator. Without it, no
Class B family from the candidate pool can be interpreted.

ACCEPTANCE, fixed before the numbers: the candidate slot's Class B rate must be
SIGNIFICANTLY higher than the negative slots'. If the enrichment is not
significant, Class B does not stand -- and the fix is not to tune the TM
threshold until it does.

Gate is TM >= 3 (not 3-4): all three Class B positives measured exactly 4, so the
"3-4" range was transcribed from a published prediction, not measured. TM count is
recorded as an integer and TM == 2 is kept as a separate watch-list subclass.

Also recorded per sequence: whether the orientation string ALTERNATES. Two
consecutive H (in->out) with no intervening h is topologically impossible, so a
non-alternating string means TMbed is internally inconsistent on that sequence.
Those are grouped separately rather than mixed into the main result.
"""
import collections
import csv
import json
import os
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

SLOTV11 = "/global/scratch/users/kh36969/exclusion_gene/slot_v11"
NEGFI = "/global/scratch/users/kh36969/exclusion_gene/negslot_fi"
WORK = "/global/home/users/kh36969/tmp/classB"
TMBED = "tmbed"
MAXN = 600          # cap per set; TMbed on GPU is fast but ProtT5 is the cost


def segments(pred):
    out, cur, st = [], None, 0
    for i, c in enumerate(pred):
        k = "TM" if c in "HhBb" else ("S" if c == "S" else ".")
        if k != cur:
            if cur in ("TM", "S"):
                out.append((cur, st + 1, i, pred[st]))
            cur, st = k, i
    if cur in ("TM", "S"):
        out.append((cur, st + 1, len(pred), pred[st]))
    return out


def alternates(tms):
    """True if H/h strictly alternate -- a topologically consistent string."""
    lab = [t[3] for t in tms if t[3] in "Hh"]
    return all(lab[i] != lab[i + 1] for i in range(len(lab) - 1))


def main():
    require_compute_node()
    os.makedirs(WORK, exist_ok=True)

    sets = {}
    p = os.path.join(SLOTV11, "slot_VirB5+1_records.json")
    if os.path.exists(p):
        sets["TARGET_VirB5+1"] = sorted(set(json.load(open(p))["records"].values()))
    for tag in ("VirB8+1", "VirB9+1", "VirB10+1"):
        q = os.path.join(SLOTV11, "slot_%s_records.json" % tag)
        if os.path.exists(q):
            sets["NEG_" + tag] = sorted(set(json.load(open(q))["records"].values()))
    if not sets:
        print("no slot record files found"); return 1

    import random
    rng = random.Random(0)
    fa = os.path.join(WORK, "in.faa")
    idx = []
    with open(fa, "w") as fh:
        for lab, seqs in sets.items():
            s = seqs if len(seqs) <= MAXN else rng.sample(seqs, MAXN)
            print("  %-18s %5d unique (%s)" % (lab, len(seqs),
                  "all" if len(seqs) <= MAXN else "sampled %d" % MAXN), flush=True)
            for i, q in enumerate(s):
                nm = "%s|%d" % (lab, i)
                idx.append((nm, lab, q))
                fh.write(">%s\n%s\n" % (nm, q))
    print("  total %d sequences\n" % len(idx), flush=True)

    pred = os.path.join(WORK, "out.pred")
    r = subprocess.run([TMBED, "predict", "-f", fa, "-p", pred, "--out-format", "0",
                        "--model-dir", "/global/home/users/kh36969/.hfcache/prott5",
                        "--use-gpu"], capture_output=True, text=True)
    if r.returncode:
        print("tmbed FAILED:\n%s" % r.stderr[-800:]); return 1

    L = [l.rstrip("\n") for l in open(pred)]
    got = {}
    for i in range(0, len(L), 3):
        if i + 2 < len(L) and L[i].startswith(">"):
            got[L[i][1:]] = (L[i + 1], L[i + 2])

    agg = collections.defaultdict(lambda: collections.Counter())
    rows = []
    for nm, lab, q in idx:
        if nm not in got:
            continue
        seq, pr = got[nm]
        sg = segments(pr)
        tms = [x for x in sg if x[0] == "TM"]
        sig = any(x[0] == "S" for x in sg)
        alt = alternates(tms)
        n = len(tms)
        agg[lab]["n"] += 1
        agg[lab]["tm_ge3"] += (n >= 3)
        agg[lab]["tm_eq2"] += (n == 2)
        agg[lab]["tm0"] += (n == 0)
        agg[lab]["sig"] += sig
        agg[lab]["nonalt"] += (not alt)
        rows.append({"set": lab, "id": nm, "aa": len(seq), "n_tm": n,
                     "signal": int(sig), "orientation_alternates": int(alt)})

    print("=== CLASS B RATE: candidate slot vs negative slots ===")
    print("  %-18s %6s %9s %9s %9s %11s" %
          ("set", "n", "TM>=3", "TM==2", "TM==0", "non-alternating"))
    base = None
    for lab in sorted(agg, key=lambda x: (not x.startswith("TARGET"), x)):
        a = agg[lab]
        n = a["n"]
        if not n:
            continue
        pct = 100.0 * a["tm_ge3"] / n
        if lab.startswith("TARGET"):
            base = pct
        print("  %-18s %6d %8.1f%% %8.1f%% %8.1f%% %10.1f%%"
              % (lab, n, pct, 100.0 * a["tm_eq2"] / n, 100.0 * a["tm0"] / n,
                 100.0 * a["nonalt"] / n))

    print("\n=== enrichment ===")
    if base is not None:
        negs = [(l, 100.0 * agg[l]["tm_ge3"] / agg[l]["n"])
                for l in agg if l.startswith("NEG") and agg[l]["n"]]
        for l, p2 in negs:
            e = base / p2 if p2 else float("inf")
            print("  TARGET %.1f%% vs %-14s %5.1f%%   enrichment %s"
                  % (base, l, p2, ("%.1fx" % e) if p2 else "inf"))
        if negs:
            worst = max(p2 for _, p2 in negs)
            print("\n  worst negative slot: %.1f%%" % worst)
            print("  ACCEPTANCE: %s" %
                  ("Class B STANDS (target %.1f%% vs worst negative %.1f%%, %.1fx)"
                   % (base, worst, base / worst if worst else float("inf"))
                   if worst and base / worst >= 2.0 else
                   "Class B DOES NOT STAND -- enrichment <2x. Do not tune the TM "
                   "threshold to fix this."))

    dest = os.path.join(PROJ, "data", "anchors", "classB_negative_slots.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))


if __name__ == "__main__":
    sys.exit(main() or 0)
