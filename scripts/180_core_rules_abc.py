#!/usr/bin/env python3
"""STAGE 2 -- apply Rules A/B/C to the conserved short families from script 179.

Script 179 produced families that are CONSERVED, SHORT, not an anchor, not
category 2, and not a CONJScan machine component. That is three of the four terms
of the frozen criterion. The membrane term is missing, and it is the one that
carries the mechanism: entry exclusion works by sitting in the envelope.

    ANCHOR_CANDIDATE = length <= 250 aa
                     AND (lipobox_strict OR tm_count >= 1)

anchor_mode is assigned exactly as in config/anchor_candidate_criterion.yaml:

    lipid        lipobox_strict AND tm == 0     E1  (TrbK, YddJ, Eex_IncN)
    multi_pass   NOT lipobox AND tm >= 2        E1  (TraS, EexS, EexR, ExcA)
    single_pass  NOT lipobox AND tm == 1        E4  no experimental sentinel
    ambiguous    lipobox AND tm >= 1            manual review

READ THE BACKGROUND RATE BEFORE READING ANY PASS RATE. Measured (job 26129299):
21.8% of arbitrary <=250 aa plasmid proteins pass this criterion; the slot rate is
92.2%, enrichment 4.22x. It is a SCREEN, not a classifier. A family passing here
has cleared a bar that roughly one in five random small proteins also clears.
What makes a candidate interesting is passing AND being conserved at a rate its
class's known EEx does not reach -- not passing alone.

The four known-EEx families that script 179 labelled ride along unlabelled in the
scoring and are reported at the end: they are the recall check on THIS stage. If
DUF4467 does not come out lipid-anchored, the topology layer is broken and no
other row is readable.

Orientation is NOT used. The i/o calibration failed (CrcB is dual-topology) and
the ">=30 aa periplasmic C-terminus" criterion was falsified by Marrero & Waldor,
who showed Eex's functional C-terminus is cytoplasmic. Recorded, never read.
"""
import argparse
import collections
import csv
import os
import re
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

WORK = "/global/scratch/users/kh36969/exclusion_gene/recovery/core179"
ANCH = os.path.join(PROJ, "data", "anchors")
STRICT = re.compile(r"[LVI][ASTVIG][GASN]C")
RELAXED = re.compile(r"[ILMFTV][ASTVILMF][GAS]C")
BG_PASS = 21.8          # measured background, job 26129299 -- quote with every rate


def segments(pred):
    out, cur, st = [], None, 0
    for i, c in enumerate(pred):
        k = "TM" if c in "HhBb" else ("S" if c == "S" else ".")
        if k != cur:
            if cur == "TM":
                out.append((st + 1, i))
            cur, st = k, i
    if cur == "TM":
        out.append((st + 1, len(pred)))
    return out


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["plsdb", "ice"], required=True)
    a = ap.parse_args()

    fa = os.path.join(WORK, "tmbed_in_%s.faa" % a.dataset)
    if not os.path.exists(fa):
        print("missing %s -- run script 179 first" % fa); return 1
    n_in = sum(1 for l in open(fa) if l[0] == ">")
    print("representatives to score: %d\n" % n_in, flush=True)

    pred = os.path.join(WORK, "tmbed_out_%s.pred" % a.dataset)
    if not os.path.exists(pred):
        r = subprocess.run(["tmbed", "predict", "-f", fa, "-p", pred,
                            "--out-format", "0", "--model-dir",
                            "/global/home/users/kh36969/.hfcache/prott5",
                            "--use-gpu"], capture_output=True, text=True)
        if r.returncode:
            print("tmbed failed:\n%s" % r.stderr[-800:]); return 1
    L = [l.rstrip("\n") for l in open(pred)]
    got = {}
    for i in range(0, len(L), 3):
        if i + 2 < len(L) and L[i].startswith(">"):
            got[L[i][1:]] = (L[i + 1], L[i + 2])
    print("TMbed returned %d of %d\n" % (len(got), n_in), flush=True)
    if len(got) < n_in:
        print("WARNING: %d representatives missing from the prediction" % (n_in - len(got)))

    rows = []
    for hdr, (seq, p) in got.items():
        f = hdr.split("|")
        cls, ident, rep, n_elem, prev, known = f[0], f[1], f[2], int(f[3]), float(f[4]), f[5]
        tm = len(segments(p))
        lb = bool(STRICT.search(seq[:40]))
        passes = lb or tm >= 1
        mode = ("lipid" if (lb and tm == 0) else
                "single_pass" if (not lb and tm == 1) else
                "multi_pass" if (not lb and tm >= 2) else
                "ambiguous" if lb else "none")
        rows.append({"dataset": a.dataset, "mpf_class": cls, "min_seq_id": ident,
                     "rep": rep, "n_elem": n_elem, "prevalence_pct": prev,
                     "aa": len(seq), "known": known if known != "-" else "",
                     "lipobox_strict": int(lb),
                     "lipobox_relaxed": int(bool(RELAXED.search(seq[:40]))),
                     "tm_count": tm, "anchor_mode": mode, "passes": int(passes)})

    print("=== PASS RATE by class (background on arbitrary small proteins: %.1f%%) ==="
          % BG_PASS)
    print("  %-10s %6s %8s %8s %8s %9s %9s %11s"
          % ("class", "id", "n", "PASS%", "lipid", "single", "multi", "vs bg"))
    for cls in ("T", "F", "I", "FA", "FATA", "G", "B", "C"):
        for ident in ("0.5", "0.7"):
            sel = [r for r in rows if r["mpf_class"] == "MPF_" + cls
                   and r["min_seq_id"] == ident]
            if not sel:
                continue
            n = len(sel)
            pc = 100.0 * sum(r["passes"] for r in sel) / n
            md = collections.Counter(r["anchor_mode"] for r in sel)
            print("  MPF_%-6s %6s %8d %7.1f%% %7d %8d %8d %10.2fx"
                  % (cls, ident, n, pc, md["lipid"], md["single_pass"],
                     md["multi_pass"], pc / BG_PASS))

    print("\n=== RECALL on this stage: the known-EEx families ===")
    kn = [r for r in rows if r["known"]]
    if not kn:
        print("  none present -- this stage has no recall check and its output")
        print("  must not be read as validated.")
    for r in sorted(kn, key=lambda r: (r["known"], -r["prevalence_pct"])):
        print("  %-10s %-8s %-6s %4d aa  lipo=%d tm=%d  %-12s %s"
              % (r["known"], r["mpf_class"], r["min_seq_id"], r["aa"],
                 r["lipobox_strict"], r["tm_count"], r["anchor_mode"],
                 "PASS" if r["passes"] else "*** FAILS THE CRITERION ***"))
    bad = [r for r in kn if not r["passes"]]
    if bad:
        print("\n  %d known EEx family/families FAIL the criterion they were built"
              " from." % len(bad))
        print("  The topology layer is suspect; do not read the candidate table.")

    print("\n=== top passing UNKNOWN families, by prevalence ===")
    for cls in ("G", "B", "C", "FATA", "FA", "T", "F", "I"):
        sel = sorted([r for r in rows if r["mpf_class"] == "MPF_" + cls
                      and r["min_seq_id"] == "0.5" and r["passes"] and not r["known"]],
                     key=lambda r: -r["prevalence_pct"])
        if not sel:
            continue
        print("\n  MPF_%s" % cls)
        print("    %8s %7s %6s %7s %5s %-12s" % ("prev%", "n_elem", "aa", "lipobox",
                                                 "TM", "anchor_mode"))
        for r in sel[:10]:
            print("    %7.1f%% %7d %6d %7d %5d %-12s"
                  % (r["prevalence_pct"], r["n_elem"], r["aa"],
                     r["lipobox_strict"], r["tm_count"], r["anchor_mode"]))

    dest = os.path.join(ANCH, "core_rules_abc_%s.tsv" % a.dataset)
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["mpf_class"], r["min_seq_id"],
                                                -r["prevalence_pct"])))
    print("\nwrote %s (%d rows)" % (dest, len(rows)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
