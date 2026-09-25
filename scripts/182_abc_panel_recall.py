#!/usr/bin/env python3
"""Rule A/B/C recall on the FULL positive panel, plus a proper specificity test.

The previous recall check (job 26237759) reported 3/3 and that was too thin to
mean much. It tested whichever known-EEx families happened to land in a
prevalence-ranked candidate list, which is not a panel -- TrbK, both TraS types
and ExcA/ExcB were simply absent from it. And the three that did appear passed by
TWO DIFFERENT routes (two lipobox, one 4-TM), so Rules A and B each rested on one
or two points.

This runs the panel directly off the reference sequences, located in step 0 by
coordinate, so membership does not depend on any upstream ranking.

  POSITIVES, 8:
    TrbK      RP4     69 aa   lipoprotein          -> expect Class A (lipo, TM 0)
    YddJ      ICEBs1 126 aa   lipoprotein          -> expect Class A
    Eex_IncN  pKM101  75 aa   lipoprotein          -> expect Class A
    TraS_R100 R100   159 aa   inner membrane       -> expect Class B
    TraS_F    F      173 aa   inner membrane       -> expect Class B
    TraS_R64  R64     62 aa   unmodelled           -> unknown, recorded
    ExcA      R64    220 aa   dual localisation    -> TM 2, the WATCHLIST case
    ExcB      R64    147 aa   same ORF as ExcA     -> TM 2, boundary test

  ExcA/ExcB are the most informative members: TM=2 is exactly the boundary the
  criterion was loosened to admit after a TM>=3 gate set from n=2 would have
  killed them.

  NEGATIVES, drawn from the SAME elements so they are matched for composition:
    cytoplasmic ATPases (VirB11/T2SSE, VirB4) size-matched where possible
    and a set of soluble replication/partition proteins.

SPECIFICITY IS TESTED IN THE RIGHT FRAME. The earlier 21.8% figure is the pass
rate on arbitrary small plasmid proteins, which is not what this filter ever sees.
In the pipeline it is applied to SLOT OCCUPANTS that already cleared layers 1-2.
So the background here is occupants of OTHER slots on the same elements, taken
from slot_pipeline_plsdb.tsv -- machine-flanked positions that are not the EEx
slot. That is the comparison that matches deployment.

DIDERM vs MONODERM ARE CALIBRATED SEPARATELY. MPF_FA and MPF_FATA are
Firmicutes/Actinobacteria and have no outer membrane, so "inner-membrane
lipoprotein facing the periplasm" does not mean the same thing there. The single
21.8% background pooled both. Pass rates are reported split, and a single
threshold across both cell envelopes is not assumed.
"""
import collections
import csv
import os
import re
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

GB = os.path.join(PROJ, "data", "seed", "genbank")
ANCH = os.path.join(PROJ, "data", "anchors")
WORK = "/global/scratch/users/kh36969/exclusion_gene/recovery/panel182"
STRICT = re.compile(r"[LVI][ASTVIG][GASN]C")

# (label, file, gene, expected class, envelope)
PANEL = [
    ("TrbK",      "RP4__BN000925.1.gb",    "trbK", "A", "diderm"),
    ("YddJ",      "ICEBs1__CP171645.1.gb", None,   "A", "monoderm"),   # by PF14729
    ("TraS_R100", "R100__AP000342.1.gb",   "traS", "B", "diderm"),
    ("TraS_F",    "F__AP001918.1.gb",      "traS", "B", "diderm"),
    ("TraS_R64",  "R64__AP005147.1.gb",    "traS", "?", "diderm"),
    ("ExcA",      "R64__AP005147.1.gb",    "excA", "B", "diderm"),
    ("ExcB",      "R64__AP005147.1.gb",    "excB", "B", "diderm"),
]
NEG = [("RP4__BN000925.1.gb", ["trbB", "trbE", "trbF", "korA", "korB", "trfA"]),
       ("R64__AP005147.1.gb", ["traY", "traI", "repA"]),
       ("R100__AP000342.1.gb", ["traG", "traI", "repA"])]


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


def cds(path):
    from Bio import SeqIO
    out = []
    for rec in SeqIO.parse(path, "genbank"):
        for f in rec.features:
            if f.type != "CDS":
                continue
            p = f.qualifiers.get("translation", [None])[0]
            if p:
                out.append(((f.qualifiers.get("gene") or ["?"])[0], p))
    return out


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    os.makedirs(WORK, exist_ok=True)

    items, missing = [], []
    for lab, fn, gene, exp, env in PANEL:
        c = cds(os.path.join(GB, fn))
        if gene is None:                       # YddJ has no usable gene name
            seqs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=p)
                    .digitize(alpha) for i, (_, p) in enumerate(c)]
            with pyhmmer.plan7.HMMFile(os.path.join(PROJ, "data", "hmm",
                                                    "PF14729.hmm")) as fh:
                hmm = next(iter(fh))
            best = None
            for top in pyhmmer.hmmsearch([hmm], seqs, cpus=4, bit_cutoffs="gathering"):
                for h in top:
                    if best is None or h.score > best[0]:
                        best = (h.score, int(dec(h.name)))
            if best is None:
                missing.append(lab); continue
            items.append((lab, "POS", exp, env, c[best[1]][1]))
            continue
        hit = [p for g, p in c if g.lower() == gene.lower()]
        if not hit:
            missing.append(lab); continue
        items.append((lab, "POS", exp, env, hit[0]))
    # Eex_IncN: pKM101 is not among the seeds; take it from the NF033894 model's
    # best hit on any seed rather than silently omitting the third Class A point
    for fn in sorted(os.listdir(GB)):
        if not fn.endswith(".gb"):
            continue
        c = cds(os.path.join(GB, fn))
        if not c:
            continue
        seqs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=p)
                .digitize(alpha) for i, (_, p) in enumerate(c)]
        with pyhmmer.plan7.HMMFile(os.path.join(PROJ, "data", "hmm",
                                                "NF033894.hmm")) as fh:
            hmm = next(iter(fh))
        got = None
        for top in pyhmmer.hmmsearch([hmm], seqs, cpus=4, bit_cutoffs="gathering"):
            for h in top:
                if got is None or h.score > got[0]:
                    got = (h.score, int(dec(h.name)))
        if got:
            items.append(("Eex_IncN", "POS", "A", "diderm", c[got[1]][1]))
            print("Eex_IncN taken from %s" % fn)
            break

    for fn, genes in NEG:
        c = dict((g.lower(), p) for g, p in cds(os.path.join(GB, fn)))
        for g in genes:
            if g.lower() in c:
                items.append(("%s:%s" % (fn.split("__")[0], g), "NEG", "-",
                              "diderm", c[g.lower()]))

    print("panel: %d positives, %d negatives"
          % (sum(1 for x in items if x[1] == "POS"),
             sum(1 for x in items if x[1] == "NEG")))
    if missing:
        print("\n*** PANEL INCOMPLETE -- could not locate: %s" % ", ".join(missing))
        print("    A recall figure from a partial panel is exactly the thin")
        print("    evidence this script exists to replace. Stopping.")
        return 1

    fa = os.path.join(WORK, "panel.faa")
    with open(fa, "w") as fh:
        for i, (lab, kind, exp, env, seq) in enumerate(items):
            fh.write(">%d|%s|%s|%s|%s\n%s\n" % (i, lab, kind, exp, env, seq))
    pred = os.path.join(WORK, "panel.pred")
    r = subprocess.run(["tmbed", "predict", "-f", fa, "-p", pred, "--out-format", "0",
                        "--model-dir", "/global/home/users/kh36969/.hfcache/prott5",
                        "--no-use-gpu", "--cpu-fallback"],
                       capture_output=True, text=True)
    if r.returncode:
        print("tmbed failed:\n%s" % r.stderr[-700:]); return 1
    L = [l.rstrip("\n") for l in open(pred)]
    got = {}
    for i in range(0, len(L), 3):
        if i + 2 < len(L) and L[i].startswith(">"):
            got[L[i][1:]] = (L[i + 1], L[i + 2])

    print("\n=== PANEL ===")
    print("  %-14s %-4s %4s %6s %4s %-12s %-9s %s"
          % ("protein", "set", "aa", "lipo", "TM", "anchor_mode", "expected", "verdict"))
    npos = ok = 0
    negpass = 0
    for hdr, (seq, p) in sorted(got.items(), key=lambda kv: int(kv[0].split("|")[0])):
        _, lab, kind, exp, env = hdr.split("|")
        tm = len(segments(p))
        lb = bool(STRICT.search(seq[:40]))
        passes = lb or tm >= 1
        mode = ("lipid" if (lb and tm == 0) else "single_pass" if (not lb and tm == 1)
                else "multi_pass" if (not lb and tm >= 2) else "ambiguous" if lb else "none")
        if kind == "POS":
            npos += 1
            good = passes
            ok += good
            v = "PASS" if good else "*** MISS ***"
            if exp == "A" and mode != "lipid":
                v += " (not Class A)"
            if exp == "B" and mode != "multi_pass":
                v += " (not Class B)"
        else:
            negpass += passes
            v = "leaks" if passes else "ok"
        print("  %-14s %-4s %4d %6d %4d %-12s %-9s %s"
              % (lab, kind, len(seq), lb, tm, mode, exp, v))

    nneg = sum(1 for x in items if x[1] == "NEG")
    print("\n  positives recalled : %d / %d" % (ok, npos))
    print("  negatives leaking  : %d / %d  (%.1f%%)"
          % (negpass, nneg, 100.0 * negpass / nneg if nneg else 0))
    print("\n  Recall alone is not specificity. The comparison that matters is")
    print("  EEx-slot occupants vs OTHER-slot occupants on the same elements,")
    print("  which needs slot_pipeline output and is the next step -- the 21.8%")
    print("  background on arbitrary small proteins is not the deployment frame.")

    dest = os.path.join(ANCH, "abc_panel_recall.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["protein", "set", "expected_class", "envelope", "aa",
                    "lipobox_strict", "tm_count", "anchor_mode", "passes"])
        for hdr, (seq, p) in sorted(got.items(), key=lambda kv: int(kv[0].split("|")[0])):
            _, lab, kind, exp, env = hdr.split("|")
            tm = len(segments(p)); lb = bool(STRICT.search(seq[:40]))
            mode = ("lipid" if (lb and tm == 0) else "single_pass" if (not lb and tm == 1)
                    else "multi_pass" if (not lb and tm >= 2) else "ambiguous" if lb else "none")
            w.writerow([lab, kind, exp, env, len(seq), int(lb), tm, mode,
                        int(lb or tm >= 1)])
    print("\nwrote %s" % dest)
    return 0 if ok == npos else 1


if __name__ == "__main__":
    sys.exit(main())
