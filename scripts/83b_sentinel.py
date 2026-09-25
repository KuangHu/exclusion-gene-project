#!/usr/bin/env python3
"""v2 §10 sentinel for step 1a: the controls must score perfectly, or stop.

Runs the exact step-1a logic on the four PLSDB-resident controls. Each must be
VirB4-positive (entry criterion) AND VirB6-positive (slot boundary) AND
VirB5-positive (empty-slot test). If any fails, the library run is invalid and
must not be read.
"""
import os, subprocess, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import orf_caller
from Bio import SeqIO
import pyhmmer

HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
GB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                  "data", "seed", "genbank")
CONTROLS = [("RP4", "RP4__BN000925.1.gb"), ("R751", "R751__NC_001735.4.gb"),
            ("R388", "R388__NC_028464.1.gb"), ("R46", "R46__AY046276.1.gb")]

# NEGATIVE controls. A sentinel made only of positive assertions cannot catch a
# bug that calls EVERYTHING positive -- which is precisely the failure the \r
# bug produced (146/146 plasmids scored as anchor-bearing). These must come out
# VirB4-NEGATIVE.
NEGATIVES = [
    ("ColE1", "ColE1__J01566.1.gb",
     "6.6 kb mobilisable ColE1: no T4SS at all, only the mob relaxase machinery"),
    ("pKM101_frag", "pKM101__U09868.1.gb", None),   # see note below
]
ANCHORS = [("VirB4", "anchor_CagE_TrbE_VirB", "ga", 0),
           ("VirB6", "TrbL", "ga", 0),
           ("VirB5", "anchor_T4SS", "e1e-5", 150)]

alpha = pyhmmer.easel.Alphabet.amino()
models = []
for lab, fn, mode, minaa in ANCHORS:
    with pyhmmer.plan7.HMMFile(os.path.join(HMM, fn + ".hmm")) as fh:
        for h in fh:
            models.append((lab, h, mode, minaa))

print("v2 step-1a sentinel -- controls must all pass\n")
print("  %-7s %-8s %-8s %-8s  verdict" % ("control", "VirB4", "VirB6", "VirB5"))
bad = 0
for name, fn in CONTROLS:
    p = os.path.join(GB, fn)
    if not os.path.exists(p):
        print("  %-7s RECORD MISSING" % name); bad += 1; continue
    rec = next(SeqIO.parse(p, "genbank"))
    calls = orf_caller.call(str(rec.seq))
    digs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=c["aa"]).digitize(alpha)
            for i, c in enumerate(calls)]
    got = {}
    for lab, hmm, mode, minaa in models:
        kw = {"bit_cutoffs": "gathering"} if mode == "ga" else {"E": 1e-5}
        for top in pyhmmer.hmmsearch([hmm], digs, cpus=1, **kw):
            for h in top:
                nm = h.name
                i = int(nm.decode() if isinstance(nm, bytes) else nm)
                if minaa and calls[i]["aa_len"] <= minaa:
                    continue
                if lab not in got or h.score > got[lab][1]:
                    got[lab] = (i, h.score)
    ok = all(k in got for k in ("VirB4", "VirB6", "VirB5"))
    if not ok:
        bad += 1
    print("  %-7s %-8s %-8s %-8s  %s"
          % (name,
             "%.0f" % got["VirB4"][1] if "VirB4" in got else "MISS",
             "%.0f" % got["VirB6"][1] if "VirB6" in got else "MISS",
             "%.0f" % got["VirB5"][1] if "VirB5" in got else "MISS",
             "OK" if ok else "FAIL"))
# ---- negative assertions
print("\n  negative controls -- must be VirB4-NEGATIVE")
print("  %-13s %-8s  verdict" % ("control", "VirB4"))
for name, fn, note in NEGATIVES:
    p = os.path.join(GB, fn)
    if not os.path.exists(p):
        print("  %-13s RECORD MISSING" % name); continue
    rec = next(SeqIO.parse(p, "genbank"))
    calls = orf_caller.call(str(rec.seq))
    digs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=c["aa"]).digitize(alpha)
            for i, c in enumerate(calls)]
    best = None
    for lab, hmm, mode, minaa in models:
        if lab != "VirB4":
            continue
        for top in pyhmmer.hmmsearch([hmm], digs, cpus=1, bit_cutoffs="gathering"):
            for h in top:
                if best is None or h.score > best:
                    best = h.score
    if name == "pKM101_frag":
        # NOT a true negative: this is the pKM101 tra fragment and it DOES carry
        # VirB4. It is here to prove the negative test can distinguish -- if this
        # one came out negative too, the VirB4 scan would be dead rather than
        # discriminating.
        ok = best is not None
        print("  %-13s %-8s  %s  (positive by design: proves the negative test discriminates)"
              % (name, "%.0f" % best if best else "MISS", "OK" if ok else "FAIL"))
    else:
        ok = best is None
        print("  %-13s %-8s  %s  %s"
              % (name, "%.0f" % best if best else "none", "OK" if ok else "FAIL -- false positive",
                 note or ""))
    if not ok:
        bad += 1

print("\n%s" % ("PASS -- step 1a may be read" if bad == 0
                else "FAIL: %d control(s) broken; step 1a output is INVALID" % bad))
sys.exit(1 if bad else 0)
