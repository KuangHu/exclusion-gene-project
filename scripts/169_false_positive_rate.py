#!/usr/bin/env python3
"""False-positive rate of the ANCHOR_CANDIDATE criterion.

Never measured. Without it, no pass rate anywhere in this module is
interpretable -- "45.8% of the slot is lipid-anchored" means nothing until we
know what fraction of ARBITRARY proteins the criterion admits.

The criterion under test (dropping the positional term, which is what makes it a
background test):

    length <= 250 aa  AND  (lipobox_strict OR tm_count >= 1)

FOUR SETS, deliberately spanning the expected range:

  1. RANDOM       proteins from the same slot-bearing plasmids, <=250 aa,
                  excluding anything at a slot and anything hit by an anchor
                  model. This is the true background: what the criterion admits
                  from arbitrary small plasmid proteins.
  2. CYTOPLASMIC  VirB11 (T2SSE) and VirB4 ATPases -- known soluble, known
                  cytoplasmic. MUST be near zero. If they pass, the criterion is
                  measuring nothing.
  3. MEMBRANE     VirB6/TrbL, polytopic inner membrane. MUST be high. This is the
                  other-direction sanity check: a criterion that rejects known
                  membrane proteins is broken in the opposite way.
  4. SLOT         the actual slot occupants, for comparison.

The ATPases are length-filtered to <=250 aa like everything else, so they are
size-matched rather than rejected on length -- otherwise set 2 would pass
trivially and prove nothing.

WHAT A GOOD RESULT LOOKS LIKE: set 2 near zero, set 3 high, and set 1 materially
below set 4. If set 1 is close to set 4, the criterion is admitting background at
the same rate as slots and its output is not a candidate list.

There is no threshold to tune here -- the criterion is already frozen. This
measures it.
"""
import collections
import csv
import os
import random
import re
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import decontaminate as DC

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
WORK = "/global/home/users/kh36969/tmp/fpr"
STRICT = re.compile(r"[LVI][ASTVIG][GASN]C")
MAXLEN = 250
N = 600


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


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    os.makedirs(WORK, exist_ok=True)
    rng = random.Random(0)

    sof = [r for r in csv.DictReader(open(SOF), delimiter="\t")
           if r["slot_occupant_coords"]]
    want = collections.defaultdict(set)
    for r in sof:
        want[r["accession"]].add(r["slot_occupant_coords"])

    seqs, aas, owner, is_slot = [], [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                if acc not in want:
                    continue
                p = line.rstrip("\n")
                k = "%s..%s" % (f[2], f[3])
                aas.append(p); owner.append(acc); is_slot.append(k in want[acc])
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(aas) - 1).encode(), sequence=p).digitize(alpha))
    print("proteins on slot-bearing plasmids: %d (%d at a slot)"
          % (len(aas), sum(is_slot)), flush=True)

    def scan(fam):
        p = os.path.join(HMM, fam + ".hmm")
        if not os.path.exists(p):
            return set()
        with pyhmmer.plan7.HMMFile(p) as fh:
            m = next(iter(fh))
        got = set()
        for top in pyhmmer.hmmsearch([m], seqs, cpus=8, bit_cutoffs="gathering"):
            for h in top:
                got.add(int(dec(h.name)))
        return got

    anchor_any = set()
    for _, fam, _, _ in DC.ANCHOR_FAMILIES:
        anchor_any |= scan(fam)
    cyto = scan("anchor_T2SSE") | scan("anchor_CagE_TrbE_VirB")   # VirB11, VirB4
    memb = scan("TrbL")                                            # VirB6
    print("  anchor-family proteins %d | VirB11+VirB4 %d | VirB6 %d"
          % (len(anchor_any), len(cyto), len(memb)), flush=True)

    pools = {}
    pools["4_SLOT"] = [i for i in range(len(aas)) if is_slot[i] and len(aas[i]) <= MAXLEN]
    pools["1_RANDOM"] = [i for i in range(len(aas))
                         if (not is_slot[i]) and i not in anchor_any
                         and len(aas[i]) <= MAXLEN]
    pools["2_CYTOPLASMIC"] = [i for i in cyto if len(aas[i]) <= MAXLEN]
    pools["3_MEMBRANE"] = [i for i in memb if len(aas[i]) <= MAXLEN]

    fa = os.path.join(WORK, "in.faa")
    idx = []
    with open(fa, "w") as fh:
        for lab in sorted(pools):
            p = pools[lab]
            s = p if len(p) <= N else rng.sample(p, N)
            print("  %-16s pool %6d  tested %4d" % (lab, len(p), len(s)), flush=True)
            for j, i in enumerate(s):
                nm = "%s|%d" % (lab, j)
                idx.append((nm, lab, aas[i]))
                fh.write(">%s\n%s\n" % (nm, aas[i]))
    if not idx:
        print("nothing to test"); return 1

    pred = os.path.join(WORK, "out.pred")
    r = subprocess.run(["tmbed", "predict", "-f", fa, "-p", pred, "--out-format", "0",
                        "--model-dir", "/global/home/users/kh36969/.hfcache/prott5",
                        "--use-gpu"], capture_output=True, text=True)
    if r.returncode:
        print("tmbed failed:\n%s" % r.stderr[-600:]); return 1
    L = [l.rstrip("\n") for l in open(pred)]
    got = {}
    for i in range(0, len(L), 3):
        if i + 2 < len(L) and L[i].startswith(">"):
            got[L[i][1:]] = (L[i + 1], L[i + 2])

    agg = collections.defaultdict(lambda: collections.Counter())
    rows = []
    for nm, lab, seq in idx:
        if nm not in got:
            continue
        s2, p = got[nm]
        tm = len([x for x in segments(p) if x[0] == "TM"])
        lb = bool(STRICT.search(s2[:40]))
        passes = lb or tm >= 1
        mode = ("lipid" if (lb and tm == 0) else "single_pass" if (not lb and tm == 1)
                else "multi_pass" if (not lb and tm >= 2) else "ambiguous" if lb else "none")
        a = agg[lab]
        a["n"] += 1; a["pass"] += passes; a["lipobox"] += lb
        a["tm1"] += (tm == 1); a["tm2plus"] += (tm >= 2)
        a[mode] += 1
        rows.append({"set": lab, "id": nm, "aa": len(s2), "lipobox": int(lb),
                     "tm_count": tm, "anchor_mode": mode, "passes": int(passes)})

    print("\n=== FALSE-POSITIVE RATE of ANCHOR_CANDIDATE ===")
    print("  criterion: length <= %d AND (lipobox_strict OR tm_count >= 1)\n" % MAXLEN)
    print("  %-16s %6s %9s %9s %8s %9s" % ("set", "n", "PASSES", "lipobox", "TM=1", "TM>=2"))
    for lab in sorted(agg):
        a = agg[lab]; n = a["n"]
        print("  %-16s %6d %8.1f%% %8.1f%% %7.1f%% %8.1f%%"
              % (lab, n, 100.0 * a["pass"] / n, 100.0 * a["lipobox"] / n,
                 100.0 * a["tm1"] / n, 100.0 * a["tm2plus"] / n))

    print("\n=== VERDICT ===")
    g = {lab: 100.0 * agg[lab]["pass"] / agg[lab]["n"] for lab in agg if agg[lab]["n"]}
    cy, me = g.get("2_CYTOPLASMIC"), g.get("3_MEMBRANE")
    bg, sl = g.get("1_RANDOM"), g.get("4_SLOT")
    ok = True
    if cy is not None:
        print("  cytoplasmic control %.1f%%  %s" % (cy, "OK" if cy <= 20 else "TOO HIGH -- criterion admits soluble proteins"))
        ok &= cy <= 20
    if me is not None:
        print("  membrane control    %.1f%%  %s" % (me, "OK" if me >= 80 else "TOO LOW -- criterion rejects known membrane proteins"))
        ok &= me >= 80
    if bg is not None and sl is not None:
        print("  background %.1f%% vs slot %.1f%%  -> enrichment %.2fx" % (bg, sl, sl / bg if bg else float("inf")))
        print("\n  %s" % ("Criterion discriminates; background rate is %.1f%% and must be"
                          " quoted alongside any slot figure." % bg if sl > bg * 1.3 else
                          "Criterion admits background at ~the slot rate. Its output is"
                          " NOT a candidate list -- it is a size-and-membrane filter."))
    dest = os.path.join(PROJ, "data", "anchors", "false_positive_rate.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))


if __name__ == "__main__":
    sys.exit(main() or 0)
