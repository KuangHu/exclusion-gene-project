#!/usr/bin/env python3
"""Is the MPF_I slot TraY+1 outside the architecture the seeds happen to occupy?

The slot rule rests on 5/5 seeds carrying a named excA at TraY+1 with an identical
local order -- traW(400) traX(194) traY(722) excA(220) 69aa, protein lengths and
all. The seed representativeness check then showed all five sit in ONE canonical
order (#1) and ONE layout (contiguous), sampling 1 of 118 orders and 45.7% of the
admitted set. So that validation says nothing about the other 54.3%, and the
uniformity was a property of the seed set, not of MPF_I.

Building the catalogue now would apply a rule verified on 45.7% to the whole class.

Stratified by architecture, never pooled. MPF_T is precedent for why: its three
architectures were found exactly this way -- IncI2's bimodal +1/+2 and pEC4115's
-1 would both have vanished into an average. When an offset disagrees, the first
question is whether it is a CONSISTENT other offset, which is a new architecture,
rather than noise.
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    csv.field_size_limit(10 ** 8)

    surv = {r["accession"]: r for r in csv.DictReader(
        open(os.path.join(PROJ, "data", "anchors", "mpfi_architecture_survey.tsv")),
        delimiter="\t")}
    oc = collections.Counter(r["order_canonical"] for r in surv.values())
    rank = {o: i + 1 for i, (o, _) in enumerate(oc.most_common())}
    print("architecture survey: %d plasmids, %d canonical orders" % (len(surv), len(oc)))

    seqs, owner, idx, strand, aas = [], [], [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                owner.append(f[0]); idx.append(int(f[1])); strand.append(int(f[4]))
                aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    print("cache: %d proteins" % len(seqs), flush=True)

    def scan(path):
        with pyhmmer.plan7.HMMFile(path) as fh:
            m = next(iter(fh))
        best = collections.defaultdict(list)
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                best[owner[i]].append((idx[i], strand[i], h.score, i))
        return best
    tray = scan(os.path.join(CJ, "T4SS_I_traY.hmm"))
    exca = scan(os.path.join(HMM, "NF033891.hmm"))
    print("TraY on %d plasmids; ExcA on %d" % (len(tray), len(exca)), flush=True)
    best = lambda d, acc: max(d[acc], key=lambda x: x[2]) if acc in d else None

    def stratum(acc):
        r = surv.get(acc)
        if not r:
            return None
        k = rank.get(r["order_canonical"], 999)
        lay = r["layout"]
        if lay != "contiguous":
            return lay
        return "order#%d" % k if k <= 3 else "order#4+"

    rows = []
    agg = collections.defaultdict(lambda: collections.Counter())
    offs = collections.defaultdict(collections.Counter)
    for acc in surv:
        s = stratum(acc)
        if not s:
            continue
        agg[s]["n"] += 1
        ty, ex = best(tray, acc), best(exca, acc)
        if not ty:
            agg[s]["no_traY"] += 1
        if not ex:
            agg[s]["no_excA"] += 1
        if not (ty and ex):
            continue
        agg[s]["both"] += 1
        t = 1 if ty[1] == 1 else -1
        d = (ex[0] - ty[0]) * t
        offs[s][max(-9, min(9, d))] += 1
        if d == 1:
            agg[s]["at_plus1"] += 1
        rows.append({"accession": acc, "stratum": s, "offset": d,
                     "traY_bit": round(ty[2], 1), "excA_bit": round(ex[2], 1),
                     "excA_aa": len(aas[ex[3]])})

    print("\n=== ExcA position relative to TraY, STRATIFIED (never pooled) ===")
    print("%-12s %7s %9s %9s %9s %11s" % ("stratum", "n", "TraY-", "ExcA-", "both",
                                          "excA at +1"))
    order = ["order#1", "order#2", "order#3", "order#4+", "split_2", "split_3",
             "split_4", "split_5", "split_6"]
    for s in order:
        if s not in agg:
            continue
        c = agg[s]
        print("%-12s %7d %9d %9d %9d %10s"
              % (s, c["n"], c["n"] - c.get("no_traY", 0), c["n"] - c.get("no_excA", 0),
                 c["both"], "%d (%.1f%%)" % (c["at_plus1"],
                 100.0 * c["at_plus1"] / max(1, c["both"]))))
    print("\n=== offset distribution per stratum ===")
    for s in order:
        if s not in offs:
            continue
        tot = sum(offs[s].values())
        det = ", ".join("%+d:%d(%.0f%%)" % (k, v, 100.0 * v / tot)
                        for k, v in sorted(offs[s].items(), key=lambda kv: -kv[1])[:5])
        print("  %-12s n=%-6d %s" % (s, tot, det))
    print("\n  A CONSISTENT other offset is a new architecture, not noise -- that is")
    print("  how MPF_T's IncI2 (+2) and pEC4115 (-1) were found.")

    dest = os.path.join(PROJ, "data", "anchors", "mpfi_slot_by_architecture.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))


if __name__ == "__main__":
    main()
