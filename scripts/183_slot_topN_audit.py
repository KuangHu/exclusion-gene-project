#!/usr/bin/env python3
"""Audit the top slots: are they EEx-shaped, or just hypervariable machine parts?

The gate passed 4/4 with all known EEx slots in the top 4. That is only
meaningful if the OTHER top slots are not simply pilins. VirB2 and VirB5 are
short, hypervariable and membrane-embedded -- the same layer-2 signature as an
EEx -- and this project already recorded that Pfam GA fails on them
systematically. A pilin below GA is INVISIBLE to the recogniser and gets passed
through as an unannotated slot occupant. If the non-EEx members of the top 4 are
pilin slots, layer 2 is measuring HYPERVARIABILITY, not exclusion, and
conjugative elements contain plenty of hypervariable things.

THREE CHECKS, all on the calibration classes:

  1. WHAT ELSE IS IN THE TOP. Occupants of the top slots are re-scanned against
     every CONJScan profile at a RELAXED E<=1e-3, not GA. Anything that lights up
     sub-GA is a machine component the GA-based recogniser missed. Pfam is not
     used for this -- the class-specific CONJScan profiles are what model pilins.

  2. DENOMINATOR SENSITIVITY. The gate's 77 / 61 / 60 / 19 come from
     MIN_DEFINED >= 30. If relaxing it to 10 triples the denominator and the
     known EEx fall out of the top 4, then the rank was produced by the
     pre-filter, not by the ranking. Ranks are recomputed at 10 / 30 / 60.

     MPF_FA's denominator of 19 is called out separately: rank 2 of 19 is about
     the 10th percentile, which on its own is close to uninformative, and the
     4/4 should not be read as four equally strong results.

  3. HOW BLIND WAS BLIND. Recorded, not measured: the RANKING never saw the EEx
     labels, but the SLOT CONSTRUCTION -- which profiles count as machine, the
     3-gene block-break rule, MIN_DEFINED -- was written by someone who already
     knew where the four EEx sit. The ranking is blind; the scaffold is not. On
     G/B/C there is no known EEx to tune against, so any implicit tuning here
     will not transfer, and a sharp drop in G/B/C performance should be blamed
     here first.
"""
import collections
import csv
import gzip
import math
import os
import random
import statistics
import sys
import xml.etree.ElementTree as ET

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import decontaminate as DC

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CJ = os.path.join(PROJ, "data", "CONJScan", "profiles")
DEFS = os.path.join(PROJ, "data", "CONJScan", "definitions", "Plasmids")
HMM = os.path.join(PROJ, "data", "hmm")
REL = os.path.join(PROJ, "data", "release")
ANCH = os.path.join(PROJ, "data", "anchors")
MAXGAP = 3
MAX_PAIRS = 30
TOPN = 8
RELAXED_E = 1e-3
CLASSES = ["T", "F", "I", "FA"]
EEX = [("TIGR04359", "TrbK"), ("PF10624", "TraS_R100"), ("NF033891", "ExcA"),
       ("NF033894", "Eex_IncN"), ("NF041429", "EexR/S"), ("PF14729", "DUF4467")]
GENERIC = ["T4SS_virb4", "T4SS_I_traU", "T4SS_t4cp1", "T4SS_t4cp2", "T4SS_tcpA"]


def class_profiles(cls):
    names = {g.get("name") for g in
             ET.parse(os.path.join(DEFS, "T4SS_type%s.xml" % cls)).getroot().iter("gene")}
    return sorted(n for n in names if not n.startswith("T4SS_MOB"))


def ident(a, b):
    best = 0.0
    for off in range(-6, 7, 3):
        x, y = (a, b[off:]) if off >= 0 else (a[-off:], b)
        n = min(len(x), len(y))
        if n < 20:
            continue
        best = max(best, sum(1 for p, q in zip(x, y) if p == q) / n)
    return 100.0 * best


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    rng = random.Random(0)

    members = {}
    for c, p in (("T", "v1.1/catalogue_MPF_T_v1.1.tsv"),
                 ("F", "v1/catalogue_MPF_F_v1.tsv"),
                 ("I", "v1.1/catalogue_MPF_I_v1.tsv")):
        members[c] = {l.split("\t")[0] for i, l in
                      enumerate(open(os.path.join(REL, p))) if i}
    fa = list(csv.DictReader(open(os.path.join(REL, "v1.1/catalogue_MPF_FA_v1.tsv")),
                             delimiter="\t"))
    members["FA"] = {r["accession"] for r in fa if r["fa_subclass"] == "FA"}
    want = set().union(*members.values())

    prot, owner, gidx = [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                if f[0] in want:
                    owner.append(f[0]); gidx.append(int(f[1]))
                    prot.append(line.rstrip("\n"))
    print("proteins: %d over %d plasmids\n" % (len(prot), len(want)), flush=True)
    seqs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=p).digitize(alpha)
            for i, p in enumerate(prot)]

    def scan(path, E=None):
        if not os.path.exists(path):
            return {}
        with pyhmmer.plan7.HMMFile(path) as fh:
            hmm = next(iter(fh))
        out = {}
        kw = {"E": E} if E else {"bit_cutoffs": "gathering"}
        try:
            it = pyhmmer.hmmsearch([hmm], seqs, cpus=16, **kw)
        except Exception:
            it = pyhmmer.hmmsearch([hmm], seqs, cpus=16, E=1e-5)
        for top in it:
            for h in top:
                i = int(dec(h.name))
                if i not in out or h.score > out[i]:
                    out[i] = h.score
        return out

    allprof = sorted({p for c in CLASSES for p in class_profiles(c)} | set(GENERIC))
    prof_ga, prof_rx = {}, collections.defaultdict(list)
    for j, p in enumerate(allprof):
        for i in scan(os.path.join(CJ, p + ".hmm")):
            prof_ga.setdefault(i, p)
        for i, sc in scan(os.path.join(CJ, p + ".hmm"), E=RELAXED_E).items():
            prof_rx[i].append((sc, p))
        if (j + 1) % 20 == 0:
            print("  profiles %d/%d" % (j + 1, len(allprof)), flush=True)
    for _, famf, _, _ in DC.ANCHOR_FAMILIES:
        for i in scan(os.path.join(HMM, famf + ".hmm")):
            prof_ga.setdefault(i, "anchor_" + famf)
    eex_at = {}
    for mdl, nm in EEX:
        for i in scan(os.path.join(HMM, mdl + ".hmm")):
            eex_at[i] = nm
    print("\nGA-recognised %d | sub-GA relaxed hits %d | known EEx %d\n"
          % (len(prof_ga), len(prof_rx), len(eex_at)), flush=True)

    byel = collections.defaultdict(list)
    for i in range(len(prot)):
        byel[owner[i]].append(i)
    for e in byel:
        byel[e].sort(key=lambda i: gidx[i])

    out_rows = []
    for cls in CLASSES:
        els = members[cls]
        cprof = set(class_profiles(cls)) | set(GENERIC)
        slot_occ = collections.defaultdict(list)
        slot_def = collections.Counter()
        for e in els:
            idxs = byel.get(e, [])
            mach = sorted((gidx[i], prof_ga[i]) for i in idxs
                          if i in prof_ga and prof_ga[i] in cprof)
            if len(mach) < 3:
                continue
            pos = {gidx[i]: i for i in idxs}
            for k in range(len(mach) - 1):
                (g1, p1), (g2, p2) = mach[k], mach[k + 1]
                if g2 - g1 - 1 > MAXGAP:
                    continue
                slot_def[(p1, p2)] += 1
                for g in range(g1 + 1, g2):
                    if g in pos and pos[g] not in prof_ga:
                        slot_occ[(p1, p2)].append(pos[g])

        def score(min_def):
            rows = []
            for key, occ in slot_occ.items():
                if slot_def[key] < min_def:
                    continue
                lens = [len(prot[i]) for i in occ]
                smp = occ if len(occ) <= MAX_PAIRS else rng.sample(occ, MAX_PAIRS)
                ids = [ident(prot[x], prot[y])
                       for ii, x in enumerate(smp) for y in smp[ii + 1:]]
                rows.append({"key": key, "occ": occ,
                             "occupancy": 100.0 * len({owner[i] for i in occ}) / slot_def[key],
                             "mpid": statistics.mean(ids) if ids else 100.0,
                             "cv": (statistics.pstdev(lens) / statistics.mean(lens))
                                   if len(lens) > 1 else 0.0,
                             "med": statistics.median(lens),
                             "eex": collections.Counter(eex_at[i] for i in occ
                                                        if i in eex_at)})
            if not rows:
                return []
            dims = ["occupancy", "mpid", "cv"]
            mu = {d: statistics.mean(r[d] for r in rows) for d in dims}
            sd = {d: (statistics.pstdev([r[d] for r in rows]) or 1.0) for d in dims}
            ref = [r for r in rows if r["eex"]]
            if not ref:
                return []
            cen = {d: statistics.mean((r[d] - mu[d]) / sd[d] for r in ref) for d in dims}
            for r in rows:
                z = {d: (r[d] - mu[d]) / sd[d] for d in dims}
                r["dist"] = math.sqrt(sum((z[d] - cen[d]) ** 2 for d in dims))
            rows.sort(key=lambda r: r["dist"])
            return rows

        print("=" * 78)
        print("MPF_%s -- DENOMINATOR SENSITIVITY" % cls)
        print("  %-10s %8s %s" % ("min_defined", "n_slots", "rank of each known EEx slot"))
        for md in (10, 30, 60):
            rr = score(md)
            if not rr:
                print("  %-10d %8s no slot carries a known EEx" % (md, "-")); continue
            kn = [(j + 1, r["eex"].most_common(1)[0][0])
                  for j, r in enumerate(rr) if r["eex"]]
            print("  %-10d %8d %s" % (md, len(rr),
                  ", ".join("%s rank %d" % (n, j) for j, n in kn[:4])))

        rows = score(30)
        if not rows:
            print(); continue
        print("\n  TOP %d SLOTS -- what is actually in them" % TOPN)
        print("  %-4s %-28s %6s %6s %6s %5s  %s"
              % ("rank", "slot", "occ%", "mpid", "med_aa", "nEEx", "sub-GA machine hits"))
        for j, r in enumerate(rows[:TOPN]):
            sub = collections.Counter()
            for i in r["occ"]:
                if i in prof_rx:
                    sub[max(prof_rx[i])[1]] += 1
            tot = len(r["occ"])
            top = ", ".join("%s %.0f%%" % (p.replace("T4SS_", ""), 100.0 * n / tot)
                            for p, n in sub.most_common(3)) or "none"
            print("  %-4d %-28s %5.1f%% %6.1f %6d %5d  %s"
                  % (j + 1, "%s>%s" % r["key"], r["occupancy"], r["mpid"], r["med"],
                     sum(r["eex"].values()), top))
            out_rows.append({"mpf_class": "MPF_" + cls, "rank": j + 1,
                             "slot": "%s>%s" % r["key"], "occupancy_pct": round(r["occupancy"], 1),
                             "mean_pid": round(r["mpid"], 1), "len_median": r["med"],
                             "n_occupants": tot, "n_known_eex": sum(r["eex"].values()),
                             "known_eex": ";".join(r["eex"]),
                             "subGA_machine_pct": round(100.0 * sum(sub.values()) / tot, 1),
                             "subGA_top": top})
        print()

    print("=" * 78)
    print("READING")
    print("  A top slot whose occupants light up a class-specific CONJScan profile")
    print("  below GA is a MACHINE COMPONENT the recogniser missed -- most likely a")
    print("  pilin. If the non-EEx members of the top 4 are those, layer 2 is")
    print("  scoring hypervariability rather than exclusion.")
    dest = os.path.join(ANCH, "slot_topN_audit.tsv")
    if out_rows:
        with open(dest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(out_rows)
        print("\nwrote %s (%d rows)" % (dest, len(out_rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
