#!/usr/bin/env python3
"""D.1 -- do the 2014 specificity statistics separate EEx from machine components?

THIS RUNS BEFORE ANY FAMILY IS BUILT. Guglielmini 2014 built its profile set by
taking +/-20 genes around every VirB4, clustering at >=30% identity / >50% overlap,
and keeping families that are SPECIFIC to one T4SS class -- high f_own, low
f_other. That criterion is designed for components belonging to one machine.

A successful entry-exclusion family is the opposite shape: it spreads across
plasmid types, because exclusion specificity is what is under selection, not class
membership. This project already recorded the same trap once -- the `in-slot/total`
ratio selected for RARITY rather than function and rejected a known positive.

So before screening anything, put the four KNOWN EEx through the same statistics
as known machine components and see whether the two occupy separable regions. If
they do not, this workflow cannot find EEx no matter how well it repairs the
annotation gap, and that has to be known first, not discovered afterwards.

STATISTICS, per family, over +/-20 genes around VirB4 (2014's window):

  f_own        hits inside the OWN class's virB4 window        machine: high
  f_other      hits inside ANOTHER class's virB4 window        machine: low
  f_none       hits not near any VirB4                         machine: low
  occupancy    own-class elements carrying it                  machine: high
  pos_entropy  entropy of the offset-from-virB4 distribution   machine: LOW
  mean_pid     mean pairwise identity among hits
  div_ratio    mean_pid(family) / mean_pid(virB4), same elements

pos_entropy is not in the 2014 paper. It asks whether a family sits at a FIXED
position in the window rather than merely somewhere inside it, which is the
stronger signature of an operon component. div_ratio normalises divergence against
VirB4 measured on the same elements: ~1 means co-evolving with the machine, <<1
means evolving independently or under diversifying selection.

NONE OF THESE IS A GATE HERE. All are recorded. The point is the geometry: where
the EEx land relative to the machine components.

REFERENCE SETS
  EEx          TrbK(T) TraS_R100(F) ExcA(I) DUF4467(FA) Eex_IncN(T) EexR/S(F)
  machine      per-class CONJScan profiles, which by construction are what 2014
               kept -- they define the "machine" region of the space
  missed       traN, traK, traG -- known machine genes CONJScan misses in some
               elements; they should sit WITH the machine, and are the positive
               control for step 7 later
"""
import collections
import csv
import math
import os
import statistics
import sys
import xml.etree.ElementTree as ET

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CJ = os.path.join(PROJ, "data", "CONJScan", "profiles")
DEFS = os.path.join(PROJ, "data", "CONJScan", "definitions", "Plasmids")
HMM = os.path.join(PROJ, "data", "hmm")
REL = os.path.join(PROJ, "data", "release")
ANCH = os.path.join(PROJ, "data", "anchors")
WINDOW = 20                    # 2014's +/-20 genes around virB4
VIRB4 = ["T4SS_virb4", "T4SS_I_traU"]
CLASSES = ["T", "F", "I", "FA"]

# family -> (model file, model dir, own class, label)
EEX = [("TIGR04359", HMM, "T", "TrbK"), ("PF10624", HMM, "F", "TraS_R100"),
       ("NF033891", HMM, "I", "ExcA"), ("PF14729", HMM, "FA", "DUF4467"),
       ("NF033894", HMM, "T", "Eex_IncN"), ("NF041429", HMM, "F", "EexR/S")]
# machine genes CONJScan misses in some elements -- the step-7 positive control
MISSED = [("T4SS_F_traN", CJ, "F", "traN"), ("T4SS_F_traK", CJ, "F", "traK"),
          ("T4SS_F_traG", CJ, "F", "traG")]


def class_profiles(cls):
    names = {g.get("name") for g in
             ET.parse(os.path.join(DEFS, "T4SS_type%s.xml" % cls)).getroot().iter("gene")}
    return sorted(n for n in names if not n.startswith("T4SS_MOB"))


def pid(a, b):
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
    import random
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
    cls_of = {}
    for c in CLASSES:
        for e in members[c]:
            cls_of.setdefault(e, c)          # first wins; overlaps are rare
    want = set(cls_of)
    print("elements: %s  total %d" % ({c: len(members[c]) for c in CLASSES}, len(want)))

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
    print("proteins: %d\n" % len(prot), flush=True)
    seqs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=p).digitize(alpha)
            for i, p in enumerate(prot)]
    ngene = collections.Counter(owner)

    def scan(path, E=1e-3):
        if not os.path.exists(path):
            return {}
        with pyhmmer.plan7.HMMFile(path) as fh:
            hmm = next(iter(fh))
        out = {}
        for top in pyhmmer.hmmsearch([hmm], seqs, cpus=16, E=E):
            for h in top:
                i = int(dec(h.name))
                if i not in out or h.score > out[i]:
                    out[i] = h.score
        return out

    # ---- virB4 position per element ----------------------------------------
    v4 = {}
    for m in VIRB4:
        for i, sc in scan(os.path.join(CJ, m + ".hmm")).items():
            e = owner[i]
            if e not in v4 or sc > v4[e][1]:
                v4[e] = (gidx[i], sc, i)
    print("elements with a located VirB4: %d of %d (%.1f%%)\n"
          % (len(v4), len(want), 100.0 * len(v4) / len(want)), flush=True)

    # VirB4's own divergence, per class, as the div_ratio denominator
    v4pid = {}
    for c in CLASSES:
        idx = [v4[e][2] for e in members[c] if e in v4]
        s = idx if len(idx) <= 40 else rng.sample(idx, 40)
        vals = [pid(prot[x], prot[y]) for a, x in enumerate(s) for y in s[a + 1:]]
        v4pid[c] = statistics.mean(vals) if vals else float("nan")
    print("VirB4 mean pairwise identity: %s\n"
          % {c: round(v4pid[c], 1) for c in CLASSES}, flush=True)

    def stats(hits, own):
        """hits: protein index -> score."""
        if not hits:
            return None
        offs, nown, nother, nnone = [], 0, 0, 0
        els = set()
        for i in hits:
            e = owner[i]
            els.add(e)
            if e not in v4:
                nnone += 1
                continue
            d = gidx[i] - v4[e][0]
            n = ngene[e]
            if d > n / 2:
                d -= n                      # circular
            elif d < -n / 2:
                d += n
            if abs(d) <= WINDOW:
                if cls_of.get(e) == own:
                    nown += 1
                else:
                    nother += 1
                offs.append(d)
            else:
                nnone += 1
        tot = len(hits)
        h = 0.0
        if offs:
            cnt = collections.Counter(offs)
            for v in cnt.values():
                p = v / len(offs)
                h -= p * math.log2(p)
        idx = list(hits)
        s = idx if len(idx) <= 40 else rng.sample(idx, 40)
        vals = [pid(prot[x], prot[y]) for a, x in enumerate(s) for y in s[a + 1:]]
        mp = statistics.mean(vals) if vals else float("nan")
        occ = 100.0 * len(els & members[own]) / len(members[own])
        return {"n_hits": tot, "n_elem": len(els),
                "f_own": 100.0 * nown / tot, "f_other": 100.0 * nother / tot,
                "f_none": 100.0 * nnone / tot, "occupancy": occ,
                "pos_entropy": h, "mean_pid": mp,
                "div_ratio": mp / v4pid[own] if v4pid.get(own) else float("nan"),
                "med_offset": statistics.median(offs) if offs else float("nan")}

    rows = []
    print("=" * 96)
    print("D.1 CALIBRATION -- where do known EEx sit relative to known machine parts?")
    print("=" * 96)
    print("  %-14s %-5s %-9s %7s %7s %7s %7s %7s %7s %7s %7s"
          % ("family", "class", "kind", "n_hits", "f_own", "f_oth", "f_none",
             "occ%", "posH", "mean_id", "div_r"))

    def emit(label, own, kind, path):
        st = stats(scan(path), own)
        if not st:
            print("  %-14s %-5s %-9s  no hits" % (label, own, kind)); return
        print("  %-14s %-5s %-9s %7d %6.1f%% %6.1f%% %6.1f%% %6.1f%% %7.2f %7.1f %7.2f"
              % (label, own, kind, st["n_hits"], st["f_own"], st["f_other"],
                 st["f_none"], st["occupancy"], st["pos_entropy"], st["mean_pid"],
                 st["div_ratio"]), flush=True)
        rows.append(dict(family=label, own_class="MPF_" + own, kind=kind, **st))

    for mdl, d, own, lab in EEX:
        emit(lab, own, "EEx", os.path.join(d, mdl + ".hmm"))
    print()
    for mdl, d, own, lab in MISSED:
        emit(lab, own, "missed", os.path.join(d, mdl + ".hmm"))
    print()
    for c in CLASSES:
        for p in class_profiles(c):
            if p in VIRB4:
                continue
            emit(p.replace("T4SS_", ""), c, "machine", os.path.join(CJ, p + ".hmm"))

    # ---- the separation question -------------------------------------------
    print("\n" + "=" * 96)
    print("SEPARATION")
    print("=" * 96)
    for k in ("f_own", "f_other", "f_none", "pos_entropy", "occupancy", "div_ratio"):
        e = [r[k] for r in rows if r["kind"] == "EEx" and not math.isnan(r[k])]
        m = [r[k] for r in rows if r["kind"] == "machine" and not math.isnan(r[k])]
        if not e or not m:
            continue
        print("  %-12s EEx median %8.2f (range %6.2f-%7.2f, n=%d) | machine median"
              " %8.2f (range %6.2f-%7.2f, n=%d)"
              % (k, statistics.median(e), min(e), max(e), len(e),
                 statistics.median(m), min(m), max(m), len(m)))
    print("\n  A statistic separates only if the EEx range and the machine range")
    print("  barely overlap. Any that does not is a field to record, never a filter.")

    dest = os.path.join(ANCH, "D1_calibration.tsv")
    if rows:
        with open(dest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(rows)
        print("\nwrote %s (%d rows)" % (dest, len(rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
