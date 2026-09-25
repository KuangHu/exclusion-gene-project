#!/usr/bin/env python3
"""The anchor ladder: per-family detection failure rate over the VirB4+ set.

VirB5 fails on 21.4% of VirB4+ plasmids. That number only means something against
the rate for anchors of the same kind, and until now only four of the ten VirB
families had been measured (VirB4/5/6/8/9). VirB1, VirB2, VirB3, VirB10 and
VirB11 were never scanned, so the "~2% baseline" rested on three points.

This completes the ladder. Every family in the CURATED anchor set is scanned at
its own declared threshold -- GA for most, `E<=1e-5 AND aa>150` for VirB5, whose
scores are bimodal on IncP and which GA misses outright.

Interpretation, stated before the numbers exist so it cannot be fitted to them:

  * VirB5 near the other families    -> 21.4% is ordinary, the slot-empty test is
                                        no weaker than any other anchor
  * VirB5 an outlier by several fold -> family-coverage bias specific to PF07996,
                                        and every "empty slot" inherits it

A high failure rate on the SLOT-BOUNDING anchors (VirB5, VirB6) is not symmetric
with a high rate on annotation-only anchors: a missed VirB5 turns an occupied slot
into no slot at all, which is the unfilterable direction.

Red line (v2 1.3): the exclusion families never appear here. This scans anchors.
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))

HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
FORBIDDEN = {"TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"}


def read_anchor_set():
    rows = []
    p = os.path.join(PROJ, "data", "anchors", "anchor_set.tsv")
    with open(p) as fh:
        lines = [l for l in fh if not l.startswith("#") and l.strip()]
    r = csv.DictReader(lines, delimiter="\t")
    for d in r:
        if d["pfam_acc"] in FORBIDDEN or d["pfam_name"] in FORBIDDEN:
            raise SystemExit("RED LINE: exclusion family %s in the anchor set"
                             % d["pfam_acc"])
        rows.append(d)
    return rows


def hmm_path(name):
    for cand in ("anchor_%s.hmm" % name, "%s.hmm" % name):
        p = os.path.join(HMM, cand)
        if os.path.exists(p):
            return p
    raise SystemExit("no HMM file for %r in %s" % (name, HMM))


def parse_threshold(t):
    """-> (evalue_or_None, min_aa). None evalue means use the model's GA."""
    t = t.strip()
    if t == "GA":
        return None, 0
    ev, aa = None, 0
    for part in t.split("AND"):
        part = part.strip()
        if part.startswith("E<="):
            ev = float(part[3:])
        elif part.startswith("aa>"):
            aa = int(part[3:])
    if ev is None:
        raise SystemExit("cannot parse threshold %r" % t)
    return ev, aa


def main():
    from assertions import require_compute_node
    require_compute_node()          # A11: no heavy scans on a login node
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()

    anchors = read_anchor_set()
    print("anchor families in the curated set: %d" % len(anchors))

    # --- load the cached proteins -------------------------------------------
    seqs, owner, aalen = [], [], []
    accs = set()
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        name = None
        with open(os.path.join(CACHE, fn)) as fh:
            for line in fh:
                if line[0] == ">":
                    name = line[1:].rstrip("\n")
                else:
                    f = name.split("|")
                    accs.add(f[0])
                    owner.append(f[0])
                    aalen.append(int(f[5]))
                    seqs.append(pyhmmer.easel.TextSequence(
                        name=str(len(seqs)).encode(),
                        sequence=line.rstrip("\n")).digitize(alpha))
    N = len(accs)
    if not N:
        raise SystemExit("CDS cache is empty -- run 93_cds_cache.py first")
    print("cached proteins: %d over %d VirB4+ accessions" % (len(seqs), N))

    # --- scan ----------------------------------------------------------------
    out = []
    present = {}
    for d in anchors:
        ev, minaa = parse_threshold(d["threshold"])
        with pyhmmer.plan7.HMMFile(hmm_path(d["pfam_name"])) as fh:
            model = next(iter(fh))
        kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
        hit_accs, nhits = set(), 0
        for top in pyhmmer.hmmsearch([model], seqs, cpus=a.cpus, **kw):
            for h in top:
                nm = h.name
                i = int(nm.decode() if isinstance(nm, bytes) else nm)
                if minaa and aalen[i] <= minaa:
                    continue
                hit_accs.add(owner[i])
                nhits += 1
        key = "%s/%s" % (d["virb"], d["pfam_name"])
        present[key] = hit_accs
        out.append({"anchor": d["virb"], "pfam_name": d["pfam_name"],
                    "pfam_acc": d["pfam_acc"], "threshold": d["threshold"],
                    "role": d["role"], "n_proteins_hit": nhits,
                    "n_plasmids_with_hit": len(hit_accs),
                    "detection_pct": round(100.0 * len(hit_accs) / N, 2),
                    "failure_pct": round(100.0 * (N - len(hit_accs)) / N, 2)})
        print("  %-6s %-16s %-20s hit on %5d/%5d  failure %5.1f%%"
              % (d["virb"], d["pfam_name"], d["threshold"], len(hit_accs), N,
                 100.0 * (N - len(hit_accs)) / N))

    dest = os.path.join(PROJ, "data", "anchors", "anchor_ladder.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(out)

    # --- the ladder ----------------------------------------------------------
    print("\n" + "=" * 74)
    print("LADDER -- failure rate over %d VirB4+ plasmids, ranked" % N)
    print("=" * 74)
    for r in sorted(out, key=lambda x: -x["failure_pct"]):
        bar = "#" * int(r["failure_pct"] / 2)
        print("  %-6s %-16s %-24s %6.2f%%  %s"
              % (r["anchor"], r["pfam_name"], r["role"], r["failure_pct"], bar))

    ann = [r for r in out if r["role"].startswith("annotation")]
    if ann:
        f = sorted(r["failure_pct"] for r in ann)
        med = f[len(f) // 2]
        print("\n  annotation-only anchors (n=%d): median failure %.2f%%, "
              "range %.2f-%.2f%%" % (len(ann), med, f[0], f[-1]))
        vb5 = [r for r in out if r["anchor"] == "VirB5"]
        if vb5 and med > 0:
            print("  VirB5 %.2f%% = %.1fx the annotation-anchor median"
                  % (vb5[0]["failure_pct"], vb5[0]["failure_pct"] / med))

    # --- how many anchors does a typical VirB4+ plasmid carry? --------------
    per = collections.Counter()
    virb_only = {k: v for k, v in present.items() if k.startswith("VirB")}
    for acc in accs:
        per[sum(1 for s in virb_only.values() if acc in s)] += 1
    print("\n  anchors detected per plasmid (of %d VirB models):" % len(virb_only))
    for k in sorted(per):
        print("    %2d anchors  %5d plasmids  (%.1f%%)"
              % (k, per[k], 100.0 * per[k] / N))
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
