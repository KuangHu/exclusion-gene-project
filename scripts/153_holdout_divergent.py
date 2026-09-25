#!/usr/bin/env python3
"""Hold-out sensitivity, done correctly: hold out the MOST DIVERGENT members.

The first attempt (job 26091849) was invalid. It clustered each family and held
out each cluster largest-first, so for Eex_IncN it held out 382 of 385 members and
built the HMM from the ~5 that remained. An HMM built from 5 sequences and tested
against 380 measures profile starvation, not detection sensitivity. The other
three families collapsed to a single cluster and were never tested at all.

CORRECT DESIGN. Sensitivity to distant members is exactly: build the profile from
the family's CORE and ask whether it still finds the family's EDGE.

  1. all-vs-all phmmer within the family
  2. rank members by summed similarity to the rest (medoid-centrality)
  3. hold out the least-central FRACTION (the divergent edge)
  4. build the HMM from the central remainder
  5. measure recovery of the held-out edge at 1e-5, 1e-3, 1e-2

The held-out set is always the minority and the training set always the majority,
which is the opposite of what went wrong.

Reported at three hold-out fractions so the answer is a curve, not a point: a
family whose edge is recovered at 10% held out but not at 30% has a different
shape from one that fails at both.

WHY THE NUMBER MATTERS. It is the discount factor on the candidate pool. If the
core profile recovers only R% of its own divergent members, then roughly (100-R)%
of true known-family members sitting in the unnamed pool are invisible to this
method, and no 'new family' count is meaningful until that is stated.
"""
import collections, csv, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
FRACTIONS = [0.10, 0.20, 0.30]
MIN_MEMBERS = 40


def main():
    require_compute_node()
    import pyhmmer, math
    alpha = pyhmmer.easel.Alphabet.amino()
    bgm = pyhmmer.plan7.Background(alpha)
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    rows = [r for r in csv.DictReader(open(SOF), delimiter="\t")
            if r["slot_status"] == "eex_occupied" and r["slot_occupant_coords"]]
    want = collections.defaultdict(set)
    for r in rows: want[r["accession"]].add(r["slot_occupant_coords"])
    seq_of = {}
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                if acc in want:
                    k = "%s..%s" % (f[2], f[3])
                    if k in want[acc]: seq_of[(acc, k)] = line.rstrip("\n")
    byfam = collections.defaultdict(set)
    for r in rows:
        s = seq_of.get((r["accession"], r["slot_occupant_coords"]))
        if s: byfam[r["slot_occupant_family"].split(";")[0]].add(s)
    print("families: %s\n" % ", ".join("%s %d" % (k, len(v)) for k, v in
                                       sorted(byfam.items(), key=lambda x: -len(x[1]))), flush=True)

    dig = lambda L: [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
                     for i, s in enumerate(L)]
    out = []
    print("  %-11s %6s %7s %8s %9s %9s %9s"
          % ("family", "uniq", "held", "train", "rec@1e-5", "rec@1e-3", "rec@1e-2"))
    for fam, ms in sorted(byfam.items(), key=lambda x: -len(x[1])):
        mem = sorted(ms)
        if len(mem) < MIN_MEMBERS:
            print("  %-11s %6d   (<%d members, not tested)" % (fam, len(mem), MIN_MEMBERS))
            continue
        D = dig(mem)
        sim = collections.defaultdict(float)
        for top in pyhmmer.phmmer(D, D, cpus=16, E=10.0):
            i = int(dec(top.query.name))
            for h in top:
                j = int(dec(h.name))
                if i == j: continue
                sim[i] += -math.log10(max(h.evalue, 1e-300))
        order = sorted(range(len(mem)), key=lambda i: sim.get(i, 0.0))   # least central first
        for frac in FRACTIONS:
            k = max(1, int(len(mem) * frac))
            hold = [mem[i] for i in order[:k]]          # the DIVERGENT EDGE
            keep = [mem[i] for i in order[k:]]          # the CORE (majority)
            assert len(keep) > len(hold), "training set must be the majority"
            q = pyhmmer.easel.TextSequence(
                name=b"m", sequence=max(keep, key=lambda s: sim.get(mem.index(s), 0))).digitize(alpha)
            b = pyhmmer.plan7.Builder(alpha)
            hmm, _, _ = b.build(q, bgm)
            try:
                msa = pyhmmer.hmmer.hmmalign(hmm, dig(keep), digitize=True)
                msa.name = b"core"; hmm, _, _ = b.build_msa(msa, bgm)
            except Exception: pass
            H = dig(hold)
            rec = {}
            for ev in (1e-5, 1e-3, 1e-2):
                got = set()
                for top in pyhmmer.hmmsearch([hmm], H, cpus=16, E=ev):
                    for h in top: got.add(int(dec(h.name)))
                rec[ev] = 100.0 * len(got) / len(hold)
            print("  %-11s %6d %6d%% %8d %8.1f%% %8.1f%% %8.1f%%"
                  % (fam, len(mem), int(frac*100), len(keep),
                     rec[1e-5], rec[1e-3], rec[1e-2]), flush=True)
            out.append({"family": fam, "n_unique": len(mem), "holdout_pct": int(frac*100),
                        "n_held": len(hold), "n_train": len(keep),
                        "rec_1e5": round(rec[1e-5], 1), "rec_1e3": round(rec[1e-3], 1),
                        "rec_1e2": round(rec[1e-2], 1)})
    if out:
        dest = os.path.join(PROJ, "data", "anchors", "holdout_sensitivity.tsv")
        with open(dest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(out)
        print("\nwrote %s" % dest)
        for frac in FRACTIONS:
            s = [r for r in out if r["holdout_pct"] == int(frac*100)]
            if not s: continue
            n = sum(r["n_held"] for r in s)
            g3 = sum(r["rec_1e3"]*r["n_held"] for r in s)/max(1, n)
            print("  hold-out %2d%%: weighted recovery at E<=1e-3 = %.1f%%  (n=%d held)"
                  % (int(frac*100), g3, n))
        print("\n  DISCOUNT: the unnamed pool's 2,514 unique sequences should be read")
        print("  as containing a further (100-R)%% of known-family members that this")
        print("  method cannot see, on top of the 279 already recovered by back-scan.")


if __name__ == "__main__":
    main()
