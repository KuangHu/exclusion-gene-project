#!/usr/bin/env python3
"""Stratify the 9,412 unnamed slot occupants BEFORE any functional screen.

The number 9,412 has no interpretable meaning until it is split. Five things live
in it:

  1. known genes the HMM missed        -> for THIS pipeline this is category 2:
     we call genes ab initio (Pyrodigal) and never used /gene qualifiers, so the
     R100-traS failure mode (RefSeq drops /gene, keeps /product) is not ours.
     PLSDB ships no GenBank/GFF at all -- there is no /product field to scan.
  2. distant members of known families -> THIS SCRIPT MEASURES THIS
  3. non-exclusion genes in the slot   -> pilus accessory, T4SS assembly, regulators
  4. annotation artefacts              -> ORF-calling errors, pseudogenes
  5. genuinely new families            -> what we want

TWO MEASUREMENTS, both HMM-only and immediately runnable:

A. RELAXED BACK-SCAN. The six named-family HMMs against the unnamed occupants at
   progressively looser cutoffs (GA, then E<=1e-3, 1e-2, 0.1). Whatever is caught
   is category 2 and leaves the candidate pool.

B. HOLD-OUT SENSITIVITY -- the diagnostic that makes A interpretable. For each
   named family, cluster its members, hold one cluster out, build the HMM from
   the rest, and measure recovery of the held-out cluster. This is the detection
   sensitivity to DISTANT members, measured rather than assumed.

   Why it matters: TrbK is 41% identical between RP4 and R751 while the rest of
   Tra2 is 75-92%. If hold-out recovery is only 60%, then ~40% of the remaining
   unnamed pool is still known-family members and the candidate count must be
   discounted accordingly.

   Using the SAME HMM that assigned a family to test recovery of that family is
   circular and gives 100% trivially. Leave-cluster-out is the non-circular form.
"""
import collections, csv, math, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
EEX = {"TIGR04359": "TrbK_RP4", "NF033894": "Eex_IncN", "NF041429": "EexR",
       "NF033891": "ExcA", "PF10624": "TraS", "PF14729": "DUF4467"}
CUTS = [("GA", None), ("E<=1e-5", 1e-5), ("E<=1e-3", 1e-3),
        ("E<=1e-2", 1e-2), ("E<=0.1", 0.1), ("E<=1", 1.0)]


def main():
    require_compute_node()
    import pyhmmer, igraph
    alpha = pyhmmer.easel.Alphabet.amino()
    bg = pyhmmer.plan7.Background(alpha)
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    rows = [r for r in csv.DictReader(open(SOF), delimiter="\t")
            if r["slot_status"] in ("eex_occupied", "candidate")]
    want = {}
    for r in rows:
        c = r["slot_occupant_coords"]
        if c: want.setdefault(r["accession"], set()).add(c)
    print("slot rows: %d (named %d, unnamed %d)"
          % (len(rows), sum(1 for r in rows if r["slot_status"] == "eex_occupied"),
             sum(1 for r in rows if r["slot_status"] == "candidate")), flush=True)

    seq_of = {}
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                if acc not in want: continue
                key = "%s..%s" % (f[2], f[3])
                if key in want[acc]:
                    seq_of[(acc, key)] = line.rstrip("\n")
    named, unnamed = [], []
    fam_of = {}
    for r in rows:
        s = seq_of.get((r["accession"], r["slot_occupant_coords"]))
        if not s: continue
        if r["slot_status"] == "eex_occupied":
            named.append(s); fam_of[s] = r["slot_occupant_family"].split(";")[0]
        else:
            unnamed.append(s)
    un = sorted(set(unnamed)); na = sorted(set(named))
    print("resolved sequences: named %d unique (%d records) | unnamed %d unique (%d records)"
          % (len(na), len(named), len(un), len(unnamed)), flush=True)

    dig = lambda L: [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
                     for i, s in enumerate(L)]
    U = dig(un)

    if os.environ.get("SKIP_A"): print("\n(part A skipped -- already measured)")
    print("\n=== A. RELAXED BACK-SCAN: the six named HMMs vs %d unique unnamed ===" % len(un))
    print("  %-10s %s" % ("cutoff", "  ".join("%-10s" % v for v in EEX.values())))
    caught = collections.defaultdict(set)
    for lab, ev in CUTS:
        line = "  %-10s" % lab
        for acc, fam in EEX.items():
            p = os.path.join(HMM, acc + ".hmm")
            if not os.path.exists(p): line += " %-11s" % "-"; continue
            with pyhmmer.plan7.HMMFile(p) as fh: m = next(iter(fh))
            kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
            got = set()
            for top in pyhmmer.hmmsearch([m], U, cpus=16, **kw):
                for h in top: got.add(int(dec(h.name)))
            caught[lab] |= got
            line += " %-11d" % len(got)
        print(line + "  | union %d (%.1f%%)"
              % (len(caught[lab]), 100.0*len(caught[lab])/len(un)), flush=True)
    print("\n  Anything caught here is a DISTANT MEMBER of a known family, not a new one.")

    print("\n=== B. HOLD-OUT SENSITIVITY (leave-cluster-out, non-circular) ===")
    byfam = collections.defaultdict(list)
    for s in na: byfam[fam_of[s]].append(s)
    print("  %-11s %6s %8s %9s %10s %10s" % ("family", "uniq", "clusters", "held-out",
                                             "rec@1e-5", "rec@1e-3"))
    summary = []
    for fam, mem in sorted(byfam.items(), key=lambda x: -len(x[1])):
        mem = sorted(set(mem))
        if len(mem) < 20: 
            print("  %-11s %6d   (too few to hold out)" % (fam, len(mem))); continue
        D = dig(mem)
        best = {}
        for top in pyhmmer.phmmer(D, D, cpus=16, E=10.0):
            i = int(dec(top.query.name))
            for h in top:
                j = int(dec(h.name))
                if i == j: continue
                k = (min(i, j), max(i, j))
                if k not in best or h.evalue < best[k]: best[k] = h.evalue
        g = igraph.Graph(n=len(mem), edges=[(u, v) for (u, v), e in best.items() if e <= 1e-5])
        part = g.community_leiden(objective_function="CPM", resolution=0.01, n_iterations=-1)
        cl = collections.defaultdict(list)
        for i, c in enumerate(part.membership): cl[c].append(i)
        cl = {k: v for k, v in cl.items() if len(v) >= 2}
        if len(cl) < 2:
            print("  %-11s %6d %8d   (only one cluster; no hold-out possible)"
                  % (fam, len(mem), len(cl))); continue
        rec_ga = rec_e3 = tot = 0
        for c, idxs in sorted(cl.items(), key=lambda x: -len(x[1]))[:6]:
            keep = [mem[i] for i in range(len(mem)) if i not in set(idxs)]
            hold = [mem[i] for i in idxs]
            if len(keep) < 5: continue
            q = pyhmmer.easel.TextSequence(name=b"m", sequence=max(keep, key=len)).digitize(alpha)
            b = pyhmmer.plan7.Builder(alpha)
            hmm, _, _ = b.build(q, bg)
            ds = dig(keep)
            try:
                msa = pyhmmer.hmmer.hmmalign(hmm, ds, digitize=True)
                msa.name = b"ho"; hmm, _, _ = b.build_msa(msa, bg)
            except Exception: pass
            H = dig(hold)
            # A freshly built HMM has NO gathering cutoff, so bit_cutoffs="gathering"
            # raises MissingCutoffs. The hold-out is therefore scored at E-value
            # thresholds only: 1e-5 as the strict setting (the clustering threshold
            # used throughout this project) and 1e-3 as the relaxed one used in
            # part A. "rec_ga" is retained as the 1e-5 column.
            for ev, which in ((1e-5, "strict"), (1e-3, "relaxed")):
                got = set()
                for top in pyhmmer.hmmsearch([hmm], H, cpus=16, E=ev):
                    for h in top: got.add(int(dec(h.name)))
                if which == "strict": rec_ga += len(got)
                else: rec_e3 += len(got)
            tot += len(hold)
        if tot:
            print("  %-11s %6d %8d %9d %9.1f%% %9.1f%%"
                  % (fam, len(mem), len(cl), tot, 100.0*rec_ga/tot, 100.0*rec_e3/tot))
            summary.append((fam, tot, rec_ga, rec_e3))
    if summary:
        T = sum(x[1] for x in summary); G = sum(x[2] for x in summary); E = sum(x[3] for x in summary)
        print("\n  OVERALL hold-out recovery: E<=1e-5 %.1f%%  |  E<=1e-3 %.1f%%  (n=%d held out)"
              % (100.0*G/T, 100.0*E/T, T))
        print("\n  READING: recovery is the sensitivity to DISTANT family members.")
        print("  If it is R%%, then roughly (100-R)%% of true known-family members in the")
        print("  unnamed pool are still being missed, and the candidate count must be")
        print("  discounted by that much BEFORE any 'new family' claim.")
    dest = os.path.join(PROJ, "data", "anchors", "unnamed_stratification.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cutoff", "n_unnamed_caught", "pct_of_unnamed"])
        for lab, _ in CUTS:
            w.writerow([lab, len(caught[lab]), round(100.0*len(caught[lab])/len(un), 2)])
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
