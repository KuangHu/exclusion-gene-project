#!/usr/bin/env python3
"""Full Pfam-A scan of the unnamed remainder -- the test that should have come first.

"Unnamed" has so far meant only: no hit to the SIX exclusion families. It has
never meant no hit to ANYTHING. A slot occupant can be a well-characterised pilus
accessory, T4SS assembly factor or regulator with a perfectly good Pfam domain,
and nothing in this pipeline has yet asked.

That is category 3, and Pfam-A (30,134 models, local) answers it directly.
InterProScan is not installed, but its Pfam component is what would do the work
anyway.

Direction: hmmscan (each sequence against the whole model library), not hmmsearch,
because the query set is small and the library is large.

Gathering thresholds only (A13).

This bears hardest on cluster 22 -- 397 unique sequences over 1,196 independent
Mash clusters, median 176 aa. That prevalence is wrong for an exclusion gene:
TrbK's mature form is 47 aa and Eex ~75 aa, and an exclusion gene is not the
DEFAULT occupant of its own slot. A conserved assembly component fits better, and
if so it will carry a Pfam domain.
"""
import collections, csv, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
PFAM = "/global/scratch/users/kh36969/funcannot_dbs/pfam/Pfam-A.hmm"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"]


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    rows = [r for r in csv.DictReader(open(SOF), delimiter="\t")
            if r["slot_status"] == "candidate" and r["slot_occupant_coords"]]
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
    acc_of = collections.defaultdict(set)
    for r in rows:
        s = seq_of.get((r["accession"], r["slot_occupant_coords"]))
        if s: acc_of[s].add(r["accession"])
    uniq = sorted(acc_of)
    dig = lambda L: [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
                     for i, s in enumerate(L)]
    U = dig(uniq); known = set()
    for f in EEX:
        p = os.path.join(HMM, f + ".hmm")
        if not os.path.exists(p): continue
        with pyhmmer.plan7.HMMFile(p) as fh: m = next(iter(fh))
        for top in pyhmmer.hmmsearch([m], U, cpus=8, E=1e-3):
            for h in top: known.add(int(dec(h.name)))
    rem = [uniq[i] for i in range(len(uniq)) if i not in known]
    print("remainder to scan against Pfam-A: %d unique sequences" % len(rem), flush=True)

    Q = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
         for i, s in enumerate(rem)]
    hit = collections.defaultdict(list)
    with pyhmmer.plan7.HMMFile(PFAM) as fh:
        n = 0
        for top in pyhmmer.hmmer.hmmscan(Q, fh, cpus=8, bit_cutoffs="gathering"):
            i = int(dec(top.query.name))
            for h in top:
                hit[i].append((dec(h.name), dec(h.accession or b""), round(h.score, 1)))
            n += 1
            if n % 500 == 0: print("  scanned %d/%d" % (n, len(rem)), flush=True)
    withhit = len(hit)
    print("\n=== Pfam-A coverage of the unnamed remainder ===")
    print("  sequences with >=1 Pfam domain at GA : %d / %d (%.1f%%)"
          % (withhit, len(rem), 100.0*withhit/len(rem)))
    print("  with NO Pfam domain at all           : %d (%.1f%%)"
          % (len(rem)-withhit, 100.0*(len(rem)-withhit)/len(rem)))
    fam = collections.Counter(h[0] for v in hit.values() for h in v)
    print("\n=== most common Pfam domains found ===")
    print("  %-24s %-12s %7s" % ("domain", "accession", "n_seqs"))
    acc_of_fam = {}
    for v in hit.values():
        for nm2, ac, _ in v: acc_of_fam[nm2] = ac
    for nm2, c in fam.most_common(20):
        print("  %-24s %-12s %7d" % (nm2, acc_of_fam.get(nm2, ""), c))
    dest = os.path.join(PROJ, "data", "anchors", "pfam_remainder.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["seq_idx", "aa_len", "n_records", "pfam_domains"])
        for i, s in enumerate(rem):
            w.writerow([i, len(s), len(acc_of[s]),
                        ";".join("%s(%s)" % (h[0], h[2]) for h in hit.get(i, []))])
    print("\nwrote %s" % dest)
    print("\n  Sequences with a Pfam domain are NOT new families -- they are")
    print("  characterised proteins that happen to sit in the slot (category 3).")


if __name__ == "__main__":
    main()
