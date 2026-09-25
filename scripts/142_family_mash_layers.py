#!/usr/bin/env python3
"""Three-layer counts (§8) for the lipobox>90% families: record / Mash cluster / unique.

Round 1 reported in-slot RECORDS but never the accessions, so the Mash layer was
missing and §8 requires all three. Until this runs, 514 in-slot records cannot be
quoted: IncP cloning is unmeasured, and IncI2 already showed 181 records
collapsing to 25 Mash clusters.

Mash clusters: single linkage on pairs at d <= 0.007, from pairs_d010.tsv.gz.
That threshold achieves DEDUPLICATION ONLY, not independence -- the graph
percolates at d = 0.016-0.018, so any biologically meaningful threshold would
connect the whole database. Cluster-level hypergeometric tests still violate
independence; cross-family enrichment is valid only at PTU/Inc-group level.

cl39 gets particular attention: 28 hits, barely over the >=20 line, and this run
already established that ratio 1.000 is trivial on small families (the 12 target
families at ratio 1.0 have median hits_total = 5).
"""
import collections, csv, gzip, json, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
OUT = "/global/scratch/users/kh36969/exclusion_gene/slot_v11"
PAIRS = "/global/scratch/users/kh36969/exclusion_gene/mash/pairs_d010.tsv.gz"
CAT = os.path.join(PROJ, "data", "release", "v1.1", "catalogue_MPF_T_v1.1.tsv")
D = 0.007
LIPO_CUT = 90.0
MIN_HITS = 20


def main():
    require_compute_node()
    import pyhmmer, re
    alpha = pyhmmer.easel.Alphabet.amino()
    bg = pyhmmer.plan7.Background(alpha)
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    rx = re.compile(r"[LVI][ASTVIG][GASN]C")

    r0 = {c["cluster"]: c for c in csv.DictReader(
        open(os.path.join(PROJ, "data", "anchors", "round0_VirB5p1_v11.tsv")), delimiter="\t")}
    r1 = [x for x in csv.DictReader(
        open(os.path.join(PROJ, "data", "anchors", "round1_VirB5p1_v11.tsv")), delimiter="\t")
        if int(x["hits_total"]) >= MIN_HITS and float(x["lipobox_pct"]) > LIPO_CUT]
    print("families at lipobox>%.0f%%, hits>=%d: %s"
          % (LIPO_CUT, MIN_HITS, ", ".join("cl" + x["cluster"] for x in r1)), flush=True)

    meta = json.load(open(os.path.join(OUT, "slot_VirB5+1_records.json")))
    recs = meta["records"]                      # accession -> slot occupant sequence
    slot_seqs = set(recs.values())
    cat = {l.split("\t")[0] for i, l in enumerate(open(CAT)) if i}

    seqs, owner, aas = [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                owner.append(nm.split("|")[0]); p = line.rstrip("\n"); aas.append(p)
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(), sequence=p).digitize(alpha))
    print("cache: %d proteins\n" % len(seqs), flush=True)

    fam_acc, fam_seq = {}, {}
    for x in r1:
        c = r0[x["cluster"]]
        med = c["medoid_seq"]
        q = pyhmmer.easel.TextSequence(name=b"m", sequence=med).digitize(alpha)
        b = pyhmmer.plan7.Builder(alpha)
        hmm, _, _ = b.build(q, bg)
        mem = [s for s in (c.get("member_seqs") or "").split(";") if s]
        if len(mem) >= 2:
            ds = [pyhmmer.easel.TextSequence(name=("m%d" % k).encode(),
                  sequence=s).digitize(alpha) for k, s in enumerate(mem)]
            try:
                msa = pyhmmer.hmmer.hmmalign(hmm, ds, digitize=True)
                msa.name = ("cl%s" % x["cluster"]).encode()
                hmm, _, _ = b.build_msa(msa, bg)
            except Exception:
                pass
        acc, sq = set(), set()
        for top in pyhmmer.hmmsearch([hmm], seqs, cpus=16, E=1e-05):
            for h in top:
                i = int(dec(h.name))
                if aas[i] in slot_seqs and owner[i] in cat:
                    acc.add(owner[i]); sq.add(aas[i])
        fam_acc[x["cluster"]] = acc; fam_seq[x["cluster"]] = sq
        print("  cl%-5s in-slot records %4d  unique seqs %4d"
              % (x["cluster"], len(acc), len(sq)), flush=True)

    need = set().union(*fam_acc.values()) if fam_acc else set()
    print("\nloading Mash pairs at d<=%.3f for %d accessions ..." % (D, len(need)), flush=True)
    par = {a: a for a in need}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x

    kept = 0
    with gzip.open(PAIRS, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3: continue
            try:
                if float(f[2]) > D: continue
            except ValueError:
                continue
            a, b2 = f[0], f[1]
            if a in par and b2 in par:
                ra, rb = find(a), find(b2)
                if ra != rb: par[ra] = rb
                kept += 1
    print("  edges at d<=%.3f among these accessions: %d" % (D, kept))

    print("\n=== THREE-LAYER COUNTS (§8) ===")
    print("  %-8s %9s %9s %9s  %s" % ("family", "records", "Mash cl", "unique", "collapse"))
    rows = []
    for x in r1:
        c = x["cluster"]
        acc = fam_acc[c]
        cl = len({find(a) for a in acc}) if acc else 0
        u = len(fam_seq[c])
        rows.append({"cluster": c, "in_slot_records": len(acc), "mash_clusters": cl,
                     "unique_seqs": u, "lipobox_pct": x["lipobox_pct"],
                     "in_slot_over_total": x["in_slot_over_total"],
                     "hits_total": x["hits_total"], "median_hit_aa": x["median_hit_aa"],
                     "seeds": x["seeds"]})
        print("  cl%-6s %9d %9d %9d  %.1fx  %s"
              % (c, len(acc), cl, u, len(acc) / max(1, cl), x["seeds"] or ""))
    allacc = set().union(*fam_acc.values()) if fam_acc else set()
    print("\n  UNION over the %d families: %d records / %d Mash clusters / %d unique"
          % (len(r1), len(allacc), len({find(a) for a in allacc}) if allacc else 0,
             len(set().union(*fam_seq.values())) if fam_seq else 0))
    dest = os.path.join(PROJ, "data", "anchors", "family_three_layer_v11.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
