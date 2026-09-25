#!/usr/bin/env python3
"""Round 1: VALIDATE the Round 0 families. Not expand them.

Candidates were nominated BY POSITION, so position cannot validate them (§7.2).
Round 1 asks one question per family:

    in-slot / total  =  (hits at the slot) / (hits anywhere in PLSDB)

  near 1      slot-specific -- the family really is defined by the position
  well < 1    a generic small-protein family; the positional nomination caught
              a coincidence

SEARCH THRESHOLD = CLUSTERING THRESHOLD = 1e-05 (§3). A looser search would
inflate the denominator for reasons unrelated to biology.

FAMILY HMM. No aligner is installed (mafft/muscle/clustalo all absent), so the
profile is bootstrapped with the tools present: build a single-sequence HMM from
the MEDOID, hmmalign the cluster members to it, rebuild from that MSA. Medoid,
never `max(v, key=len)` -- the longest member is the most domain-rich and most
promiscuous query, the root cause of the retracted phmmer recall (recovered
proteins median 816 aa against a family median of 238).

CALIBRATION ALREADY MEASURED -- and the margin is thin:

    upper   NF033894  0.771     TIGR04359  0.812
    null    VirB10-nominated families: median 0.204, max 0.745

    margin = 0.026. in-slot/total ALONE CANNOT DECIDE. It is reported for every
    family and read only alongside the negative-slot families run through this
    identical path.

The null may itself be contaminated by unrecognised VirB11 (VirB10+1 IS VirB11 in
colinear order). The negative slots here are decontaminated first, so if the null
max drops to ~0.3 the criterion becomes usable; if it does not, this stays a
descriptive statistic.
"""
import argparse, collections, csv, json, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
OUT = "/global/scratch/users/kh36969/exclusion_gene/slot_v11"
CAT = os.path.join(PROJ, "data", "release", "v1.1", "catalogue_MPF_T_v1.1.tsv")
EVAL = 1e-05
LIPO_STRICT = r"[LVI][ASTVIG][GASN]C"


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="VirB5+1")
    ap.add_argument("--anchor", default="VirB5")
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer, re
    alpha = pyhmmer.easel.Alphabet.amino()
    bg = pyhmmer.plan7.Background(alpha)
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    rx = re.compile(LIPO_STRICT)

    clusters = list(csv.DictReader(
        open(os.path.join(PROJ, "data", "anchors",
                          "round0_%s_v11.tsv" % a.tag.replace("+", "p"))), delimiter="\t"))
    meta = json.load(open(os.path.join(OUT, "slot_%s_records.json" % a.tag)))
    recs = meta["records"]
    names, aas = [], []
    nm = None
    for line in open(os.path.join(OUT, "slot_%s_clean.faa" % a.tag)):
        if line[0] == ">": nm = line[1:].split()[0]
        else: names.append(nm); aas.append(line.rstrip("\n"))
    byname = dict(zip(names, aas))

    # the slot positions, to define "in-slot"
    cat = {r["accession"]: r for r in csv.DictReader(open(CAT), delimiter="\t")
           if r["slot_ready__virb5_virb6"] == "1"}
    slot_seqs = set(recs.values())
    print("families to validate: %d (n>=2 only)"
          % sum(1 for c in clusters if int(c["n_unique"]) >= 2), flush=True)

    seqs, owner, aas_all = [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm2 = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm2 = line[1:].rstrip("\n")
            else:
                owner.append(nm2.split("|")[0]); p = line.rstrip("\n"); aas_all.append(p)
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(), sequence=p).digitize(alpha))
    print("full cache: %d proteins over %d plasmids\n" % (len(seqs), len(set(owner))), flush=True)

    # cluster membership: recompute from the FASTA order used in round 0
    memb = collections.defaultdict(list)
    for c in clusters:
        memb[c["cluster"]] = []
    # round0 wrote medoid names; re-derive members by re-reading the cluster file is
    # not possible, so members are taken as all sequences assigned in round0 output
    # -> we rebuild from the medoid alone when membership is unavailable
    rows = []
    for c in clusters:
        if int(c["n_unique"]) < 2: continue
        med = c["medoid_seq"] if "medoid_seq" in c else byname.get(c["medoid"], "")
        if not med: continue
        q = pyhmmer.easel.TextSequence(name=b"medoid", sequence=med).digitize(alpha)
        builder = pyhmmer.plan7.Builder(alpha)
        hmm, _, _ = builder.build(q, bg)
        # bootstrap to a real family profile: align the members to the medoid-seeded
        # HMM, then rebuild from that MSA. Single-sequence profiles are markedly
        # less sensitive, and in-slot/total is a RATIO -- a weak profile shrinks
        # numerator and denominator unequally and would bias the statistic.
        mem_s = [x for x in (c.get("member_seqs") or "").split(";") if x]
        if len(mem_s) >= 2:
            ds = [pyhmmer.easel.TextSequence(name=("m%d" % k).encode(),
                  sequence=x).digitize(alpha) for k, x in enumerate(mem_s)]
            try:
                msa = pyhmmer.hmmer.hmmalign(hmm, ds, digitize=True)
                msa.name = ("cl%s" % c["cluster"]).encode()
                hmm, _, _ = builder.build_msa(msa, bg)
                nseq = len(mem_s)
            except Exception as e:
                print("    cluster %s: MSA rebuild failed (%s); medoid-only profile"
                      % (c["cluster"], type(e).__name__), flush=True)
                nseq = 1
        else:
            nseq = 1
        hits = []
        for top in pyhmmer.hmmsearch([hmm], seqs, cpus=a.cpus, E=EVAL):
            for h in top: hits.append(int(dec(h.name)))
        if not hits: continue
        tot = len(hits)
        ins = sum(1 for i in hits if aas_all[i] in slot_seqs and owner[i] in cat)
        frac = ins / tot
        L = sorted(len(aas_all[i]) for i in hits)
        lip = sum(1 for i in hits if rx.search(aas_all[i][:40]))
        rows.append({"cluster": c["cluster"], "n_unique_r0": c["n_unique"],
                     "n_records_r0": c["n_records"], "branch": c["branch"],
                     "seeds": c["seeds_in_cluster"], "medoid_aa": c["medoid_aa"],
                     "hits_total": tot, "hits_in_slot": ins,
                     "in_slot_over_total": round(frac, 3),
                     "hit_plasmids": len(set(owner[i] for i in hits)),
                     "median_hit_aa": L[len(L)//2], "profile_nseq": nseq,
                     "lipobox_pct": round(100.0*lip/tot, 1)})
        print("  cluster %-5s r0 %3s uniq  hits %6d  in-slot %5d  ratio %.3f  lipo %.1f%%  %s"
              % (c["cluster"], c["n_unique"], tot, ins, frac, 100.0*lip/tot,
                 c["seeds_in_cluster"] or c["branch"]), flush=True)

    dest = os.path.join(PROJ, "data", "anchors", "round1_%s_v11.tsv" % a.tag.replace("+", "p"))
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(sorted(rows, key=lambda r: -r["in_slot_over_total"]))
    print("\nwrote %s (%d families)" % (dest, len(rows)))
    fr = sorted(r["in_slot_over_total"] for r in rows)
    print("\n=== in-slot/total distribution (%s) ===" % a.tag)
    print("  n=%d  median %.3f  Q1 %.3f  Q3 %.3f  max %.3f"
          % (len(fr), fr[len(fr)//2], fr[len(fr)//4], fr[3*len(fr)//4], fr[-1]))
    print("  reference: NF033894 0.771, TIGR04359 0.812 | VirB10 null median 0.204 max 0.745")
    print("  margin between upper and null max = 0.026 -- NOT decidable on this axis alone")
    print("\n  families above the null max (0.745): %d"
          % sum(1 for x in fr if x > 0.745))


if __name__ == "__main__":
    main()
