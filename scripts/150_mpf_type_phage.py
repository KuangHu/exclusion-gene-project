#!/usr/bin/env python3
"""Type the four MPF classes in prophages. One SLURM array task per shard.

Same frozen criteria and same frozen Pyrodigal settings as the plasmid release
and the ICE run, so all three datasets are directly comparable.

PROPHAGES ARE THE NEGATIVE CONTROL. Phages are not conjugative; a prophage has no
reason to carry an MPF system. Expected near-zero.

A non-zero result is informative in a specific way rather than a discovery:

  * TIGER/Islander call mobile islands by integrase + att sites. An ICE satisfies
    those signals too, so an MPF-positive "prophage" is most likely a MISANNOTATED
    ICE -- this search doubles as QC on the prophage set.
  * Genuine phage-plasmids (P1-like) and a few phage-inducible islands do carry
    conjugation or transfer modules.

So the reading is NOT "prophages are conjugative". It is: how much of the prophage
set is actually something else, and does it concentrate where an ICE would.

Reported per shard; combine across shards before drawing any rate.
"""
import collections, csv, gzip, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import orf_caller

SEQ = "/global/scratch/users/kh36969/phage_SSAP_SSB/seq_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
OUT = "/global/scratch/users/kh36969/exclusion_gene/mge_search/phage"
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"]
NAME = {"TIGR04359": "TrbK_RP4", "NF033894": "Eex_IncN", "NF041429": "EexR",
        "NF033891": "ExcA", "PF10624": "TraS", "PF14729": "DUF4467"}


def main():
    require_compute_node()
    shard = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    os.makedirs(OUT, exist_ok=True)
    fa = os.path.join(SEQ, "phage.shard%03d.fna.gz" % shard)
    if not os.path.exists(fa):
        print("no such shard: %s" % fa); return 1

    owner, aas, ncds, slen = [], [], collections.Counter(), {}
    n = 0
    with gzip.open(fa, "rt") as fh:
        name, buf = None, []
        def flush():
            nonlocal name, buf
            if name and buf:
                s = "".join(buf).upper()
                slen[name] = len(s)
                for g in orf_caller.call(s):
                    p = (g["aa"] or "").rstrip("*")
                    if len(p) < 30: continue
                    owner.append(name); aas.append(p); ncds[name] += 1
            name, buf = None, []
        for line in fh:
            if line[0] == ">":
                flush(); name = line[1:].split()[0]; buf = []; n += 1
            else:
                buf.append(line.strip())
        flush()
    print("shard %03d: %d sequences, %d proteins" % (shard, n, len(aas)), flush=True)
    seqs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=p).digitize(alpha)
            for i, p in enumerate(aas)]

    def scan(path):
        if not os.path.exists(path): return set()
        with pyhmmer.plan7.HMMFile(path) as fh: m = next(iter(fh))
        got = set()
        for top in pyhmmer.hmmsearch([m], seqs, cpus=16, bit_cutoffs="gathering"):
            for h in top: got.add(owner[int(dec(h.name))])
        return got

    T = scan(os.path.join(HMM, "anchor_CagE_TrbE_VirB.hmm"))
    F = scan(os.path.join(HMM, "mpff_TraC_F_IV.hmm"))
    I = scan(os.path.join(CJ, "T4SS_I_traU.hmm"))
    vb4 = scan(os.path.join(CJ, "T4SS_virb4.hmm"))
    faf = sorted(f[:-4] for f in os.listdir(CJ)
                 if (f.startswith("FA_") or f.startswith("FATA_")) and f.endswith(".hmm"))
    nprof = collections.Counter()
    for p in faf:
        for a in scan(os.path.join(CJ, p + ".hmm")): nprof[a] += 1
    prev = T | F | I
    FA = {a for a in vb4 if a not in prev and nprof[a] >= 2}
    eex = {NAME[f]: scan(os.path.join(HMM, f + ".hmm")) for f in EEX}
    relax = set()
    for p in ("T4SS_MOBB","T4SS_MOBC","T4SS_MOBF","T4SS_MOBH","T4SS_MOBP1",
              "T4SS_MOBP2","T4SS_MOBP3","T4SS_MOBQ","T4SS_MOBT","T4SS_MOBV"):
        relax |= scan(os.path.join(CJ, p + ".hmm"))

    rows = []
    for acc in sorted(set(owner)):
        cls = [c for c, s in (("MPF_T", T), ("MPF_F", F), ("MPF_I", I), ("MPF_FA", FA))
               if acc in s]
        ex = [k for k, v in eex.items() if acc in v]
        if not cls and not ex and acc not in relax:
            continue                      # keep the output to informative rows
        rows.append({"phage_id": acc, "shard": shard, "length": slen.get(acc, ""),
                     "n_cds_called": ncds[acc],
                     "mpf_classes": ";".join(cls) or "none",
                     "n_fa_fata_profiles": nprof.get(acc, 0),
                     "relaxase": int(acc in relax),
                     "exclusion_family_hit": ";".join(ex)})
    dest = os.path.join(OUT, "phage_mpf.shard%03d.tsv" % shard)
    cols = ["phage_id","shard","length","n_cds_called","mpf_classes",
            "n_fa_fata_profiles","relaxase","exclusion_family_hit"]
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    nm = sum(1 for r in rows if r["mpf_classes"] != "none")
    print("  MPF-positive %d / %d = %.3f%%" % (nm, n, 100.0*nm/max(1, n)))
    print("  relaxase %d | exclusion-family %d"
          % (sum(r["relaxase"] for r in rows),
             sum(1 for r in rows if r["exclusion_family_hit"])))
    print("  classes: %s" % dict(collections.Counter(
        r["mpf_classes"] for r in rows if r["mpf_classes"] != "none")))
    print("wrote %s (%d informative rows of %d sequences)" % (dest, len(rows), n))
    with open(os.path.join(OUT, "shard%03d.n" % shard), "w") as fh:
        fh.write("%d\t%d\n" % (n, nm))


if __name__ == "__main__":
    sys.exit(main() or 0)
