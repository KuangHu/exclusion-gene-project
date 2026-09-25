#!/usr/bin/env python3
"""MPF typing of CURATED COMPLETE phage genomes (RefSeq viral).

WHY THIS SET RATHER THAN INPHARED. INPHARED's S3 bucket returns 403 from this
cluster for both listing and direct keys. RefSeq viral is one 178 MB file over
plain HTTPS, curated, and contains the canonical controls -- lambda, T4, Mu, P22
and P1 (the phage-plasmid).

WHY THIS IS A BETTER NEGATIVE CONTROL THAN THE PROPHAGE SET. The prophage run
(job 26081820, 16/32 shards) gave 0.165% MPF-positive, and those positives were
median 123 kb with 95% relaxase -- an ICE signature. The interpretation was that
TIGER/Islander call mobile islands from integrase + att sites, which ICEs satisfy
too, so the positives are misannotated ICEs rather than conjugative phages.

THAT INTERPRETATION MAKES A PREDICTION: complete, curated phage genomes -- which
are not integrase-flanked islands excised from a chromosome -- should be at or
near ZERO.

  prophages (island calls)   0.165%   <- predicted to be ICE contamination
  RefSeq phage genomes       ~0%      <- this run tests that

If RefSeq phages come back at a similar 0.16%, the interpretation is wrong and the
signal is something about phages, not about island calling.

Same frozen criteria, same frozen Pyrodigal settings, gathering thresholds only.
"""
import collections, csv, gzip, os, re, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import orf_caller

FA = "/global/scratch/users/kh36969/phage_refseq/phage_subset.fna.gz"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"]
NAME = {"TIGR04359": "TrbK_RP4", "NF033894": "Eex_IncN", "NF041429": "EexR",
        "NF033891": "ExcA", "PF10624": "TraS", "PF14729": "DUF4467"}
MOB = ["T4SS_MOB" + x for x in ("B","C","F","H","P1","P2","P3","Q","T","V")]


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    owner, aas, ncds, slen, desc = [], [], collections.Counter(), {}, {}
    n = 0
    with gzip.open(FA, "rt") as fh:
        name, buf = None, []
        def flush():
            nonlocal name, buf
            if name and buf:
                s = "".join(buf).upper(); slen[name] = len(s)
                for g in orf_caller.call(s):
                    p = (g["aa"] or "").rstrip("*")
                    if len(p) < 30: continue
                    owner.append(name); aas.append(p); ncds[name] += 1
            name, buf = None, []
        for line in fh:
            if line[0] == ">":
                flush(); h = line[1:].rstrip("\n")
                name = h.split()[0]; desc[name] = h; buf = []; n += 1
            else: buf.append(line.strip())
        flush()
    print("phage genomes %d ; proteins %d" % (n, len(aas)), flush=True)
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
    FA_ = {a for a in vb4 if a not in (T | F | I) and nprof[a] >= 2}
    relax = set()
    for p in MOB: relax |= scan(os.path.join(CJ, p + ".hmm"))
    eex = {NAME[f]: scan(os.path.join(HMM, f + ".hmm")) for f in EEX}

    rows = []
    for acc in sorted(set(owner)):
        cls = [c for c, s in (("MPF_T", T), ("MPF_F", F), ("MPF_I", I), ("MPF_FA", FA_))
               if acc in s]
        ex = [k for k, v in eex.items() if acc in v]
        rows.append({"accession": acc, "length": slen.get(acc, ""),
                     "n_cds": ncds[acc], "mpf_classes": ";".join(cls) or "none",
                     "n_fa_fata": nprof.get(acc, 0), "relaxase": int(acc in relax),
                     "exclusion_family_hit": ";".join(ex),
                     "description": desc.get(acc, "")[:120]})
    dest = os.path.join(PROJ, "data", "anchors", "mpf_in_refseq_phage.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    nm = sum(1 for r in rows if r["mpf_classes"] != "none")
    nr = sum(r["relaxase"] for r in rows)
    ne = sum(1 for r in rows if r["exclusion_family_hit"])
    print("\n=== RESULT: curated complete phage genomes ===")
    print("  genomes            %6d" % len(rows))
    print("  MPF-positive       %6d  (%.3f%%)" % (nm, 100.0*nm/len(rows)))
    print("  relaxase-positive  %6d  (%.3f%%)" % (nr, 100.0*nr/len(rows)))
    print("  exclusion family   %6d  (%.3f%%)" % (ne, 100.0*ne/len(rows)))
    print("\n=== the prediction this run tests ===")
    print("  prophage island calls   0.165%%  (attributed to misannotated ICEs)")
    print("  RefSeq phage genomes    %.3f%%" % (100.0*nm/len(rows)))
    print("  -> %s" % ("CONSISTENT: complete phage genomes are far below the island calls,"
                       " supporting the ICE-contamination reading"
                       if 100.0*nm/len(rows) < 0.08 else
                       "NOT consistent: the signal is not explained by island calling alone"))
    if nm:
        print("\n=== the MPF-positive phage genomes ===")
        for r in rows:
            if r["mpf_classes"] != "none":
                print("  %-14s %7s bp  %-12s relax=%d  %s"
                      % (r["accession"], r["length"], r["mpf_classes"],
                         r["relaxase"], r["description"][:70]))
    print("\n=== canonical controls ===")
    for pat in ("phage lambda", "phage P1", "phage T4", "phage Mu", "phage P22"):
        for r in rows:
            if re.search(pat, r["description"], re.I):
                print("  %-16s %-14s %7s bp  MPF=%s relax=%d"
                      % (pat, r["accession"], r["length"], r["mpf_classes"], r["relaxase"]))
                break
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
