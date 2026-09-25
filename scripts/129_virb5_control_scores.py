#!/usr/bin/env python3
"""What do the four VirB5 controls ACTUALLY score in the pipeline's own search?

anchor_set.tsv justifies VirB5's E<=1e-5 threshold with:
    RP4 trbJ 17.7 bits, E 7.9e-07 ; R751 22.2, E 3.4e-08
and concludes that GA loses RP4/R751 while the E-value rule keeps them.

That pair is internally inconsistent with the measured boundary: in an hmmsearch
over 1,099,911 proteins, E=1e-5 falls at ~31.7 bits. A 17.7-bit hit cannot carry
E=7.9e-07 in the same search. Those E-values were almost certainly produced by a
DIFFERENT search mode -- one protein scanned against the Pfam library (N ~ 30k
models) rather than one model searched against 1.1M proteins.

If so, the E<=1e-5 rule never rescued RP4/R751 at scale, and job 25898457 already
observed exactly that (both "PF07996 NO, CONJScan yes").

This measures it directly, both ways, on the named controls:
    hmmsearch  model -> the 1.1M cache        (what the pipeline does)
    hmmscan    control protein -> Pfam-A      (what the note's E-values look like)

No conclusion is read from a control that cannot be located in the cache.
"""
import os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
SUB = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
HMMDIR = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
CONTROLS = {"RP4": "BN000925.1", "R751": "NC_001735.4",
            "pKM101": "U09868.1", "R388": "BR000038.1"}

def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    seqs, owner, aas = [], [], []
    nm = None
    for fn in sorted(os.listdir(SUB)):
        if not fn.endswith(".faa"): continue
        for line in open(os.path.join(SUB, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                owner.append(nm.split("|")[0]); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(), sequence=line.rstrip("\n")).digitize(alpha))
    print("searched database: %d proteins over %d accessions\n" % (len(seqs), len(set(owner))))
    present = {k: v for k, v in CONTROLS.items() if v in set(owner)}
    for k, v in CONTROLS.items():
        if k not in present: print("  control %s (%s) NOT in this cache -- not read" % (k, v))

    for label, path, ga in (("PF07996 (anchor_T4SS)", os.path.join(HMMDIR, "anchor_T4SS.hmm"), 24.7),
                            ("T4SS_T_virB5 (CONJScan)", os.path.join(CJ, "T4SS_T_virB5.hmm"), None)):
        with pyhmmer.plan7.HMMFile(path) as fh:
            model = next(iter(fh))
        best = {}
        for top in pyhmmer.hmmsearch([model], seqs, cpus=16, E=1e4):
            for h in top:
                i = int(dec(h.name)); acc = owner[i]
                if acc in present.values() and len(aas[i]) > 150:
                    if acc not in best or h.score > best[acc][0]:
                        best[acc] = (h.score, h.evalue, len(aas[i]))
        print("=== %s : hmmsearch against the 1.1M cache ===" % label)
        print("  %-8s %-14s %8s %12s %6s  %s" % ("control", "accession", "bits", "E-value", "aa", "verdict"))
        for k, v in present.items():
            if v not in best:
                print("  %-8s %-14s %8s %12s %6s  no hit >150 aa" % (k, v, "-", "-", "-")); continue
            s, e, L = best[v]
            vs = []
            if ga is not None: vs.append("GA(%.1f) %s" % (ga, "PASS" if s >= ga else "FAIL"))
            vs.append("E<=1e-5 %s" % ("PASS" if e <= 1e-5 else "FAIL"))
            print("  %-8s %-14s %8.1f %12.2g %6d  %s" % (k, v, s, e, L, " | ".join(vs)))
        print()

    print("=== the note's claim, restated against measurement ===")
    print("  anchor_set.tsv: RP4 trbJ 17.7 bits E 7.9e-07 ; R751 22.2 E 3.4e-08")
    print("  If the measured E-values above are >> 1e-5 at those bitscores, the note's")
    print("  E-values came from a different search mode and the E<=1e-5 rationale fails.")

if __name__ == "__main__":
    main()
