#!/usr/bin/env python3
"""Type the four MPF classes in the ICE/ciMGE database.

Applies the FROZEN entry criteria from data/release/v1.1 to a new dataset:

    MPF_T   PF03135  anchor_CagE_TrbE_VirB   (VirB4)
    MPF_F   PF11130  mpff_TraC_F_IV          (TraC, F-specific)
    MPF_I   T4SS_I_traU
    MPF_FA  T4SS_virb4 AND NOT(T/F/I) AND >=2 FA/FATA profiles

Genes are called with the pipeline's frozen Pyrodigal settings (ANON=True,
MIN_GENE=90, CLOSED=True, TRANSL_TABLE=11, enforced by A1) so results are
directly comparable to the plasmid release. Search is hmmsearch at GATHERING
thresholds only (A13 -- no E-value cutoffs).

THE CONTROLS COME FREE WITH THIS DATASET:

  ICE   2,634  self-transmissible   -> SHOULD carry MPF
  IME   1,490  mobilizable          -> relaxase but NO MPF (this is what makes
                                        it an IME rather than an ICE)
  AICE     86  FtsK/SpoIIIE         -> NO T4SS at all
  ICEBs1    1  NC_000964:529362..549932 -> must type MPF_FA (our own seed)

Relaxase is present in ~95% of BOTH ICE and IME, so it cannot separate them.
MPF should. If ICE >> IME does not hold, the typing is wrong and nothing else
here is read.
"""
import collections, csv, gzip, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import orf_caller

ICE = "/global/scratch/users/kh36969/ice_db"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
OUT = "/global/scratch/users/kh36969/exclusion_gene/mge_search"
EEX = ["TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"]
NAME = {"TIGR04359": "TrbK_RP4", "NF033894": "Eex_IncN", "NF041429": "EexR",
        "NF033891": "ExcA", "PF10624": "TraS", "PF14729": "DUF4467"}


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    os.makedirs(OUT, exist_ok=True)

    meta = {}
    for r in csv.DictReader(open(os.path.join(ICE, "subset", "ice_elements.tsv")),
                            delimiter="\t"):
        meta[r["element_id"]] = r
    print("ICE elements in table: %d" % len(meta), flush=True)

    seqs, owner, aas, ncds = [], [], [], collections.Counter()
    n = 0
    with gzip.open(os.path.join(ICE, "seq", "ice.shard000.fna.gz"), "rt") as fh:
        name, buf = None, []
        def flush():
            nonlocal name, buf
            if name and buf:
                s = "".join(buf).upper()
                for g in orf_caller.call(s):
                    p = g["aa"] if isinstance(g, dict) else g
                    if not p or p == "*": continue
                    p = p.rstrip("*")
                    if len(p) < 30: continue
                    owner.append(name); aas.append(p); ncds[name] += 1
            name, buf = None, []
        for line in fh:
            if line[0] == ">":
                flush(); name = line[1:].split()[0]; buf = []; n += 1
            else:
                buf.append(line.strip())
        flush()
    print("sequences %d ; proteins %d" % (n, len(aas)), flush=True)
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
    print("\nentry criteria: MPF_T %d | MPF_F %d | MPF_I %d | T4SS_virb4 %d"
          % (len(T), len(F), len(I), len(vb4)), flush=True)

    fa = sorted(f[:-4] for f in os.listdir(CJ)
                if (f.startswith("FA_") or f.startswith("FATA_")) and f.endswith(".hmm"))
    nprof = collections.Counter()
    for p in fa:
        for a in scan(os.path.join(CJ, p + ".hmm")): nprof[a] += 1
    prev = T | F | I
    FA = {a for a in vb4 if a not in prev and nprof[a] >= 2}
    print("MPF_FA (virb4 AND NOT(T/F/I) AND >=2 FA/FATA): %d" % len(FA), flush=True)

    eex = {}
    for f in EEX:
        eex[NAME[f]] = scan(os.path.join(HMM, f + ".hmm"))
    print("exclusion families: %s"
          % ", ".join("%s %d" % (k, len(v)) for k, v in eex.items() if v), flush=True)

    rows = []
    for acc in sorted(set(owner)):
        m = meta.get(acc, {})
        cls = [c for c, s in (("MPF_T", T), ("MPF_F", F), ("MPF_I", I), ("MPF_FA", FA))
               if acc in s]
        ex = [k for k, v in eex.items() if acc in v]
        rows.append({"element_id": acc, "element_type": m.get("element_type", "?"),
                     "length": m.get("length", ""), "n_cds_called": ncds[acc],
                     "has_relaxase": m.get("has_relaxase", ""),
                     "relaxase_families": m.get("relaxase_families", ""),
                     "mpf_classes": ";".join(cls) or "none",
                     "n_mpf_classes": len(cls),
                     "n_fa_fata_profiles": nprof.get(acc, 0),
                     "exclusion_family_hit": ";".join(ex),
                     "phylum": m.get("phylum", ""), "genus": m.get("genus", "")})
    dest = os.path.join(PROJ, "data", "anchors", "mpf_in_ice.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    print("\n=== THE CONTROL: MPF by element_type ===")
    print("  %-10s %7s %9s %9s %9s" % ("type", "n", "any MPF", "%", "relaxase%"))
    for t in ("ICE", "IME", "AICE", "NA", "Conj_reg"):
        sel = [r for r in rows if r["element_type"] == t]
        if not sel: continue
        k = sum(1 for r in sel if r["mpf_classes"] != "none")
        rx = sum(1 for r in sel if r["has_relaxase"] == "1")
        print("  %-10s %7d %9d %8.1f%% %8.1f%%"
              % (t, len(sel), k, 100.0*k/len(sel), 100.0*rx/len(sel)))
    print("\n  PREDICTION: ICE >> IME (self-transmissible vs mobilizable).")
    print("  Relaxase is ~95%% in BOTH and cannot separate them.")

    print("\n=== MPF class distribution ===")
    cc = collections.Counter(r["mpf_classes"] for r in rows)
    for k, v in cc.most_common(10):
        print("  %-24s %6d (%.1f%%)" % (k, v, 100.0*v/len(rows)))

    print("\n=== POSITIVE CONTROL: ICEBs1 ===")
    for r in rows:
        if r["element_id"].startswith("NC_000964_529362"):
            print("  %s  len %s  MPF=%s  FA/FATA profiles=%d  eex=%s"
                  % (r["element_id"], r["length"], r["mpf_classes"],
                     r["n_fa_fata_profiles"], r["exclusion_family_hit"] or "none"))
            print("  -> expected MPF_FA: %s"
                  % ("PASS" if "MPF_FA" in r["mpf_classes"] else "FAIL"))
    print("\n=== exclusion genes found in ICEs ===")
    ec = collections.Counter(x for r in rows for x in r["exclusion_family_hit"].split(";") if x)
    for k, v in ec.most_common(): print("  %-12s %5d" % (k, v))
    print("\nwrote %s (%d rows)" % (dest, len(rows)))


if __name__ == "__main__":
    main()
