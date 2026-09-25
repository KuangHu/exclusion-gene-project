#!/usr/bin/env python3
"""Unified PLSDB annotation for the two studied classes: MPF_T (trb) and MPF_F (tra).

One row per plasmid, both column blocks present. A plasmid admitted by only one
entry criterion carries `not_admitted` in the other block -- it is not dropped, and
the two blocks are never merged into a single "class" call.

Row count 18,381 = 7,436 (MPF_T) + 11,120 (MPF_F) - 175 admitted by both.

**dual_system is not bookkeeping.** 144 of those 175 are tier A in BOTH
catalogues: two complete, independently scored T4SS on one plasmid. MOB-suite
calls 143 of them MPF_F and 31 MPF_T -- its single-valued mpf_type cannot
represent a dual system, so it reports one and misses the other. That is exactly
why the classes are kept as separate blocks rather than resolved to one label.

Per anchor family: gene_id, start, end, strand, bitscore. Coordinates were already
complete in both catalogues (zero malformed); bitscore is the only field being
added -- from anchor_matrix.tsv for MPF_T, by rescan for MPF_F.

family_call_mode is GA throughout. An E-value exception for NF041429 was
considered and is NOT needed: its GA of 200 on a 142-residue model looked
prohibitive, but SXT EexR1 scores 274.7 (E 2.5e-86). The model fires. EexR's
absence from MPF_F slots is therefore not a threshold artefact, and the most
likely explanation -- SXT/R391-family elements being chromosomal ICEs and
under-represented in a plasmid database -- is a hypothesis, not a measurement.

NOT an operon annotation: HMM hits and gene positions, no promoter, terminator or
co-transcription evidence.
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
MPFF_FAMS = ["TraU", "TraH", "TraF", "TraC_F_IV", "TraG_N", "TrbI", "P-loop_TraG",
             "TraE", "TrbC_Ftype", "F_T4SS_TraN", "TrwB_AAD_bind", "TraV", "TraG-D_C"]


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    csv.field_size_limit(10 ** 8)
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    T = {r["accession"]: r for r in csv.DictReader(
        open(os.path.join(PROJ, "data", "release", "v1", "catalogue_MPF_T_v1.tsv")), delimiter="\t")}
    F = {r["accession"]: r for r in csv.DictReader(
        open(os.path.join(PROJ, "data", "release", "v1", "catalogue_MPF_F_v1.tsv")), delimiter="\t")}
    AM = {r["accession"]: r for r in csv.DictReader(
        open(os.path.join(PROJ, "data", "anchors", "anchor_matrix.tsv")), delimiter="\t")}
    TANCH = [c for c in next(iter(T.values()))
             if c.startswith(("VirB", "VirD", "TraN")) and not c.endswith("_score")]
    print("MPF_T %d rows (%d anchors), MPF_F %d rows (%d anchors)"
          % (len(T), len(TANCH), len(F), len(MPFF_FAMS)), flush=True)

    # --- bitscores for MPF_F (the only missing field) -----------------------
    seqs, owner, idx, aas = [], [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                owner.append(f[0]); idx.append(int(f[1])); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    fscore = collections.defaultdict(dict)
    for fam in MPFF_FAMS:
        p = None
        for c in ("mpff_%s.hmm" % fam, "anchor_%s.hmm" % fam, "%s.hmm" % fam):
            if os.path.exists(os.path.join(HMM, c)):
                p = os.path.join(HMM, c); break
        if not p:
            continue
        with pyhmmer.plan7.HMMFile(p) as fh:
            model = next(iter(fh))
        for top in pyhmmer.hmmsearch([model], seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                cur = fscore[fam].get(owner[i])
                if cur is None or h.score > cur[1]:
                    fscore[fam][owner[i]] = (idx[i], h.score)
        print("  bitscores %-16s %6d" % (fam, len(fscore[fam])), flush=True)

    def tup(cell, score):
        """ga:idx:start..end:strand  ->  gene_id;start;end;strand;bitscore"""
        if not cell or cell == "none":
            return ""
        _, gi, coords, strand = cell.split(":")
        s, e = coords.split("..")
        return "%s;%s;%s;%s;%s" % (gi, s, e, strand, score if score != "" else "NA")

    rows = []
    for acc in sorted(set(T) | set(F)):
        t, f = T.get(acc), F.get(acc)
        r = {"accession": acc,
             "plsdb_mpf_type": (t or f).get("mpf_class") or (t or f).get("mob_mpf_type", "-"),
             "replicons": (t or f).get("replicons", "-"),
             "n_cds": (t or f).get("n_cds", ""),
             "dual_system": int(bool(t) and bool(f)),
             "dual_system_both_tierA": int(bool(t) and bool(f)
                                           and t["tier"] == "A" and f["tier"] == "A"),
             "family_call_mode": "GA"}
        # --- MPF_T block
        r["mpft_admitted"] = int(bool(t))
        if t:
            am = AM.get(acc, {})
            for c in TANCH:
                r["mpft_" + c] = tup(t.get(c), am.get(c + "_score", ""))
            r.update({"mpft_tier": t["tier"],
                      "mpft_core_completeness_ga": t["core_completeness_ga"],
                      "mpft_architecture": t["architecture"],
                      "mpft_slot_ready__virb5_virb6": t.get("slot_ready__virb5_virb6", ""),
                      "mpft_fused_virB3_virB4": t.get("fused_virB3_virB4", ""),
                      "mpft_virB7_status": t.get("virB7_status", ""),
                      "mpft_n_gap_genes_le10kb": t.get("n_gap_genes_le10kb", ""),
                      "mpft_carries_novel84": t.get("carries_novel84_slot_lipoprotein", "")})
        else:
            for c in TANCH:
                r["mpft_" + c] = ""
            for k in ("mpft_tier", "mpft_core_completeness_ga", "mpft_architecture",
                      "mpft_slot_ready__virb5_virb6", "mpft_fused_virB3_virB4",
                      "mpft_virB7_status", "mpft_n_gap_genes_le10kb", "mpft_carries_novel84"):
                r[k] = "not_admitted"
        # --- MPF_F block
        r["mpff_admitted"] = int(bool(f))
        if f:
            for c in MPFF_FAMS:
                sc = fscore.get(c, {}).get(acc)
                r["mpff_" + c] = tup(f.get(c), round(sc[1], 1) if sc else "")
            r.update({"mpff_tier": f["tier"],
                      "mpff_core_completeness_fspecific": f["core_completeness_fspecific"],
                      "mpff_layout": f["layout"],
                      "mpff_anchor_span_bp": f["anchor_span_bp"],
                      "mpff_slot_ready__traG": f["slot_ready__traG"],
                      "mpff_traN_present": f["traN_present"],
                      "mpff_slot_occupant": ";".join([f["slot_coords"], str(f["slot_aa_len"]),
                                                      f["slot_eex_family"] or "unnamed",
                                                      str(f["slot_lipobox"])]) if f["slot_ready__traG"] == "1" else ""})
        else:
            for c in MPFF_FAMS:
                r["mpff_" + c] = ""
            for k in ("mpff_tier", "mpff_core_completeness_fspecific", "mpff_layout",
                      "mpff_anchor_span_bp", "mpff_slot_ready__traG",
                      "mpff_traN_present", "mpff_slot_occupant"):
                r[k] = "not_admitted"
        rows.append(r)

    cols = list(rows[0].keys())
    dest = os.path.join(PROJ, "data", "annotation_MPF_T_F_regen.tsv")  # A9: see 106
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           restval="")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d rows x %d cols)" % (dest, len(rows), len(cols)))
    print("  MPF_T only %d, MPF_F only %d, both %d"
          % (sum(1 for r in rows if r["mpft_admitted"] and not r["mpff_admitted"]),
             sum(1 for r in rows if r["mpff_admitted"] and not r["mpft_admitted"]),
             sum(1 for r in rows if r["dual_system"])))
    print("  dual_system tier A in both: %d" % sum(r["dual_system_both_tierA"] for r in rows))


if __name__ == "__main__":
    main()
