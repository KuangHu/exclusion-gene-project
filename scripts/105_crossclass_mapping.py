#!/usr/bin/env python3
"""Map anchor positions across MPF classes, by homology rather than by memory.

The two empirically-read anchor sets share only 2 Pfam families of 27:

    PF03743  TrbI            MPF_T VirB10   =  MPF_F traB
    PF06986  F_T4SS_TraN     MPF_T TraN     =  MPF_F traN

So a class-agnostic catalogue cannot be assembled from the readouts alone. It needs
a POSITION MAPPING layer, and the honest way to build it is the same way the anchor
sets were built: measure it. Not from recalled textbook correspondences, which is
exactly the kind of claim this project has had to retract repeatedly.

Method: take the protein that each anchor family actually hits on each seed, then
phmmer the MPF_F anchor proteins against the MPF_T anchor proteins. A reciprocal
best hit at E <= 1e-3 is evidence of positional correspondence; anything weaker is
left unresolved rather than filled in.

Evidence classes are kept distinct in the output, because they are not
interchangeable:
    shared_family   the same Pfam hits both classes  (strongest)
    phmmer_rbh      reciprocal best hit between seed proteins
    unresolved      no homology found -- the column stays EMPTY, not guessed
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SEEDS = {
    "MPF_T": [("RP4", "RP4__BN000925.1.gb"), ("R751", "R751__NC_001735.4.gb"),
              ("pKM101", "pKM101__U09868.1.gb"), ("R388", "R388__NC_028464.1.gb")],
    "MPF_F": [("F", "F__AP001918.1.gb"), ("R100", "R100__AP000342.1.gb"),
              ("IncC", "pVCR94__CP033514.1.gb"), ("SXT", "SXT__KJ817376.1.gb"),
              ("R27", "R27__AF250878.1.gb")],
}


def read_anchor_set(path, acc_col="pfam_acc", name_col="pfam_name", lab_col=None):
    lines = [l for l in open(path) if not l.startswith("#") and l.strip()]
    out = []
    for r in csv.DictReader(lines, delimiter="\t"):
        out.append((r[acc_col], r[name_col], r.get(lab_col or "", "") or r.get("virb", "")))
    return out


def hmm_for(name):
    for cand in ("anchor_%s.hmm" % name, "mpff_%s.hmm" % name, "%s.hmm" % name):
        p = os.path.join(HMM, cand)
        if os.path.exists(p):
            return p
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=8)
    a = ap.parse_args()
    import pyhmmer
    from Bio import SeqIO
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    T = read_anchor_set(os.path.join(PROJ, "data", "anchors", "anchor_set.tsv"), lab_col="virb")
    F = read_anchor_set(os.path.join(PROJ, "data", "anchors", "anchor_set_MPF_F.tsv"), lab_col="role")
    print("MPF_T anchors %d, MPF_F anchors %d" % (len(T), len(F)))

    # proteins each anchor family actually hits, per class
    def harvest(cls, anchors):
        seqs, meta = [], []
        for sname, fn in SEEDS[cls]:
            p = os.path.join(PROJ, "data", "seed", "genbank", fn)
            if not os.path.exists(p):
                print("  missing seed %s" % fn); continue
            rec = next(SeqIO.parse(p, "genbank"))
            for f in rec.features:
                if f.type != "CDS":
                    continue
                tr = (f.qualifiers.get("translation") or [None])[0]
                if not tr:
                    continue
                gene = ((f.qualifiers.get("gene") or [""])[0] or
                        (f.qualifiers.get("locus_tag") or [""])[0])
                meta.append((sname, gene, len(tr)))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(), sequence=tr).digitize(alpha))
        got = {}
        for accn, name, lab in anchors:
            p = hmm_for(name)
            if not p:
                print("  no HMM for %s (%s)" % (name, accn)); continue
            with pyhmmer.plan7.HMMFile(p) as fh:
                model = next(iter(fh))
            best = {}
            for top in pyhmmer.hmmsearch([model], seqs, cpus=a.cpus,
                                         bit_cutoffs="gathering"):
                for h in top:
                    i = int(dec(h.name))
                    s = meta[i][0]
                    if s not in best or h.score > best[s][1]:
                        best[s] = (i, h.score)
            if best:
                got[(accn, name, lab)] = {s: (meta[i], seqs[i]) for s, (i, _) in best.items()}
        return got
    print("\nharvesting MPF_T anchor proteins ...", flush=True)
    tp = harvest("MPF_T", T)
    print("harvesting MPF_F anchor proteins ...", flush=True)
    fp = harvest("MPF_F", F)
    print("  MPF_T families with proteins: %d ; MPF_F: %d" % (len(tp), len(fp)))

    # phmmer MPF_F anchor proteins vs MPF_T anchor proteins
    tlist, tkey = [], []
    for k, d in tp.items():
        for s, (m, sq) in d.items():
            tkey.append((k, s, m)); tlist.append(sq)
    flist, fkey = [], []
    for k, d in fp.items():
        for s, (m, sq) in d.items():
            fkey.append((k, s, m)); flist.append(sq)
    redig = lambda L, t: [pyhmmer.easel.TextSequence(
        name=("%s%d" % (t, i)).encode(), sequence=x.textize().sequence).digitize(alpha)
        for i, x in enumerate(L)]
    print("\nphmmer %d MPF_F anchor proteins vs %d MPF_T anchor proteins ..."
          % (len(flist), len(tlist)), flush=True)
    pair = collections.defaultdict(list)
    for top in pyhmmer.phmmer(redig(flist, "f"), redig(tlist, "t"), cpus=a.cpus, E=1.0):
        qi = int(dec(top.query.name)[1:])
        for h in top:
            ti = int(dec(h.name)[1:])
            pair[(fkey[qi][0], tkey[ti][0])].append((h.evalue, h.score,
                                                     fkey[qi][1], tkey[ti][1]))
    # shared families are correspondences BY IDENTITY -- they must be entered
    # directly, not left to phmmer. The first version missed both of them.
    rows0, seen0 = [], set()
    tacc = {k[0]: k for k in tp}
    for fk in fp:
        if fk[0] in tacc:
            tk = tacc[fk[0]]
            rows0.append({"virb_position": tk[2] or tk[1], "mpf_T_family": tk[0],
                          "mpf_T_name": tk[1], "mpf_F_family": fk[0],
                          "mpf_F_name": fk[1], "evidence": "shared_family",
                          "best_evalue": "identity", "bit": "", "seed_pair": "-"})
            seen0.add(tk[0])
    print("\n=== SHARED FAMILIES (correspondence by identity) ===")
    for r in rows0:
        print("  %-9s %-20s  MPF_T position %s" % (r["mpf_T_family"], r["mpf_T_name"],
                                                   r["virb_position"]))

    # RECIPROCAL best hit. One-directional best hits let a generic P-loop NTPase
    # cross-hit through (PF10412 -> PF00437 VirB11 at bit 11.0, while the same
    # query's real best is PF02534 VirD4 at bit 164.9).
    fbest, tbest = {}, {}
    for (fk, tk), v in pair.items():
        e = min(x[0] for x in v)
        if fk not in fbest or e < fbest[fk][0]:
            fbest[fk] = (e, tk)
        if tk not in tbest or e < tbest[tk][0]:
            tbest[tk] = (e, fk)

    print("\n=== CROSS-CLASS CORRESPONDENCE (reciprocal best hit, E <= 1e-3) ===")
    print("%-9s %-18s -> %-9s %-18s %10s %6s %s"
          % ("MPF_F", "family", "MPF_T", "position", "best E", "bit", "seed pair"))
    rows = list(rows0)
    seen = set(seen0)
    for (fk, tk), v in sorted(pair.items(), key=lambda kv: min(x[0] for x in kv[1])):
        e, sc, fs, ts = min(v, key=lambda x: x[0])
        if e > 1e-3:
            continue
        if fbest.get(fk, (None, None))[1] != tk or tbest.get(tk, (None, None))[1] != fk:
            continue                          # not reciprocal -- likely a fold cross-hit
        if tk[0] in seen0:
            continue
        print("%-9s %-18s -> %-9s %-18s %10.1e %6.1f %s/%s"
              % (fk[0], fk[1][:18], tk[0], (tk[2] or tk[1])[:18], e, sc, fs, ts))
        rows.append({"virb_position": tk[2] or tk[1], "mpf_T_family": tk[0],
                     "mpf_T_name": tk[1], "mpf_F_family": fk[0], "mpf_F_name": fk[1],
                     "evidence": "shared_family" if fk[0] == tk[0] else "phmmer_hit",
                     "best_evalue": "%.1e" % e, "bit": round(sc, 1),
                     "seed_pair": "%s/%s" % (fs, ts)})
        seen.add(tk[0])
    print("\n=== MPF_T positions with NO MPF_F correspondence found ===")
    for accn, name, lab in T:
        if accn not in seen:
            print("  %-9s %-20s %-12s  -> unresolved (column left EMPTY, not guessed)"
                  % (accn, name, lab))
            rows.append({"virb_position": lab, "mpf_T_family": accn, "mpf_T_name": name,
                         "mpf_F_family": "", "mpf_F_name": "", "evidence": "unresolved",
                         "best_evalue": "", "bit": "", "seed_pair": ""})
    dest = os.path.join(PROJ, "data", "anchors", "crossclass_position_map.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["virb_position", "mpf_T_family", "mpf_T_name",
                                           "mpf_F_family", "mpf_F_name", "evidence",
                                           "best_evalue", "bit", "seed_pair"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
