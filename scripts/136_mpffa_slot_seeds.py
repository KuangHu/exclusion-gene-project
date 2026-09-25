#!/usr/bin/env python3
"""MPF_FA §3: the slot, measured on the conG neighbourhood of all five seeds.

On ICEBs1 the slot is conG+3 and is occupied by conJ, which carries PF14729
(DUF4467) at 99.1 bits -- job 25974613, CDS indices conG 18, cwlT 19,
ACGWZ1_02860 20, conJ 21.

Two candidate anchors, because neither covers all five seeds:

    FA_orf15  marks conG  -- fires on 2 of 5 (ICEBs1, pAM373)
    FA_orf14  marks cwlT  -- fires on 3 of 5 (ICEBs1, pLS20, pAD1)
              cwlT is conG+1, so the same slot is cwlT+2

Both are reported per seed with their offsets. Neither is promoted to "the"
anchor here; that needs the class-scale census the slot control does not yet have.

WHAT THIS CANNOT DO: PF14729 fires on ICEBs1 ONLY. No exclusion family fires on
pLS20, pCF10, pAD1 or pAM373. So on four of five seeds the slot occupant can be
located positionally but CANNOT be confirmed as an exclusion gene by any named
family. That is the single-element control weakness, restated at the slot.
"""
import collections, csv, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
GB = os.path.join(PROJ, "data", "seed", "genbank")
SEEDS = [("ICEBs1", "ICEBs1__CP171645.1"), ("pLS20", "pLS20__AB615352.1"),
         ("pCF10", "pCF10__AY855841.2"), ("pAD1", "pAD1__CP046109.1"),
         ("pAM373", "pAM373__AE002565.1")]
EEX = ["PF14729", "PF10624", "NF033891", "NF033894", "NF041429", "TIGR04359"]
ANCHORS = [("FA_orf15", "conG", 3), ("FA_orf14", "cwlT", 2)]
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
import re as _re
LIPO = _re.compile(r"[LVI][ASTVIG][GASN]C")


def maxkd(s, w=19):
    if len(s) < w:
        return round(sum(KD.get(c, 0) for c in s) / max(1, len(s)), 2)
    return round(max(sum(KD.get(c, 0) for c in s[i:i+w]) / w
                     for i in range(len(s)-w+1)), 2)


def main():
    require_compute_node()
    import pyhmmer
    from Bio import SeqIO
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner, meta, aas = [], [], [], []
    genes = collections.defaultdict(list)
    for name, fn in SEEDS:
        rec = next(SeqIO.parse(os.path.join(GB, fn + ".gb"), "genbank"))
        k = 0
        for f in rec.features:
            if f.type != "CDS": continue
            p = f.qualifiers.get("translation", [None])[0]
            if not p: continue
            g = (f.qualifiers.get("gene") or f.qualifiers.get("locus_tag") or ["?"])[0]
            pr = (f.qualifiers.get("product") or [""])[0]
            owner.append(name); aas.append(p)
            meta.append((k, g, int(f.location.start)+1, int(f.location.end),
                         1 if f.location.strand == 1 else -1, pr))
            genes[name].append(len(owner)-1); k += 1
            seqs.append(pyhmmer.easel.TextSequence(
                name=str(len(seqs)-0).encode(), sequence=p).digitize(alpha))
    for i, s in enumerate(seqs): s.name = str(i).encode()

    miss = [p for p in EEX if not os.path.exists(os.path.join(HMM, p + ".hmm"))]
    if miss:
        raise SystemExit("MISSING named control models: %s" % ", ".join(miss))

    def scan(path):
        with pyhmmer.plan7.HMMFile(path) as fh: m = next(iter(fh))
        best = {}
        for top in pyhmmer.hmmsearch([m], seqs, cpus=4, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                if owner[i] not in best or h.score > best[owner[i]][1]:
                    best[owner[i]] = (i, round(h.score, 1))
        return best

    prof = {}
    for f in sorted(os.listdir(CJ)):
        if f.endswith(".hmm") and (f.startswith("FA_") or f.startswith("FATA_")):
            prof[f[:-4]] = scan(os.path.join(CJ, f))
    eex = {p: scan(os.path.join(HMM, p + ".hmm")) for p in EEX}

    # what profile/family marks each protein
    mark = collections.defaultdict(list)
    for p, d in prof.items():
        for acc, (i, sc) in d.items(): mark[i].append("%s:%.0f" % (p, sc))
    for p, d in eex.items():
        for acc, (i, sc) in d.items(): mark[i].append("**%s:%.0f**" % (p, sc))

    rows = []
    for aname, gname, off in ANCHORS:
        print("\n=== anchor %s (marks %s); slot = %s+%d ==="
              % (aname, gname, gname, off))
        d = prof.get(aname, {})
        for sname, _ in SEEDS:
            if sname not in d:
                print("  %-8s anchor ABSENT -- slot not locatable from this anchor" % sname)
                continue
            ai, asc = d[sname]
            k0, g0, s0, e0, st0, _ = meta[ai]
            idxmap = {meta[i][0]: i for i in genes[sname]}
            print("  %-8s anchor %s at CDS %d (%.1f bits)" % (sname, g0, k0, asc))
            for o in range(0, off + 2):
                j = idxmap.get(k0 + o * (1 if st0 == 1 else -1))
                if j is None:
                    print("      %+d  (no CDS)" % o); continue
                k, g, s, e, stx, pr = meta[j]
                tagp = ",".join(mark.get(j, [])) or "-"
                flag = "  <== SLOT" if o == off else ""
                print("      %+d  %-16s %4d aa  %6d..%-6d  %-28s %s%s"
                      % (o, g, len(aas[j]), s, e, pr[:28], tagp, flag))
                if o == off:
                    rows.append({"seed": sname, "anchor": aname, "anchor_gene": g0,
                                 "offset": off, "slot_gene": g, "slot_aa": len(aas[j]),
                                 "slot_start": s, "slot_end": e,
                                 "slot_strand": "+" if stx == 1 else "-",
                                 "slot_product": pr[:60], "slot_maxkd": maxkd(aas[j]),
                                 "slot_lipobox": int(bool(LIPO.search(aas[j][:40]))),
                                 "slot_family": ";".join(m for m in mark.get(j, [])
                                                         if m.startswith("**")) or ""})
    print("\n=== named-family confirmation at the slot ===")
    for r in rows:
        print("  %-8s %-9s %s+%d -> %-16s %4d aa  family %s"
              % (r["seed"], r["anchor"], r["anchor_gene"], r["offset"],
                 r["slot_gene"], r["slot_aa"], r["slot_family"] or "NONE"))
    n = sum(1 for r in rows if r["slot_family"])
    print("\n  %d of %d located slots carry a named exclusion family." % (n, len(rows)))
    print("  On every other seed the slot is POSITIONAL ONLY -- located, not confirmed.")
    dest = os.path.join(PROJ, "data", "anchors", "mpffa_slot_seeds.tsv")
    if rows:
        with open(dest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(rows)
        print("\nwrote %s (%d rows)" % (dest, len(rows)))


if __name__ == "__main__":
    main()
