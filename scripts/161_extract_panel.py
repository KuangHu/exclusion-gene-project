#!/usr/bin/env python3
"""Extract the P1-P8 / N1-N4 test panel from the local GenBank seeds.

The element YAMLs carry gene names and coordinates but no stored sequences
(`sequence: None` throughout, `NOT_RECOVERED` for IncC eexC and pKPC_UVA01).
So the panel is pulled from the GenBank records directly.

MATCHING IS BY GENE/PRODUCT TEXT, which is exactly the failure mode already hit
three times in this project (the P1 regex matching Mycoplasma phage P1, PGAP
dropping R100 traS's /gene, the pED208 sfx name collision). So every hit is
reported with its accession, coordinates, length and product string, and nothing
is used without the length matching the literature value.

Literature lengths to check against:
  TrbK   RP4    69 aa precursor, 47 aa mature
  YddJ   ICEBs1 126 aa (conJ, DUF4467)
  TraS   F      173 aa
  TraT   F      244 aa
"""
import os
import re
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GB = os.path.join(PROJ, "data", "seed", "genbank")
OUT = os.path.join(PROJ, "data", "anchors", "panel.faa")

# id, seed file, regex over gene+product+locus_tag, role, expected aa (None = unknown)
PANEL = [
    ("P1_TrbK_RP4",     "RP4__BN000925.1",      r"\btrbK\b",                "positive", 69),
    ("P2_YddJ_ICEBs1",  "ICEBs1__CP171645.1",   r"\bconJ\b|\byddJ\b|DUF4467", "positive", 126),
    ("P3_Eex_pKM101",   "pKM101__U09868.1",     r"\beex\b",                 "positive", None),
    ("P4_EexS_SXT",     "SXT__AY034138.1",      r"\beex\b|\bs0", "positive", None),
    ("P5_EexR_R391",    "R391__AY090559.1",     r"\beex\b",                 "positive", None),
    ("P6_TraS_F",       "F__AP001918.1",        r"\btraS\b",                "positive", 173),
    ("P7_ExcA_R64",     "R64__AP005147.1",      r"\bexcA\b",                "positive", None),
    ("N1_TraT_F",       "F__AP001918.1",        r"\btraT\b",                "negative", 244),
]


def main():
    from Bio import SeqIO
    rows, seqs = [], []
    for pid, seed, pat, role, exp in PANEL:
        p = os.path.join(GB, seed + ".gb")
        if not os.path.exists(p):
            print("  %-18s SEED MISSING (%s)" % (pid, seed))
            continue
        rec = next(SeqIO.parse(p, "genbank"))
        rx = re.compile(pat, re.I)
        hits = []
        for ft in rec.features:
            if ft.type != "CDS":
                continue
            txt = " ".join(sum((ft.qualifiers.get(k, []) for k in
                                ("gene", "product", "locus_tag", "note")), []))
            if rx.search(txt):
                aa = ft.qualifiers.get("translation", [""])[0]
                if not aa:
                    continue
                g = (ft.qualifiers.get("gene") or ft.qualifiers.get("locus_tag") or ["?"])[0]
                pr = (ft.qualifiers.get("product") or [""])[0]
                hits.append((g, len(aa), int(ft.location.start) + 1,
                             int(ft.location.end), pr, aa))
        if not hits:
            print("  %-18s NO MATCH in %s" % (pid, rec.id))
            continue
        # prefer the hit whose length matches the literature value
        best = None
        if exp:
            exact = [h for h in hits if h[1] == exp]
            if exact:
                best = exact[0]
        if best is None:
            best = hits[0]
        g, n, s, e, pr, aa = best
        flag = ""
        if exp and n != exp:
            flag = "  *** LENGTH MISMATCH: expected %d aa ***" % exp
        print("  %-18s %-14s %-12s %4d aa  %d..%d  %-34s%s"
              % (pid, rec.id, g, n, s, e, pr[:34], flag))
        if len(hits) > 1:
            print("      (%d candidates matched; others: %s)"
                  % (len(hits), ", ".join("%s/%daa" % (h[0], h[1]) for h in hits[1:4])))
        rows.append((pid, role, rec.id, g, n, s, e, pr, exp, flag != ""))
        seqs.append((pid, aa))

    with open(OUT, "w") as fh:
        for pid, aa in seqs:
            fh.write(">%s\n%s\n" % (pid, aa))
    print("\nwrote %s (%d sequences)" % (OUT, len(seqs)))
    miss = [p[0] for p in PANEL if p[0] not in {s[0] for s in seqs}]
    if miss:
        print("  NOT RECOVERED: %s" % ", ".join(miss))
    bad = [r[0] for r in rows if r[9]]
    if bad:
        print("  LENGTH MISMATCH (do not use without checking): %s" % ", ".join(bad))


if __name__ == "__main__":
    main()
