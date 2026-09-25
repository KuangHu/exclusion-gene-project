#!/usr/bin/env python3
"""Gene architecture of all eight MPF classes, read off the prototypes.

Everything so far has counted profiles PER ELEMENT. That says which components a
class has, not how they are arranged. This walks each class's reference element in
coordinate order and prints what sits where, so the architecture is read off the
sequence rather than recalled from a paper.

One element per class -- the published prototype. This is an ARCHITECTURE readout,
not a population statistic: n=1 per class by construction, and nothing here may be
quoted as "MPF_x operons look like this" without the population behind it. What it
IS good for is locating the exclusion gene relative to the machinery, which is the
thing the slot work depends on and which no count can show.

Exclusion families are scanned LAST and POST-HOC, exactly as everywhere else in
this project. They are never anchors and never enter a class definition:

    TIGR04359  TrbK      RP4          lipoprotein, C-term 8-aa truncation kills it
    PF10624    TraS      F            blocks transfer after mating-pair formation
    NF033891   ExcA      R64          donor TraY is its target (Sakuma 2013)
    NF033894   Eex_IncN  pKM101       single lipid-attachment-motif protein
    NF041429   EexR/S    R391/SXT     Marrero & Waldor
    PF14729    DUF4467   ICEBs1       ConG E288K resistance; a DUF on the LITERATURE

Genes with no profile hit are collapsed into gaps, since the point is the spacing
between recognised components, not a full annotation.
"""
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import orf_caller

CJ = os.path.join(PROJ, "data", "CONJScan", "profiles")
DEFS = os.path.join(PROJ, "data", "CONJScan", "definitions", "Plasmids")
HMM = os.path.join(PROJ, "data", "hmm")
GB = os.path.join(PROJ, "data", "seed", "genbank")
PILOT = os.path.join(PROJ, "data", "seed", "pilot_mpf")

# class -> (label, reference element, path)
SEEDS = [
    ("T",    "RP4",        os.path.join(GB, "RP4__BN000925.1.gb")),
    ("F",    "F",          os.path.join(GB, "F__AP001918.1.gb")),
    ("I",    "R64",        os.path.join(GB, "R64__AP005147.1.gb")),
    ("FA",   "ICEBs1",     os.path.join(GB, "ICEBs1__CP171645.1.gb")),
    ("FATA", "pCF10",      os.path.join(GB, "pCF10__AY855841.2.gb")),
    ("G",    "ICEHin1056", os.path.join(PILOT, "AJ627386.gb")),
    ("B",    "CTnDOT",     os.path.join(PILOT, "AF289050.gb")),
    ("C",    "pCC7120a",   os.path.join(PILOT, "BA000020.gb")),
]
GENERIC = ["T4SS_virb4", "T4SS_I_traU", "T4SS_t4cp1", "T4SS_t4cp2", "T4SS_tcpA"]
MOB = ["T4SS_MOB" + x for x in ("B", "C", "F", "H", "P1", "P2", "P3", "Q", "T", "V")]
EEX = [("TIGR04359", "TrbK"), ("PF10624", "TraS"), ("NF033891", "ExcA"),
       ("NF033894", "Eex_IncN"), ("NF041429", "EexR/S"), ("PF14729", "DUF4467")]
MAXGAP = 3          # genes; wider than this is printed as a break


def class_profiles(cls):
    import xml.etree.ElementTree as ET
    names = {g.get("name") for g in
             ET.parse(os.path.join(DEFS, "T4SS_type%s.xml" % cls)).getroot().iter("gene")}
    return sorted(names)


def proteins(path):
    """(index, start, end, strand, aa) in coordinate order."""
    from Bio import SeqIO
    out = []
    recs = list(SeqIO.parse(path, "genbank"))
    got = False
    for rec in recs:
        for f in rec.features:
            if f.type != "CDS":
                continue
            p = f.qualifiers.get("translation", [None])[0]
            if not p:
                continue
            got = True
            out.append((int(f.location.start) + 1, int(f.location.end),
                        1 if f.location.strand == 1 else -1, p))
    if not got:                      # no annotation -- call them
        orf_caller.assert_anon_only(orf_caller.ANON)
        for rec in recs:
            for g in orf_caller.call(str(rec.seq).upper()):
                out.append((g["start"], g["end"], g["strand"], g["aa"]))
    out.sort(key=lambda x: x[0])
    return [(i,) + t for i, t in enumerate(out)]


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    rows = []
    for cls, label, path in SEEDS:
        if not os.path.exists(path):
            print("MISSING seed for MPF_%s: %s" % (cls, path)); return 1
        genes = proteins(path)
        seqs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=g[4])
                .digitize(alpha) for i, g in enumerate(genes)]

        best = {}        # gene index -> (bitscore, label)
        def run(name, disp, path_):
            if not os.path.exists(path_):
                return False
            with pyhmmer.plan7.HMMFile(path_) as fh:
                hmm = next(iter(fh))
            for top in pyhmmer.hmmsearch([hmm], seqs, cpus=4, bit_cutoffs="gathering"):
                for h in top:
                    i = int(dec(h.name))
                    if i not in best or h.score > best[i][0]:
                        best[i] = (h.score, disp)
            return True

        for p in class_profiles(cls) + GENERIC + MOB:
            run(p, p.replace("T4SS_", "").replace(cls + "_", ""),
                os.path.join(CJ, p + ".hmm"))
        # exclusion families LAST, and they overwrite -- we want them visible
        eex_at = {}
        for mdl, nm in EEX:
            hit_before = set(best)
            if not run(mdl, "**" + nm + "**", os.path.join(HMM, mdl + ".hmm")):
                print("  WARNING: %s missing, MPF_%s readout is incomplete" % (mdl, cls))
            for i in set(best) - hit_before:
                eex_at[i] = nm
            for i in hit_before:
                if best[i][1] == "**" + nm + "**":
                    eex_at[i] = nm

        idx = sorted(best)
        print("\n" + "=" * 74)
        print("MPF_%-5s  %-12s  %d CDS, %d recognised" % (cls, label, len(genes), len(idx)))
        print("=" * 74)
        prev = None
        line = []
        for i in idx:
            if prev is not None and i - prev > MAXGAP:
                print("   " + "  ".join(line)); line = []
                print("   ---- %d genes ----" % (i - prev - 1))
            line.append("[%s]" % best[i][1])
            prev = i
        if line:
            print("   " + "  ".join(line))

        if eex_at:
            for i, nm in sorted(eex_at.items()):
                nb = [best[j][1] for j in idx if 0 < abs(j - i) <= 2]
                print("\n   EXCLUSION GENE: %s at gene %d (%d aa), neighbours: %s"
                      % (nm, i, len(genes[i][4]), ", ".join(nb) or "none within 2"))
        else:
            print("\n   EXCLUSION GENE: none of the six families hits this element")

        for i in idx:
            rows.append({"mpf_class": "MPF_" + cls, "element": label,
                         "gene_index": i, "start": genes[i][1], "end": genes[i][2],
                         "strand": genes[i][3], "aa": len(genes[i][4]),
                         "profile": best[i][1].strip("*"),
                         "is_exclusion_family": int(i in eex_at)})

    dest = os.path.join(PROJ, "data", "anchors", "architecture_all8.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
