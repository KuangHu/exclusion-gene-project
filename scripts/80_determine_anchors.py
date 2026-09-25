#!/usr/bin/env python3
"""Pipeline v2 §1.1 -- determine the anchor set EMPIRICALLY.

Do not copy an enumerated accession list. Scan the four MPF_T positive controls
against the complete Pfam-A at GA thresholds, and read off which family hits each
known T4SS component. The anchor set obtained this way carries its own
positive-control validation, and it yields VirB5's accession as a byproduct
(blocker #2 in §9).

Also records which components have NO usable Pfam family -- those are the blind
spots for slot definition.

RED LINE (§1.3): exclusion-gene families must never enter the anchor set. Any hit
to TIGR04359 / NF033894 / NF041429 / NF033891 / PF10624 / PF14729 is reported
separately under POST-HOC ONLY and excluded from the anchors. Letting them in
would be label leakage -- the scan would "discover" what defined the set.

RED LINE (§1.4): anchors are keyed by VirB number / Pfam family, never by IncP
gene name. trbA has no VirB homolog and would bias the set toward IncP; trbK is
the target itself.
"""
import csv
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from Bio import SeqIO
from Bio.Seq import Seq

import orf_caller

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GB = os.path.join(PROJ, "data", "seed", "genbank")
PFAM = "/global/scratch/users/kh36969/funcannot_dbs/pfam/Pfam-A.hmm"
OUT = os.path.join(PROJ, "data", "anchors")

CONTROLS = [
    ("RP4",    "BN000925.1",  "IncP-alpha"),
    ("R751",   "NC_001735.4", "IncP-beta"),
    ("pKM101", "U09868.1",    "IncN"),
    ("R388",   "NC_028464.1", "IncW"),
]

# known T4SS components per system, from the curated element records.
# Used ONLY to label Pyrodigal calls for read-off -- never as anchors.
COMPONENTS = {
    "RP4":    ["trbB", "trbC", "trbD", "trbE", "trbF", "trbG", "trbH", "trbI",
               "trbJ", "trbL", "trbM", "trbN", "trbO", "trbP"],
    "R751":   ["trbB", "trbC", "trbD", "trbE", "trbF", "trbG", "trbH", "trbI",
               "trbJ", "trbL", "trbM", "trbN", "trbO", "trbP"],
    "pKM101": ["traL", "traM", "traA", "traB", "traC", "traD", "traN", "traE",
               "traO", "traF", "traG"],
    "R388":   ["trwB", "trwC", "trwD", "trwE", "trwF", "trwG", "trwH", "trwI",
               "trwJ", "trwK", "trwL", "trwM", "trwN"],
}

# VirB equivalence, from the curated element records (§1.4: this is the key,
# not the gene name)
VIRB = {
    ("RP4", "trbB"): "VirB11", ("RP4", "trbC"): "VirB2", ("RP4", "trbD"): "VirB3",
    ("RP4", "trbE"): "VirB4", ("RP4", "trbF"): "VirB8", ("RP4", "trbG"): "VirB9",
    ("RP4", "trbJ"): "VirB5", ("RP4", "trbL"): "VirB6", ("RP4", "trbN"): "VirB1",
    ("R751", "trbB"): "VirB11", ("R751", "trbC"): "VirB2", ("R751", "trbD"): "VirB3",
    ("R751", "trbE"): "VirB4", ("R751", "trbF"): "VirB8", ("R751", "trbG"): "VirB9",
    ("R751", "trbJ"): "VirB5", ("R751", "trbL"): "VirB6", ("R751", "trbN"): "VirB1",
    ("pKM101", "traL"): "VirB1", ("pKM101", "traM"): "VirB2", ("pKM101", "traA"): "VirB3",
    ("pKM101", "traB"): "VirB4", ("pKM101", "traC"): "VirB5", ("pKM101", "traD"): "VirB6",
    ("pKM101", "traN"): "VirB7", ("pKM101", "traE"): "VirB8", ("pKM101", "traO"): "VirB9",
    ("pKM101", "traF"): "VirB10", ("pKM101", "traG"): "VirB11",
}

EXCLUSION_FAMILIES = {"PF10624", "PF14729", "TIGR04359", "NF033894",
                      "NF041429", "NF033891"}


def annotated_genes(rec):
    """gene name -> (start, end, strand) from the deposited annotation."""
    out = {}
    for f in rec.features:
        if f.type != "CDS":
            continue
        g = ((f.qualifiers.get("gene") or [""])[0]).strip()
        if g:
            out.setdefault(g, (int(f.location.start) + 1, int(f.location.end),
                               1 if f.location.strand == 1 else -1))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    prots = []          # (id, aa)
    meta = {}           # id -> dict
    for eid, acc, inc in CONTROLS:
        path = os.path.join(GB, "%s__%s.gb" % (eid, acc))
        rec = next(SeqIO.parse(path, "genbank"))
        seq = str(rec.seq).upper()
        ann = annotated_genes(rec)
        calls = orf_caller.call(seq)
        for i, c in enumerate(calls):
            # label this Pyrodigal call with the annotated gene it matches
            label = ""
            for g, (s, e, st) in ann.items():
                if c["start"] == s and c["end"] == e and c["strand"] == st:
                    label = g
                    break
            pid = "%s|%d" % (eid, i)
            prots.append((pid, c["aa"]))
            meta[pid] = {"element": eid, "inc": inc, "gene": label,
                         "virb": VIRB.get((eid, label), ""),
                         "start": c["start"], "end": c["end"],
                         "strand": c["strand"], "aa_len": c["aa_len"]}
        print("  %-8s %-14s %3d Pyrodigal calls, %d matched to an annotated gene"
              % (eid, acc, len(calls),
                 sum(1 for p, _ in prots if p.startswith(eid + "|") and meta[p]["gene"])))

    print("\nscanning %d proteins against the complete Pfam-A at GA thresholds..."
          % len(prots))
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "p.faa")
        with open(fa, "w") as fh:
            for i, s in prots:
                fh.write(">%s\n%s\n" % (i, s))
        dom = os.path.join(OUT, "controls_vs_pfamA.domtbl")
        r = subprocess.run(["hmmscan", "--cut_ga", "--cpu", "8",
                            "--domtblout", dom, PFAM, fa],
                           capture_output=True, text=True)
        if not os.path.exists(dom):
            print("hmmscan failed:\n" + (r.stderr or r.stdout)[-800:])
            return 1

    hits = {}
    for line in open(dom):
        if line.startswith("#"):
            continue
        f = line.split()
        fam_name, fam_acc, pid, score = f[0], f[1], f[3], float(f[13])
        acc_short = fam_acc.split(".")[0]
        cur = hits.get(pid)
        if cur is None or score > cur[2]:
            hits[pid] = (fam_name, acc_short, score)

    rows = []
    for pid, m in meta.items():
        fam = hits.get(pid, ("", "", 0.0))
        rows.append({**m, "protein_id": pid, "pfam_name": fam[0],
                     "pfam_acc": fam[1], "score": round(fam[2], 1)})
    with open(os.path.join(OUT, "control_calls_labelled.tsv"), "w", newline="") as fh:
        cols = ["element", "inc", "gene", "virb", "start", "end", "strand",
                "aa_len", "pfam_name", "pfam_acc", "score", "protein_id"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["element"], r["start"])))

    # ---- read off the anchor set, keyed by VirB number
    print("\n" + "=" * 78)
    print("ANCHOR SET, read off the positive controls (keyed by VirB number)")
    print("=" * 78)
    by_virb = {}
    for r in rows:
        if not r["virb"]:
            continue
        by_virb.setdefault(r["virb"], []).append(r)
    order = ["VirB1", "VirB2", "VirB3", "VirB4", "VirB5", "VirB6", "VirB7",
             "VirB8", "VirB9", "VirB10", "VirB11"]
    anchors = {}
    print("  %-8s %-34s %s" % ("VirB", "Pfam family (acc)", "observed in"))
    print("  " + "-" * 72)
    for v in order:
        rs = by_virb.get(v, [])
        fams = {}
        for r in rs:
            if r["pfam_acc"]:
                fams.setdefault((r["pfam_name"], r["pfam_acc"]), []).append(
                    "%s/%s" % (r["element"], r["gene"]))
        if not fams:
            seen = ", ".join("%s/%s" % (r["element"], r["gene"]) for r in rs)
            print("  %-8s %-34s %s" % (v, "*** NO Pfam HIT (blind spot) ***", seen or "-"))
            continue
        for (nm, ac), where in sorted(fams.items(), key=lambda x: -len(x[1])):
            leak = ac in EXCLUSION_FAMILIES
            tag = "  <-- EXCLUSION FAMILY, POST-HOC ONLY" if leak else ""
            print("  %-8s %-34s %s%s" % (v, "%s (%s)" % (nm, ac),
                                         ", ".join(where), tag))
            if not leak:
                anchors.setdefault(v, set()).add((nm, ac))

    print("\n" + "=" * 78)
    print("BLOCKER #2 -- VirB5 / TrbJ accession")
    print("=" * 78)
    if "VirB5" in anchors:
        for nm, ac in anchors["VirB5"]:
            print("  RESOLVED: VirB5 = %s (%s)" % (nm, ac))
    else:
        print("  UNRESOLVED -- no Pfam family hits VirB5 in these controls.")
        for r in by_virb.get("VirB5", []):
            print("     %s/%s  %d aa  best Pfam: %s"
                  % (r["element"], r["gene"], r["aa_len"], r["pfam_acc"] or "none"))

    # writes the RAW read-off only. data/anchors/anchor_set.tsv is CURATED and
    # hand-maintained -- this script once wrote there and clobbered it, losing the
    # threshold column, the VirB5 rationale and the VirB7 row.
    with open(os.path.join(OUT, "anchor_set_RAW_readoff.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["virb", "pfam_name", "pfam_acc", "role"])
        for v in order:
            for nm, ac in sorted(anchors.get(v, [])):
                role = ("entry_criterion" if v == "VirB4" else
                        "slot_downstream_boundary" if v == "VirB6" else
                        "slot_empty_test" if v == "VirB5" else "annotation")
                w.writerow([v, nm, ac, role])
    print("\nwrote %s" % os.path.relpath(os.path.join(OUT, "anchor_set_RAW_readoff.tsv"), PROJ))
    print("wrote %s" % os.path.relpath(os.path.join(OUT, "control_calls_labelled.tsv"), PROJ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
