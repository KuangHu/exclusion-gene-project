#!/usr/bin/env python3
"""Per-plasmid anchor matrix -- the tra gene cluster catalogue's core table.

One row per accession, one column per anchor family, recording HOW the anchor was
found (ga / none for now; phmmer and sixframe columns are filled by later passes)
plus its gene index and coordinates. Deliberately NOT called an operon table: this
is HMM hits and positions, with no promoter, terminator or co-transcription
evidence.

Also answers two questions that decide how the catalogue is tiered:

  * Of the ~10% of plasmids that do NOT have all 7 core anchors, WHICH anchor is
    missing? Concentrated in one or two positions means it is still a model
    problem; spread out means genuinely incomplete systems.

  * VirB1 only fell from 19.5% to 13.7% under phmmer recall, unlike VirB2
    (20.9 -> 3.5) and VirB3 (9.0 -> 3.2). VirB1 is a peptidoglycan hydrolase
    (SLT, PF01464) and is genuinely OPTIONAL in many T4SS. So: among plasmids
    lacking VirB1, are the other six core anchors present? If yes, that is real
    absence, not detection failure, and VirB1 must not be a tiering requirement.
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
import decontaminate as DC

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CORE = ["VirD4", "VirB4", "VirB6", "VirB8", "VirB9", "VirB10", "VirB11"]


def main():
    from assertions import require_compute_node
    require_compute_node()          # A11: no heavy scans on a login node
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=8)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner, idx, start, end, strand, aas = [], [], [], [], [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                owner.append(f[0]); idx.append(int(f[1]))
                start.append(int(f[2])); end.append(int(f[3]))
                strand.append(int(f[4])); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    accs = sorted(set(owner))
    N = len(accs)
    ncds = collections.Counter(owner)
    print("cached proteins %d over %d accessions" % (len(seqs), N), flush=True)

    hit = {}
    for name, fam, ev, minaa in DC.ANCHOR_FAMILIES:
        p = os.path.join(HMM, fam + ".hmm")
        if not os.path.exists(p):
            continue
        with pyhmmer.plan7.HMMFile(p) as fh:
            model = next(iter(fh))
        kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
        best = {}
        for top in pyhmmer.hmmsearch([model], seqs, cpus=a.cpus, **kw):
            for h in top:
                i = int(dec(h.name))
                if minaa and len(aas[i]) <= minaa:
                    continue
                acc = owner[i]
                if acc not in best or h.score > best[acc][1]:
                    best[acc] = (i, h.score)
        hit[name] = best
        print("  %-8s %5d plasmids" % (name, len(best)), flush=True)
    anchors = list(hit)

    rows = []
    for acc in accs:
        r = {"accession": acc, "n_cds": ncds[acc]}
        ncore = 0
        for nm in anchors:
            b = hit[nm].get(acc)
            if b:
                i, sc = b
                r[nm] = "ga:%d:%d..%d:%s" % (idx[i], start[i], end[i],
                                             "+" if strand[i] == 1 else "-")
                r[nm + "_score"] = round(sc, 1)
                if nm in CORE:
                    ncore += 1
            else:
                r[nm] = "none"
                r[nm + "_score"] = ""
        r["core_completeness"] = "%d/%d" % (ncore, len(CORE))
        r["n_core"] = ncore
        rows.append(r)
    cols = (["accession", "n_cds", "core_completeness", "n_core"] +
            [c for nm in anchors for c in (nm, nm + "_score")])
    dest = os.path.join(PROJ, "data", "anchors", "anchor_matrix.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("wrote %s (%d rows)" % (dest, len(rows)))

    print("\n=== WHICH core anchor is missing in the non-7/7 plasmids? ===")
    incomplete = [r for r in rows if r["n_core"] < len(CORE)]
    print("  plasmids below 7/7: %d (%.1f%%)" % (len(incomplete), 100.0 * len(incomplete) / N))
    miss = collections.Counter()
    for r in incomplete:
        for nm in CORE:
            if r[nm] == "none":
                miss[nm] += 1
    for nm, c in miss.most_common():
        print("    %-8s missing in %5d of them (%.1f%%)" % (nm, c, 100.0 * c / len(incomplete)))
    pat = collections.Counter(
        tuple(sorted(nm for nm in CORE if r[nm] == "none")) for r in incomplete)
    print("\n  most common missing-anchor COMBINATIONS:")
    for p2, c in pat.most_common(8):
        print("    %-46s %5d" % (",".join(p2) or "(none)", c))
    print("\n  Concentrated in one or two anchors -> still a model problem.")
    print("  Spread across combinations -> genuinely incomplete systems.")

    print("\n=== Is VirB1 absence real, or detection failure? ===")
    nov1 = [r for r in rows if r.get("VirB1") == "none"]
    print("  plasmids without VirB1 at GA: %d (%.1f%%)" % (len(nov1), 100.0 * len(nov1) / N))
    other = [nm for nm in CORE]
    full = sum(1 for r in nov1 if all(r[nm] != "none" for nm in other))
    print("  of those, ALL %d core anchors present: %d (%.1f%%)"
          % (len(other), full, 100.0 * full / max(1, len(nov1))))
    print("  VirB1 is a peptidoglycan hydrolase and is optional in many T4SS.")
    print("  A high figure here means real absence -- VirB1 must NOT gate tiering.")
    dist = collections.Counter(r["n_core"] for r in nov1)
    print("  core completeness among VirB1-negative plasmids: %s"
          % ", ".join("%d/7:%d" % (k, v) for k, v in sorted(dist.items(), reverse=True)))


if __name__ == "__main__":
    main()
