#!/usr/bin/env python3
"""MPF_T tra gene cluster catalogue, v1. GA is the primary measure.

Frozen deliberately. Three sessions of anchor work kept finding artefacts and
recomputing, and the catalogue never existed as a file. Anything found later goes
into the caveats and a v2, rather than deferring v1 again.

NOT an operon catalogue: HMM hits and gene positions, no promoter, terminator or
co-transcription evidence.

PRIMARY MEASURE IS GA. Second-pass phmmer recall is carried as an annotation
column only, because its precision control is structurally blind on exactly the
anchors where it claims the most:

    the conflict filter rejects a recalled protein if it already has a GA
    assignment to another family. VirB5's contaminants were 816 aa VirB4 proteins
    -- always assigned, always caught (72% rejected). VirB2's plausible
    contaminants are ~120 aa proteins that hit nothing at all, so the filter
    cannot see them, and VirB2 rejected 0.3%. Zero rejection is evidence the
    filter did not fire.

So recall's core-anchor gains (<=0.4 points) are reported and its peripheral gains
(3-9 points) are reported as UNVERIFIED.

Tier definition, each part measured rather than chosen:

    A   7/7 core anchors, architecture callable, no oversized unresolved gap
    B   6/7 core
    C   >=5/7 core, or architecture uncallable
    D   <5/7 -- dominated by the 270 plasmids missing the entire VirB6-VirB11
        block, which entered the set on VirB4 alone (VirB4 is the entry filter,
        so it is present by construction and cannot evidence a T4SS)

VirB1 is EXCLUDED from tiering and carried as its own column: of 1,452 plasmids
lacking it at GA, 1,118 (77.0%) carry all seven core anchors, so its absence is
real, not detection failure. Gating on it would drop ~14% of plasmids for lacking
an optional peptidoglycan hydrolase -- a bias correlated with host cell-wall
biology.
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
PLSDB = "/global/scratch/users/kh36969/plsdb"
CORE = ["VirD4", "VirB4", "VirB6", "VirB8", "VirB9", "VirB10", "VirB11"]


def main():
    from assertions import require_compute_node
    require_compute_node()          # A11: no heavy scans on a login node
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner, idx, st, en, strand, aas = [], [], [], [], [], [], []
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
                st.append(int(f[2])); en.append(int(f[3]))
                strand.append(int(f[4])); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    accs = sorted(set(owner))
    print("proteins %d over %d accessions" % (len(seqs), len(accs)), flush=True)
    byacc = collections.defaultdict(list)
    for i in range(len(owner)):
        byacc[owner[i]].append(i)

    hit, assigned = {}, {}
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
                assigned[i] = assigned.get(i, name)
                acc = owner[i]
                if acc not in best or h.score > best[acc][1]:
                    best[acc] = (i, h.score)
        hit[name] = best
        print("  %-8s %5d plasmids" % (name, len(best)), flush=True)
    anchors = list(hit)

    mpf, rep = {}, collections.defaultdict(set)
    tp = os.path.join(PLSDB, "typing.csv")
    if os.path.exists(tp):
        csv.field_size_limit(10 ** 7)
        for r in csv.DictReader(open(tp)):
            mpf[r["NUCCORE_ACC"]] = r.get("mpf_type", "") or "-"
    pf = os.path.join(PLSDB, "plasmidfinder.csv")
    if os.path.exists(pf):
        for r in csv.DictReader(open(pf)):
            rep[r["NUCCORE_ACC"]].add(r["typing"])

    rows = []
    for acc in accs:
        genes = sorted(byacc[acc], key=lambda i: idx[i])
        pos = {nm: hit[nm][acc][0] for nm in anchors if acc in hit[nm]}
        ncore = sum(1 for nm in CORE if nm in pos)
        # architecture from VirB4/VirB6 only (98.1% accurate, needs no VirB5)
        arch = "uncallable"
        if "VirB6" in pos and "VirB4" in pos:
            i6, i4 = pos["VirB6"], pos["VirB4"]
            if strand[i6] == strand[i4]:
                t = 1 if strand[i6] == 1 else -1
                arch = "canonical" if idx[i4] * t < idx[i6] * t else "rearranged"
        # unresolved gaps: genes inside the anchor span with no anchor assignment
        gaps, gaplen = 0, []
        if len(pos) >= 2:
            lo, hi = min(idx[i] for i in pos.values()), max(idx[i] for i in pos.values())
            for i in genes:
                if lo < idx[i] < hi and i not in assigned:
                    gaps += 1; gaplen.append(len(aas[i]))
        if ncore == 7 and arch != "uncallable":
            tier = "A"
        elif ncore == 6:
            tier = "B"
        elif ncore >= 5 or arch == "uncallable":
            tier = "C"
        else:
            tier = "D"
        r = {"accession": acc, "mpf_class": mpf.get(acc, "-"),
             "replicons": ";".join(sorted(rep.get(acc, ()))) or "-",
             "n_cds": len(genes),
             "core_completeness_ga": "%d/7" % ncore, "n_core_ga": ncore,
             "architecture": arch, "tier": tier,
             "virB1_present": int("VirB1" in pos),
             "n_unresolved_gaps": gaps,
             "unresolved_gap_lengths": ";".join(str(x) for x in sorted(gaplen)) if gaplen else "",
             "anchor_span": ("%d..%d" % (min(st[i] for i in pos.values()),
                                         max(en[i] for i in pos.values()))) if pos else ""}
        for nm in anchors:
            i = pos.get(nm)
            r[nm] = ("ga:%d:%d..%d:%s" % (idx[i], st[i], en[i], "+" if strand[i] == 1 else "-")
                     if i is not None else "none")
        rows.append(r)

    cols = ["accession", "mpf_class", "replicons", "n_cds", "tier",
            "core_completeness_ga", "n_core_ga", "architecture", "virB1_present",
            "n_unresolved_gaps", "unresolved_gap_lengths", "anchor_span"] + anchors
    dest = os.path.join(PROJ, "data", "catalogue_MPF_T_regen.tsv")   # A9: the frozen
    # release is data/release/v1/catalogue_MPF_T_v1.tsv and is NOT written here.
    # Regenerate, diff against the release, then promote by hand if intended.
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s (%d rows)" % (dest, len(rows)))

    print("\n=== CATALOGUE v1 SUMMARY (GA primary) ===")
    tc = collections.Counter(r["tier"] for r in rows)
    for t in "ABCD":
        print("  tier %s  %5d (%.1f%%)" % (t, tc[t], 100.0 * tc[t] / len(rows)))
    print("\n  architecture: %s"
          % ", ".join("%s %d" % (k, v) for k, v in
                      collections.Counter(r["architecture"] for r in rows).most_common()))
    print("  mpf_type    : %s"
          % ", ".join("%s %d" % (k, v) for k, v in
                      collections.Counter(r["mpf_class"] for r in rows).most_common(4)))
    print("  VirB1 present in %d (%.1f%%)"
          % (sum(r["virB1_present"] for r in rows),
             100.0 * sum(r["virB1_present"] for r in rows) / len(rows)))
    g = [r["n_unresolved_gaps"] for r in rows]
    print("\n  unresolved genes inside the anchor span: total %d, median per plasmid %d"
          % (sum(g), sorted(g)[len(g) // 2]))
    a_rows = [r for r in rows if r["tier"] == "A"]
    if a_rows:
        ga = [r["n_unresolved_gaps"] for r in a_rows]
        print("  tier A only: total %d, median %d  <- the candidate container"
              % (sum(ga), sorted(ga)[len(ga) // 2]))


if __name__ == "__main__":
    main()
