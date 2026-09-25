#!/usr/bin/env python3
"""STEP 0 -- are the four known exclusion genes actually in our data?

Blocking prerequisite for the covariation module. Every downstream statistic is
computed over called ORFs, so if the caller or the catalogue silently drops a
known EEx, the covariation scan is measuring a set that cannot contain its own
positive controls. That is the dangerous failure class already recorded here: a
successful run producing a meaningless number.

Each control has a SPECIFIC, predictable way to go missing, and they are
different failures with different consequences:

  TrbK   RP4     69 aa   -- short ORF. MIN_GENE is 90 nt = 30 aa so it is above
                            the floor, but anon-mode scoring can still drop it.
  ExcB   R64    147 aa   -- overlaps excA, translated by reinitiation. Gene
                            callers are unreliable on overlapping ORFs and will
                            usually emit only one of the pair.
  YddJ   ICEBs1 126 aa   -- lives on an ICE. If a catalogue was built from
                            plasmids only, this control is absent by construction
                            and its absence says nothing about the caller.
  TraS   F/R64          -- unmodelled (PF10624 is R100-specific, 1/14 seeds), so
                            they must be fetched by COORDINATE, never by name.

TWO TESTS, kept separate because they fail for different reasons:

  A. CALLER    run the frozen caller on the reference sequence and ask whether
               the control gene comes back. This is the one that generalises:
               a blind spot for short or overlapping ORFs applies equally to
               every undiscovered EEx in G/B/C.
  B. CATALOGUE is the reference accession in cds_cache_full at all, and if so at
               which gene index. A miss here is dataset COVERAGE, not a caller
               defect, and must not be reported as one.

MATCHING IS BY STOP COORDINATE, not by sequence identity. This project already
established that as the stable identity of a called gene: two callers that
disagree about the start still agree about the stop. Matching on identity alone
would score a correct call with a different start site as a miss.

Test A is an ASSERTION. Any control the caller fails to recover is a hard failure
and the module does not proceed.
"""
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node
import orf_caller

GB = os.path.join(PROJ, "data", "seed", "genbank")
HMM = os.path.join(PROJ, "data", "hmm")
CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"

# element -> (file, [(control, why at risk)], partners, catalogue accessions)
#
# PLSDB is keyed by REFSEQ accession. Looking a reference up by its GenBank
# accession returns "not in cache" for a plasmid that is plainly there -- R64 is
# NC_005014.1, R100 is NC_002134.1. Same shape as the two join traps already
# recorded in this project (gene names are not keys; accession_related is a
# different coordinate frame), so every reference carries its alternates here.
REFS = [
    ("RP4",    "RP4__BN000925.1.gb",    [("trbK", "69 aa short ORF")],
     ["trbL", "trbB"], ["BN000925.1", "CP152305.1", "NZ_PP591959.1"]),
    ("R64",    "R64__AP005147.1.gb",    [("excA", "reference"),
                                         ("excB", "147 aa, overlaps excA"),
                                         ("traS", "62 aa, unmodelled")],
     ["traY"], ["AP005147.1", "NC_005014.1"]),
    ("F",      "F__AP001918.1.gb",      [("traS", "173 aa, unmodelled")],
     ["traG"], ["AP001918.1", "U01159.2", "X06915.1"]),
    ("R100",   "R100__AP000342.1.gb",   [("traS", "159 aa, the ONE PF10624 hit")],
     ["traG"], ["AP000342.1", "NC_002134.1"]),
    ("ICEBs1", "ICEBs1__CP171645.1.gb", [("@PF14729", "YddJ; ICE, may be absent "
                                          "from a plasmid-built catalogue")],
     ["conG"], ["CP171645.1"]),
]


def annotated(path):
    """(gene_name, start, end, strand, aa) for every annotated CDS."""
    from Bio import SeqIO
    out, seqs = [], []
    for rec in SeqIO.parse(path, "genbank"):
        seqs.append((rec.id, str(rec.seq).upper()))
        for f in rec.features:
            if f.type != "CDS":
                continue
            p = f.qualifiers.get("translation", [None])[0]
            if not p:
                continue
            nm = (f.qualifiers.get("gene") or f.qualifiers.get("locus_tag") or ["?"])[0]
            out.append((nm, int(f.location.start) + 1, int(f.location.end),
                        1 if f.location.strand == 1 else -1, p))
    return out, seqs


def by_hmm(cds, model):
    """Locate a control that has no usable gene name, via its HMM."""
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x
    path = os.path.join(HMM, model + ".hmm")
    if not os.path.exists(path):
        return None
    seqs = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=c[4]).digitize(alpha)
            for i, c in enumerate(cds)]
    with pyhmmer.plan7.HMMFile(path) as fh:
        hmm = next(iter(fh))
    best = None
    for top in pyhmmer.hmmsearch([hmm], seqs, cpus=4, bit_cutoffs="gathering"):
        for h in top:
            if best is None or h.score > best[0]:
                best = (h.score, int(dec(h.name)))
    return cds[best[1]] if best else None


def ident(a, b):
    n = min(len(a), len(b))
    if not n:
        return 0.0
    # compare from the C-terminus: callers agree on the stop, differ on the start
    return 100.0 * sum(x == y for x, y in zip(a[-n:], b[-n:])) / n


def cache_index(accs):
    """Gene indices for these accessions in cds_cache_full, if present at all."""
    accs = set(accs)
    got = {}
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            elif nm.split("|")[0] in accs:
                f = nm.split("|")
                got[(int(f[2]), int(f[3]))] = (int(f[1]), line.rstrip("\n"))
    return got


def main():
    require_compute_node()
    orf_caller.assert_anon_only(orf_caller.ANON)
    print("frozen caller: anon=%s min_gene=%d nt (%d aa) closed=%s table=%d\n"
          % (orf_caller.ANON, orf_caller.MIN_GENE, orf_caller.MIN_GENE // 3,
             orf_caller.CLOSED, orf_caller.TRANSL_TABLE), flush=True)

    rows, failures, coverage = [], [], []
    for label, fn, controls, partners, cat_accs in REFS:
        path = os.path.join(GB, fn)
        if not os.path.exists(path):
            failures.append("%s: reference file missing" % label); continue
        cds, seqs = annotated(path)
        acc = seqs[0][0]
        called = []
        for _, s in seqs:
            called.extend(orf_caller.call(s))
        # stop coordinate is the stable key
        by_stop = {(g["stop_coord"], g["strand"]): g for g in called}
        cache = cache_index(cat_accs)
        in_cache = len(cache) > 0
        print("=" * 72)
        print("%s  (ref %s; catalogue alts %s)\n   %d annotated CDS | %d called"
              " | in cds_cache_full: %s"
              % (label, acc, ",".join(cat_accs), len(cds), len(called),
                 "YES" if in_cache else "NO"))
        print("=" * 72)
        if not in_cache:
            coverage.append("%s not in cds_cache_full under any of %s"
                            % (label, cat_accs))

        for name, why in controls:
            if name.startswith("@"):
                rec = by_hmm(cds, name[1:])
                disp = "%s(%s)" % (name[1:], rec[0] if rec else "?")
                if rec is None:
                    failures.append("%s: %s -- model found nothing to test"
                                    % (label, name)); continue
            else:
                hit = [c for c in cds if c[0].lower() == name.lower()]
                if not hit:
                    failures.append("%s: %s not annotated in the reference"
                                    % (label, name)); continue
                rec, disp = hit[0], name
            nm, st, en, sd, aa = rec
            stop = en if sd == 1 else st
            g = by_stop.get((stop, sd))
            # TEST A -- the caller
            if g is None:
                verdict = "MISSING FROM CALLER"
                failures.append("%s / %s (%d aa, %s): caller did not emit this ORF"
                                % (label, disp, len(aa), why))
                idn, gi = 0.0, ""
            else:
                idn = ident(aa, g["aa"])
                verdict = "ok" if idn >= 99.0 else "start differs (%.0f%% id)" % idn
                if idn < 99.0:
                    verdict = "CALLED, %s" % verdict
                gi = ""
            # TEST B -- the catalogue
            cell = cache.get((st, en)) if in_cache else None
            if not cell and in_cache:
                # coordinates may differ by start site; match on stop
                for (cs, ce), v in cache.items():
                    if (ce if sd == 1 else cs) == stop:
                        cell = v; break
            cat = ("gene_index %d" % cell[0]) if cell else ("--" if in_cache else "n/a")
            print("  %-14s %4d aa  caller: %-28s catalogue: %s"
                  % (disp, len(aa), verdict, cat))
            rows.append({"element": label, "accession": acc, "control": disp,
                         "aa": len(aa), "why_at_risk": why,
                         "caller_recovered": int(g is not None),
                         "identity_pct": round(idn, 1),
                         "in_cds_cache": int(in_cache),
                         "catalogue_gene_index": cell[0] if cell else ""})

        for p in partners:
            hit = [c for c in cds if c[0].lower() == p.lower()]
            if hit:
                nm, st, en, sd, aa = hit[0]
                stop = en if sd == 1 else st
                print("  %-14s %4d aa  (partner) caller: %s"
                      % (p, len(aa), "ok" if (stop, sd) in by_stop else "MISSING"))

    print("\n" + "=" * 72)
    if coverage:
        print("COVERAGE NOTES (not caller defects):")
        for c in coverage:
            print("  - %s" % c)
    dest = os.path.join(PROJ, "data", "anchors", "step0_controls.tsv")
    if rows:
        with open(dest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(rows)
        print("\nwrote %s (%d rows)" % (dest, len(rows)))

    if failures:
        print("\n*** STEP 0 FAILED -- %d control(s) unrecoverable ***" % len(failures))
        for f in failures:
            print("  %s" % f)
        print("\nThe covariation module MUST NOT proceed. A caller blind spot for")
        print("short or overlapping ORFs applies equally to every undiscovered EEx")
        print("in MPF_G/B/C, so it is a larger problem than the scan it blocks.")
        return 1
    print("\nSTEP 0 PASSED -- every control is recoverable by the frozen caller.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
