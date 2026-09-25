#!/usr/bin/env python3
"""Do the NAMED exclusion families obey Class A or Class B?

  Class A -- lipoprotein: lipobox positive, signal peptide + lipidated Cys
  Class B -- TM-anchored: exactly one M segment, starting within the N-terminal
             40 residues, followed by a C-terminal O (periplasmic) region of
             >=30 aa, whole protein 60-200 aa

This script does CLASS A now. Class B needs per-residue topology (M/I/O labels),
which the KD-window counter CANNOT supply: it has no orientation, so it cannot
distinguish the C-terminal side from the N-terminal one, and it already failed its
own control gate (83.1% positives vs 22.5% size-matched negatives, short of the
pre-registered 80/20 bar). Reporting Class B from KD windows would be inventing
topology. Deferred to a real topology predictor.

The C-terminal >=30 aa requirement is the mechanistically load-bearing one: a
single-pass protein with a 5-residue C-terminus has no periplasmic domain to
recognise the VirB6 periplasmic loop with. That is exactly the criterion a
hydrophobicity window cannot test.

Lipobox definitions, all three reported, as calibrated on the MPF_T slot:
  strict      [LVI][ASTVIG][GASN]C
  relaxed     [LVIMFWY][ASTVIGN][GASNDQ]C
  structural  Cys at position <=25 preceded by a hydrophobic 8-mer

The regex is applied to the N-terminal 40 residues only -- a lipobox is a signal
peptide feature, and scanning the whole protein would count internal cysteines.
"""
import collections, csv, os, re, sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
SOF = os.path.join(PROJ, "data", "release", "v1.1", "slot_occupant_families.tsv")
FAM = {"TIGR04359": "TrbK", "PF20084": "TrbK_Pfam", "NF033894": "Eex_IncN",
       "NF041429": "EexR", "NF033891": "ExcA", "PF10624": "TraS",
       "PF14729": "DUF4467", "PF05818": "TraT", "PrgA_Sea1": "PrgA_Sea1"}
STRICT = re.compile(r"[LVI][ASTVIG][GASN]C")
RELAXED = re.compile(r"[LVIMFWY][ASTVIGN][GASNDQ]C")
HYDRO = set("AVLIMFWYC")
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
NTERM = 40


def structural(s):
    """Cys at <=25 with a hydrophobic 8-mer immediately before it."""
    for i, c in enumerate(s[:25]):
        if c == "C" and i >= 8:
            w = s[i - 8:i]
            if sum(1 for x in w if x in HYDRO) >= 6:
                return True
    return False


def maxkd(s, w=19):
    if len(s) < w:
        return round(sum(KD.get(c, 0) for c in s) / max(1, len(s)), 2)
    return round(max(sum(KD.get(c, 0) for c in s[i:i + w]) / w
                     for i in range(len(s) - w + 1)), 2)


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    rows = [r for r in csv.DictReader(open(SOF), delimiter="\t")
            if r["slot_occupant_coords"]]
    want = collections.defaultdict(set)
    for r in rows:
        want[r["accession"]].add(r["slot_occupant_coords"])
    seq_of = {}
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                if acc in want:
                    k = "%s..%s" % (f[2], f[3])
                    if k in want[acc]:
                        seq_of[(acc, k)] = line.rstrip("\n")
    occ = sorted({s for s in seq_of.values()})
    print("slot occupants: %d unique\n" % len(occ), flush=True)

    dig = [pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=s).digitize(alpha)
           for i, s in enumerate(occ)]
    byfam = collections.defaultdict(set)
    for acc, lab in FAM.items():
        p = os.path.join(HMM, acc + ".hmm")
        if not os.path.exists(p):
            print("  (model %s absent)" % acc); continue
        with pyhmmer.plan7.HMMFile(p) as fh:
            m = next(iter(fh))
        try:
            it = pyhmmer.hmmsearch([m], dig, cpus=8, bit_cutoffs="gathering")
        except Exception:
            it = pyhmmer.hmmsearch([m], dig, cpus=8, E=1e-5)
        try:
            for top in it:
                for h in top:
                    byfam[lab].add(int(dec(h.name)))
        except pyhmmer.errors.MissingCutoffs:
            for top in pyhmmer.hmmsearch([m], dig, cpus=8, E=1e-5):
                for h in top:
                    byfam[lab].add(int(dec(h.name)))

    print("=== CLASS A (lipoprotein): lipobox in the N-terminal %d residues ===" % NTERM)
    print("  %-11s %6s %8s %9s %9s %9s %9s" %
          ("family", "n", "median aa", "strict", "relaxed", "structural", "max KD"))
    out = []
    for lab in sorted(byfam, key=lambda x: -len(byfam[x])):
        idx = sorted(byfam[lab])
        if not idx:
            continue
        S = [occ[i] for i in idx]
        L = sorted(len(x) for x in S)
        st = sum(1 for x in S if STRICT.search(x[:NTERM]))
        rx = sum(1 for x in S if RELAXED.search(x[:NTERM]))
        sr = sum(1 for x in S if structural(x))
        kd = sorted(maxkd(x) for x in S)
        n = len(S)
        print("  %-11s %6d %8d %8.1f%% %8.1f%% %9.1f%% %9.2f"
              % (lab, n, L[n // 2], 100.0 * st / n, 100.0 * rx / n,
                 100.0 * sr / n, kd[n // 2]))
        out.append({"family": lab, "n_unique": n, "median_aa": L[n // 2],
                    "lipobox_strict_pct": round(100.0 * st / n, 1),
                    "lipobox_relaxed_pct": round(100.0 * rx / n, 1),
                    "lipobox_structural_pct": round(100.0 * sr / n, 1),
                    "median_max_kd": kd[n // 2],
                    "len_q1": L[n // 4], "len_q3": L[3 * n // 4],
                    "in_60_200_pct": round(100.0 * sum(1 for x in L if 60 <= x <= 200) / n, 1)})

    print("\n=== length window 60-200 aa (the Class B size constraint) ===")
    print("  %-11s %6s %7s %7s %10s" % ("family", "Q1", "median", "Q3", "in 60-200"))
    for r in out:
        print("  %-11s %6d %7d %7d %9.1f%%"
              % (r["family"], r["len_q1"], r["median_aa"], r["len_q3"], r["in_60_200_pct"]))

    print("\n=== reading ===")
    for r in out:
        a = r["lipobox_strict_pct"] >= 50
        print("  %-11s %s" % (r["family"],
              "CLASS A (lipoprotein): strict lipobox %.0f%%" % r["lipobox_strict_pct"] if a
              else "NOT class A (strict lipobox %.0f%%) -- Class B needs topology, deferred"
                   % r["lipobox_strict_pct"]))
    dest = os.path.join(PROJ, "data", "anchors", "class_A_named.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print("\nwrote %s" % dest)
    print("\n  CLASS B NOT TESTED HERE. It requires per-residue M/I/O topology:")
    print("  exactly one M segment, starting within 40 aa, followed by a C-terminal")
    print("  O region of >=30 aa. A KD window has no orientation and cannot say")
    print("  which side is periplasmic, so it cannot test the >=30 aa C-terminal")
    print("  requirement -- the one that carries the mechanism.")


if __name__ == "__main__":
    main()
