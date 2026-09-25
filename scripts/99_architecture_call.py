#!/usr/bin/env python3
"""Determine skeleton architecture WITHOUT VirB5, then de-bias the nomination set.

The nomination rule is "eex immediately 3' of VirB5", but VirB5 fails on 27.1% of
VirB4+ plasmids -- 1,286 with `no_VirB5_only` -- and that loss is not random: it
enriches for divergent systems, i.e. exactly where a novel family would be. Round 1
would search a biased set.

The fix follows from the architecture result. In canonical order VirB4 precedes
VirB6; in the rearranged (IncI2-type) order VirB6 precedes VirB4:

    canonical  ... VirB3 VirB4 [VirB5 eex VirB6] VirB8 VirB9 VirB10 VirB11 ...
    rearranged ... [VirB6 VirB5 eex] VirB3 VirB4 VirB8 VirB9 VirB10 ...

So the SIGN of (VirB6 - VirB4) in transcription order calls the architecture using
two anchors that never need VirB5: VirB4 fails 0% by construction (it is the entry
criterion) and VirB6 fails 4.9%.

Given the architecture, the eex position follows from VirB6 alone:

    canonical  -> eex at -1 from VirB6
    rearranged -> eex at +2 from VirB6

The rule is validated against the plasmids where VirB5 IS detectable (so the
architecture is already known from the VirB5/VirB6 order) before it is applied to
the 1,286 where it is not.
"""
import argparse
import collections
import csv
import glob
import os
import re
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
LIPO = re.compile(r"[LVI][ASTVIG][GASN]C")
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
      'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,
      'Y':-1.3,'V':4.2}


def hreg(s):
    for i, c in enumerate(s[:25]):
        if c == "C" and i >= 8 and sum(KD.get(x, 0) for x in s[i-8:i]) / 8.0 >= 1.0:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner, idx, strand, aas = [], [], [], [], []
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
                strand.append(int(f[4])); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(),
                    sequence=line.rstrip("\n")).digitize(alpha))
    print("cached proteins %d over %d accessions" % (len(seqs), len(set(owner))), flush=True)
    bypos = {(owner[i], idx[i]): i for i in range(len(owner))}

    def scan(fam, ev=None, minaa=0):
        with pyhmmer.plan7.HMMFile(os.path.join(HMM, fam + ".hmm")) as fh:
            m = next(iter(fh))
        kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
        out = collections.defaultdict(list)
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, **kw):
            for h in top:
                i = int(dec(h.name))
                if minaa and len(aas[i]) <= minaa:
                    continue
                out[owner[i]].append((idx[i], strand[i], h.score, i))
        return out
    B4 = scan("anchor_CagE_TrbE_VirB")
    B5 = scan("anchor_T4SS", 1e-5, 150)
    B6 = scan("TrbL")
    print("VirB4 %d  VirB5 %d  VirB6 %d plasmids" % (len(B4), len(B5), len(B6)), flush=True)

    gap = {}
    for f in glob.glob(os.path.join(PROJ, "data", "positional", "virb5_virb6_gap_*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            gap[r["accession"]] = r.get("order", "")

    def best(d, acc):
        return max(d[acc], key=lambda x: x[2]) if acc in d else None

    # ---- architecture call from VirB4/VirB6 only ---------------------------
    call = {}
    for acc in B6:
        b6, b4 = best(B6, acc), best(B4, acc)
        if not b4 or b4[1] != b6[1]:
            continue
        t = lambda i: i * (1 if b6[1] == 1 else -1)
        call[acc] = "canonical" if t(b4[0]) < t(b6[0]) else "rearranged"

    print("\n=== VALIDATION: architecture called from VirB4/VirB6 vs the known label ===")
    tab = collections.Counter()
    for acc, lab in gap.items():
        if lab in ("canonical", "inverted") and acc in call:
            tab[(lab, call[acc])] += 1
    tot_c = tab[("canonical", "canonical")] + tab[("canonical", "rearranged")]
    tot_i = tab[("inverted", "canonical")] + tab[("inverted", "rearranged")]
    print("  known canonical -> called canonical  %5d / %5d = %.1f%%"
          % (tab[("canonical", "canonical")], tot_c,
             100.0 * tab[("canonical", "canonical")] / max(1, tot_c)))
    print("  known inverted  -> called rearranged %5d / %5d = %.1f%%"
          % (tab[("inverted", "rearranged")], tot_i,
             100.0 * tab[("inverted", "rearranged")] / max(1, tot_i)))
    acc_overall = (tab[("canonical", "canonical")] + tab[("inverted", "rearranged")]) / max(1, tot_c + tot_i)
    print("  overall accuracy %.1f%%   (uses NO VirB5)" % (100.0 * acc_overall))

    # ---- apply to the VirB5-negative plasmids ------------------------------
    novb5 = [acc for acc in B6 if acc not in B5]
    print("\n=== APPLYING to the %d plasmids with VirB6 but no VirB5 ===" % len(novb5))
    cc = collections.Counter(call.get(a2, "uncallable") for a2 in novb5)
    for k, v in cc.most_common():
        print("  %-12s %5d (%.1f%%)" % (k, v, 100.0 * v / len(novb5)))

    OFFSET = {"canonical": -1, "rearranged": +2}
    nominated, rows = {}, []
    for acc in novb5:
        arch = call.get(acc)
        if arch not in OFFSET:
            continue
        b6 = best(B6, acc)
        t = 1 if b6[1] == 1 else -1
        j = bypos.get((acc, b6[0] + OFFSET[arch] * t))
        if j is None or strand[j] != b6[1]:
            continue
        nominated.setdefault(aas[j], []).append(acc)
        rows.append({"accession": acc, "architecture": arch,
                     "offset_from_virb6": OFFSET[arch], "aa_len": len(aas[j]),
                     "lipobox_strict": int(bool(LIPO.search(aas[j][:40]))),
                     "lipobox_hregion": int(hreg(aas[j])), "aa_seq": aas[j]})
    print("  nominated: %d records, %d unique proteins" % (len(rows), len(nominated)))
    if rows:
        import statistics
        L = sorted(r["aa_len"] for r in rows)
        print("  length: median %d  Q1 %d  Q3 %d" % (statistics.median(L), L[len(L)//4], L[3*len(L)//4]))
        for k, lab in (("lipobox_strict", "strict"), ("lipobox_hregion", "h-region")):
            print("  lipobox (%-8s): records %.1f%%   unique %.1f%%"
                  % (lab, 100.0 * sum(r[k] for r in rows) / len(rows),
                     100.0 * sum(1 for s in nominated if (LIPO.search(s[:40]) if k.endswith("strict") else hreg(s))) / len(nominated)))
        band = sum(1 for r in rows if 60 <= r["aa_len"] <= 100)
        print("  in the 60-100 aa control band: %d records (%.1f%%)" % (band, 100.0 * band / len(rows)))
        eex = scan("NF033894")
        named = sum(1 for r in rows
                    if any(i2 == bypos.get((r["accession"], 0), -1) for i2 in []) )
        hits_at = 0
        for r in rows:
            if r["accession"] in eex:
                for i2, _, _, gi in eex[r["accession"]]:
                    if aas[gi] == r["aa_seq"]:
                        hits_at += 1; break
        print("  already hit by NF033894 at GA: %d records (%.1f%%)"
              % (hits_at, 100.0 * hits_at / len(rows)))
    dest = os.path.join(PROJ, "data", "positional", "nominated_no_virb5.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "architecture", "offset_from_virb6",
                                           "aa_len", "lipobox_strict", "lipobox_hregion", "aa_seq"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
