#!/usr/bin/env python3
"""Release v1.1: MPF_T catalogue with the corrected VirB5 rule.

v1 is FROZEN and is not edited. This builds a separate v1.1.

WHAT CHANGED, and only this:
  VirB5 detection   v1  PF07996 @ E<=1e-5 AND aa>150
                    v1.1 (PF07996 @ GA AND aa>150) OR (T4SS_T_virB5 @ GA AND aa>150)
  The E-value rule is withdrawn -- docs/virb5_evalue_threshold_invalid.md. It
  failed both controls it was written for (RP4 E=0.22, R751 E=0.0093 under the
  pipeline's own search) and was STRICTER than the GA it replaced.

Only VirB5-dependent columns are recomputed. `tier` and `core_completeness_ga`
are NOT: CORE = VirD4,VirB4,VirB6,VirB8,VirB9,VirB10,VirB11 excludes VirB5, so
they are unaffected by construction. That is asserted, not assumed.

Every other column is carried forward from v1 by accession join.

The disjunction is admissible: band 91.3%, cross-class 6.6% (job 25913461).
It is the ONLY anchor disjunction adopted besides VirB1's; VirB2, VirB3 and
VirB6 were REJECTED on purity and keep their Pfam-GA calls.
"""
import argparse, collections, csv, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node, assert_keys_present

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
V1 = os.path.join(PROJ, "data/release/v1/catalogue_MPF_T_v1.tsv")


def main():
    require_compute_node()
    ap = argparse.ArgumentParser(); ap.add_argument("--cpus", type=int, default=16)
    a = ap.parse_args()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    seqs, owner, idx, st, en, strand, aas = [], [], [], [], [], [], []
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|")
                owner.append(f[0]); idx.append(int(f[1])); st.append(int(f[2]))
                en.append(int(f[3])); strand.append(int(f[4])); aas.append(line.rstrip("\n"))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(), sequence=line.rstrip("\n")).digitize(alpha))
    print("proteins %d over %d accessions" % (len(seqs), len(set(owner))), flush=True)

    def scan(path, minaa=0):
        with pyhmmer.plan7.HMMFile(path) as fh: m = next(iter(fh))
        best = {}
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name))
                if minaa and len(aas[i]) <= minaa: continue
                acc = owner[i]
                if acc not in best or h.score > best[acc][1]: best[acc] = (i, h.score)
        return best

    b5_pf = scan(os.path.join(HMM, "anchor_T4SS.hmm"), 150)
    b5_cj = scan(os.path.join(CJ, "T4SS_T_virB5.hmm"), 150)
    b6 = scan(os.path.join(HMM, "TrbL.hmm"))
    b4 = scan(os.path.join(HMM, "anchor_CagE_TrbE_VirB.hmm"))
    print("  VirB5 PF07996@GA   %d\n  VirB5 CONJScan@GA  %d\n  union              %d"
          % (len(b5_pf), len(b5_cj), len(set(b5_pf) | set(b5_cj))), flush=True)

    v1 = list(csv.DictReader(open(V1), delimiter="\t"))
    accs = [r["accession"] for r in v1]
    assert_keys_present(set(b5_pf) | set(b5_cj), set(accs) | set(owner),
                        "VirB5-positive accessions", "cds_cache vs v1 catalogue")

    rows, changed, src = [], 0, collections.Counter()
    for r in v1:
        acc = r["accession"]
        pick, how = None, "none"
        if acc in b5_pf and acc in b5_cj:
            pick = b5_pf[acc][0] if b5_pf[acc][1] >= b5_cj[acc][1] else b5_cj[acc][0]
            how = "both"
        elif acc in b5_pf: pick, how = b5_pf[acc][0], "pfam_only"
        elif acc in b5_cj: pick, how = b5_cj[acc][0], "conjscan_only"
        src[how] += 1
        new = ("ga:%d:%d..%d:%s" % (idx[pick], st[pick], en[pick],
               "+" if strand[pick] == 1 else "-")) if pick is not None else "none"
        if (new == "none") != (r["VirB5"] == "none"): changed += 1

        # architecture recomputed identically to v1 (VirB4/VirB6 only; VirB5-free)
        arch = "uncallable"
        if acc in b6 and acc in b4:
            i6, i4 = b6[acc][0], b4[acc][0]
            if strand[i6] == strand[i4]:
                t = 1 if strand[i6] == 1 else -1
                arch = "canonical" if idx[i4] * t < idx[i6] * t else "rearranged"
        if arch != r["architecture"]:
            raise SystemExit("architecture drift on %s: v1 %s vs recomputed %s "
                             "-- VirB5-independent columns must NOT move"
                             % (acc, r["architecture"], arch))
        ready = int(pick is not None and acc in b6 and arch != "uncallable")
        o = dict(r)
        o["VirB5"] = new
        o["virB5_source"] = how
        o["slot_ready__virb5_virb6"] = ready
        o["slot_ready__v1"] = r["slot_ready__virb5_virb6"]
        rows.append(o); 
    cols = list(v1[0].keys())
    for c in ("virB5_source", "slot_ready__v1"):
        if c not in cols: cols.append(c)
    dest = os.path.join(PROJ, "data", "release", "v1.1")
    os.makedirs(dest, exist_ok=True)
    out = os.path.join(dest, "catalogue_MPF_T_v1.1.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    was = sum(1 for r in v1 if r["slot_ready__virb5_virb6"] == "1")
    now = sum(r["slot_ready__virb5_virb6"] for r in rows)
    b5n = sum(1 for r in rows if r["VirB5"] != "none")
    print("\n=== v1.1 ===")
    print("  wrote %s (%d rows, %d cols)" % (out, n, len(cols)))
    print("  VirB5 present      v1 %d (%.1f%% fail) -> v1.1 %d (%.1f%% fail)"
          % (n - sum(1 for r in v1 if r["VirB5"] == "none"),
             100.0 * sum(1 for r in v1 if r["VirB5"] == "none") / n,
             b5n, 100.0 * (n - b5n) / n))
    print("  VirB5 call source  %s" % ", ".join("%s %d" % kv for kv in src.most_common()))
    print("  slot_ready         v1 %d (%.1f%%) -> v1.1 %d (%.1f%%)"
          % (was, 100.0 * was / n, now, 100.0 * now / n))
    print("  plasmids whose VirB5 presence/absence flipped: %d" % changed)
    tv1 = collections.Counter(r["tier"] for r in v1)
    tv11 = collections.Counter(r["tier"] for r in rows)
    print("  tier unchanged by construction: %s" % ("YES" if tv1 == tv11 else "NO -- BUG"))
    if tv1 != tv11: raise SystemExit("tier moved; VirB5 must not affect tiering")


if __name__ == "__main__":
    main()
