#!/usr/bin/env python3
"""MPF_FA catalogue. 2,391 plasmids.

ENTRY (docs/mpffa_entry_criterion.md):
    T4SS_virb4 AND NOT (MPF_T or MPF_F or MPF_I) AND >=2 FA/FATA profiles
Cross-class is EXCLUDED BY CONSTRUCTION, not measured -- no 0% figure is quoted.

*** THERE IS NO SLOT COLUMN, DELIBERATELY. ***

The conG+3 slot is confirmed on ICEBs1 and nowhere else:
  - seed check (job 26035629): 1 of 4 locatable seeds. pAD1's +2 occupant is a
    type II toxin-antitoxin system (antitoxin MazE at +1); pLS20's and pAM373's
    are ANTI-ORIENTED to the anchor and so not plausibly co-transcribed; pCF10
    has no anchor at all.
  - class census (job 26038132): DUF4467 occurs on 11 of 2,391 admitted plasmids
    (0.5%), against TraS 2,725 / ExcA 2,213 / TrbK 224 for the other three
    classes. Within those 11 the modal offset from FA_orf15 is -1, NOT ICEBs1's
    +3 -- the one known case is not the modal case.

A slot column validated on a single element would not be a caveated column, it
would be an empty one. It is omitted rather than shipped with a warning.

TIER: the specified definition was "entry hit + slot-flanking anchors present".
The second half is impossible with no slot, so tier is renamed `profile_support`
and counts FA/FATA profiles -- the only completeness signal this class has.
It is NOT called A/B/C/D: those letters mean 7/7 structural core in MPF_T and
would import a meaning that does not exist here. Every CONJScan Gram-positive
class-specific profile is `accessory`, so there is no mandatory core to count.
"""
import collections, csv, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
ADM = os.path.join(PROJ, "data", "anchors", "mpffa_admitted.tsv")
MAND = ["T4SS_virb4", "T4SS_tcpA", "T4SS_t4cp1", "T4SS_t4cp2"] + \
       ["T4SS_MOB" + x for x in ("B", "C", "F", "H", "P1", "P2", "P3", "Q", "T", "V")]
EEX = ["PF14729", "PF10624", "NF033891", "NF033894", "NF041429", "TIGR04359"]


def main():
    require_compute_node()
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    dec = lambda x: x.decode() if isinstance(x, bytes) else x

    adm = {}
    for i, l in enumerate(open(ADM)):
        if i:
            f = l.rstrip("\n").split("\t"); adm[f[0]] = int(f[1])
    print("admitted: %d" % len(adm), flush=True)

    seqs, owner, idx, st, en, strand = [], [], [], [], [], []
    ncds = collections.Counter()
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"): continue
        nm = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">": nm = line[1:].rstrip("\n")
            else:
                f = nm.split("|"); acc = f[0]
                ncds[acc] += 1
                if acc not in adm: continue
                owner.append(acc); idx.append(int(f[1])); st.append(int(f[2]))
                en.append(int(f[3])); strand.append(int(f[4]))
                seqs.append(pyhmmer.easel.TextSequence(
                    name=str(len(seqs)).encode(), sequence=line.rstrip("\n")).digitize(alpha))
    for k, s in enumerate(seqs): s.name = str(k).encode()
    print("subset: %d proteins over %d plasmids\n" % (len(seqs), len(set(owner))), flush=True)

    def scan(path):
        with pyhmmer.plan7.HMMFile(path) as fh: m = next(iter(fh))
        best = {}
        for top in pyhmmer.hmmsearch([m], seqs, cpus=16, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name)); a = owner[i]
                if a not in best or h.score > best[a][1]: best[a] = (i, round(h.score, 1))
        return best

    fa = sorted(f[:-4] for f in os.listdir(CJ) if f.startswith("FA_") and f.endswith(".hmm"))
    fata = sorted(f[:-4] for f in os.listdir(CJ) if f.startswith("FATA_") and f.endswith(".hmm"))
    hit = {}
    for p in fa + fata + MAND:
        q = os.path.join(CJ, p + ".hmm")
        if os.path.exists(q): hit[p] = scan(q)
    for p in EEX:
        q = os.path.join(HMM, p + ".hmm")
        if os.path.exists(q): hit[p] = scan(q)
    print("scanned %d profiles\n" % len(hit), flush=True)

    rows = []
    for acc in sorted(adm):
        f_hits = [p for p in fa if acc in hit.get(p, {})]
        t_hits = [p for p in fata if acc in hit.get(p, {})]
        nf, nt = len(f_hits), len(t_hits)
        if max(nf, nt) < 2: sub = "unassigned"
        elif nf >= 2 * max(1, nt): sub = "FA"
        elif nt >= 2 * max(1, nf): sub = "FATA"
        else: sub = "unassigned"
        n = nf + nt
        r = {"accession": acc, "mpf_class_call": "MPF_FA",
             "entry_criterion": "T4SS_virb4 AND NOT(T/F/I) AND >=2 FA/FATA",
             "n_cds": ncds[acc],
             "n_fa_fata_profiles": n, "n_FA": nf, "n_FATA": nt,
             "fa_subclass": sub,
             "profile_support": ("high" if n >= 5 else "medium" if n >= 3 else "minimal"),
             "FA_profiles": ";".join(f_hits), "FATA_profiles": ";".join(t_hits)}
        for p in MAND:
            h = hit.get(p, {}).get(acc)
            r[p.replace("T4SS_", "")] = ("%d:%d..%d:%s" % (idx[h[0]], st[h[0]], en[h[0]],
                "+" if strand[h[0]] == 1 else "-")) if h else "none"
        r["virb4_bit"] = hit["T4SS_virb4"][acc][1] if acc in hit.get("T4SS_virb4", {}) else ""
        r["relaxase_profiles"] = ";".join(p.replace("T4SS_", "") for p in MAND
                                          if p.startswith("T4SS_MOB") and acc in hit.get(p, {}))
        ex = [p for p in EEX if acc in hit.get(p, {})]
        r["exclusion_family_hit"] = ";".join(ex)
        r["slot_status"] = "no_slot_defined"
        rows.append(r)

    cols = list(rows[0].keys())
    dest = os.path.join(PROJ, "data", "release", "v1.1", "catalogue_MPF_FA_v1.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    n = len(rows)
    print("wrote %s (%d rows x %d cols)" % (dest, n, len(cols)))
    print("\n  fa_subclass : %s" % ", ".join("%s %d" % kv for kv in
          collections.Counter(r["fa_subclass"] for r in rows).most_common()))
    print("  profile_support: %s" % ", ".join("%s %d" % kv for kv in
          collections.Counter(r["profile_support"] for r in rows).most_common()))
    print("  n_fa_fata_profiles: median %d, max %d"
          % (sorted(r["n_fa_fata_profiles"] for r in rows)[n//2],
             max(r["n_fa_fata_profiles"] for r in rows)))
    print("  virb4 present : %d (100%% by construction -- it is in the entry criterion)"
          % sum(1 for r in rows if r["virb4"] != "none"))
    print("  relaxase present: %d (%.1f%%)"
          % (sum(1 for r in rows if r["relaxase_profiles"]),
             100.0*sum(1 for r in rows if r["relaxase_profiles"])/n))
    print("  any exclusion family: %d (%.1f%%)"
          % (sum(1 for r in rows if r["exclusion_family_hit"]),
             100.0*sum(1 for r in rows if r["exclusion_family_hit"])/n))
    print("  slot_status   : no_slot_defined on all %d rows" % n)


if __name__ == "__main__":
    main()
