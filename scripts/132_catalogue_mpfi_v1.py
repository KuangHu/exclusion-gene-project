#!/usr/bin/env python3
"""MPF_I (IncI) catalogue.

Entry criterion `T4SS_I_traU` -- 3,732 plasmids, 96.6% ExcA recall, 0.46%
cross-class. Slot is TraY+1, validated 86-98% across strata.

SCHEMA DECISIONS (as specified):

* `slot_maxkd` is a COLUMN, not a filter. `slot_admission` is `undetermined`
  throughout. Hydrophobicity was measured at 2.3x over a 40.5% background --
  too weak to admit on, and MPF_T's lipobox (56.1% vs 2.8%) does not transfer.
  Carrying the number without a verdict is the honest form.
* excB is reported in FOUR columns, never collapsed to one call:
      excA_hit            NF033891 at GA
      excB_inframe_orfs   ALL in-frame reinitiation candidates, not just 147 aa
      excB_family_hit     any exclusion family on a translated candidate
      excB_call_mode      how excB was called, or why it was not
  excB is an in-frame reinitiation ORF inside excA, so it shares excA's stop.
  Reporting only the 147 aa seed length would be fitting the answer to the seeds.
* `slot_ready__traY` follows the v1.1 definition: anchor present AND slot
  occupant resolvable AND architecture callable. Same three-part shape as
  `slot_ready__virb5_virb6`.

RED LINES OBSERVED: exclusion families (NF033891 ExcA included) are POST-HOC
classification only -- they never enter the entry criterion, the anchors, the
architecture or the slot definition. `assign_method = positional` never feeds
`slot_status`.

NOT an operon annotation.
"""
import argparse, collections, csv, os, sys
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))
from assertions import require_compute_node

CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache_full"
HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CJ = "/global/scratch/users/kh36969/funcannot_dbs/macsy_models/CONJScan/profiles"
PLSDB_FA = "/global/scratch/users/kh36969/plsdb/sequences.fasta"
SURVEY = os.path.join(PROJ, "data", "anchors", "mpfi_architecture_survey.tsv")
ENTRY = "T4SS_I_traU"
SLOT_ANCHOR = "T4SS_I_traY"
EEX_FAMS = ["NF033891", "NF033894", "NF041429", "TIGR04359", "PF10624", "PF14729"]
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
import re as _re
LIPO = _re.compile(r"[LVI][ASTVIG][GASN]C")


def maxkd(s, w=19):
    if len(s) < w:
        return round(sum(KD.get(c, 0) for c in s) / max(1, len(s)), 2)
    return round(max(sum(KD.get(c, 0) for c in s[i:i+w]) / w
                     for i in range(len(s)-w+1)), 2)


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
    print("cache: %d proteins" % len(seqs), flush=True)

    def scan(path):
        with pyhmmer.plan7.HMMFile(path) as fh: m = next(iter(fh))
        best = {}
        for top in pyhmmer.hmmsearch([m], seqs, cpus=a.cpus, bit_cutoffs="gathering"):
            for h in top:
                i = int(dec(h.name)); acc = owner[i]
                if acc not in best or h.score > best[acc][1]: best[acc] = (i, h.score)
        return best

    fams = sorted(f[:-4] for f in os.listdir(CJ)
                  if f.startswith("T4SS_I_") and f.endswith(".hmm"))
    hit = {}
    for f in fams:
        hit[f] = scan(os.path.join(CJ, f + ".hmm"))
        print("  %-16s %5d" % (f, len(hit[f])), flush=True)
    eex = {}
    for f in EEX_FAMS:
        p = os.path.join(HMM, f + ".hmm")
        if os.path.exists(p):
            eex[f] = scan(p); print("  eex %-12s %5d" % (f, len(eex[f])), flush=True)

    admitted = sorted(hit[ENTRY])
    print("\nadmitted by %s: %d" % (ENTRY, len(admitted)), flush=True)

    surv = {r["accession"]: r for r in csv.DictReader(open(SURVEY), delimiter="\t")}
    oc = collections.Counter(r["order_canonical"] for r in surv.values())
    rank = {o: i + 1 for i, (o, _) in enumerate(oc.most_common())}

    def stratum(acc):
        r = surv.get(acc)
        if not r: return "uncallable"
        if r["layout"] != "contiguous": return r["layout"]
        k = rank.get(r["order_canonical"], 999)
        return "order#%d" % k if k <= 3 else "order#4+"

    # --- excB: in-frame reinitiation ORFs inside excA (needs nucleotides) ----
    need = {acc: eex["NF033891"][acc][0] for acc in admitted
            if "NF033891" in eex and acc in eex["NF033891"]}
    print("plasmids with ExcA (NF033891) among admitted: %d" % len(need), flush=True)
    excb = {}
    if need:
        from Bio import SeqIO
        from Bio.Seq import Seq
        seen = 0
        for rec in SeqIO.parse(PLSDB_FA, "fasta"):
            if rec.id not in need: continue
            seen += 1
            i = need[rec.id]
            sub = str(rec.seq).upper()[st[i]-1:en[i]]
            if strand[i] != 1: sub = str(Seq(sub).reverse_complement())
            cands = []
            for off in range(3, len(sub) - 3, 3):
                if sub[off:off+3] in ("ATG", "GTG", "TTG"):
                    naa = (len(sub) - off) // 3 - 1
                    if naa >= 30: cands.append((off, naa))
            excb[rec.id] = cands
            if seen == len(need): break
        print("  scanned %d of %d genomes for in-frame excB candidates"
              % (seen, len(need)), flush=True)
        if seen != len(need):
            miss = sorted(set(need) - set(excb))
            print("  NOT RETRIEVED from the nucleotide FASTA: %d -- %s"
                  % (len(miss), ", ".join(miss[:5])), flush=True)
            print("  these are reported as excB_call_mode=genome_not_retrieved,")
            print("  NEVER as 'no candidate': not looking is not a negative result.",
                  flush=True)

    rows = []
    for acc in admitted:
        sv = surv.get(acc, {})
        s = stratum(acc)
        ty = hit[SLOT_ANCHOR].get(acc)
        r = {"accession": acc, "mpf_class_call": "MPF_I",
             "entry_criterion": ENTRY,
             "n_anchors": sv.get("n_anchors", ""), "layout": sv.get("layout", "uncallable"),
             "order_canonical": sv.get("order_canonical", ""), "stratum": s,
             "anchor_span_bp": sv.get("anchor_span_bp", ""),
             "max_internal_gap_bp": sv.get("max_internal_gap_bp", "")}
        for f in fams:
            h = hit[f].get(acc)
            r[f.replace("T4SS_I_", "")] = ("%d:%d..%d:%s" % (idx[h[0]], st[h[0]], en[h[0]],
                "+" if strand[h[0]] == 1 else "-")) if h else "none"
        # slot (TraY+1) is filled in the second pass, once the per-accession
        # gene index exists -- see below
        r["traY_bit"] = round(ty[1], 1) if ty else ""
        rows.append((r, ty, s, acc))

    # slot occupant needs a per-accession gene index; build once
    byacc = collections.defaultdict(dict)
    for i in range(len(owner)):
        byacc[owner[i]][idx[i]] = i
    out = []
    for r, ty, s, acc in rows:
        occ_i = None
        if ty:
            t = 1 if strand[ty[0]] == 1 else -1
            occ_i = byacc[acc].get(idx[ty[0]] + t)
        if occ_i is not None:
            pep = aas[occ_i]
            r.update({"slot_resolved": 1,
                      "slot_coords": "%d..%d" % (st[occ_i], en[occ_i]),
                      "slot_strand": "+" if strand[occ_i] == 1 else "-",
                      "slot_aa_len": len(pep),
                      "slot_maxkd": maxkd(pep),
                      "slot_lipobox": int(bool(LIPO.search(pep[:40]))),
                      "slot_nterm30": pep[:30]})
        else:
            r.update({"slot_resolved": 0, "slot_coords": "", "slot_strand": "",
                      "slot_aa_len": "", "slot_maxkd": "", "slot_lipobox": "",
                      "slot_nterm30": ""})
        r["slot_admission"] = "undetermined"
        r["slot_ready__traY"] = int(bool(ty) and r["slot_resolved"] == 1
                                    and s != "uncallable")
        # --- exclusion block, POST-HOC ------------------------------------
        ex = eex.get("NF033891", {}).get(acc)
        r["excA_hit"] = int(bool(ex))
        r["excA_aa"] = len(aas[ex[0]]) if ex else ""
        r["excA_bit"] = round(ex[1], 1) if ex else ""
        if ex and ty:
            t = 1 if strand[ty[0]] == 1 else -1
            r["excA_offset_from_traY"] = (idx[ex[0]] - idx[ty[0]]) * t
        else:
            r["excA_offset_from_traY"] = ""
        scanned = acc in excb
        c = excb.get(acc, [])
        r["excB_scanned"] = int(scanned)
        r["excB_inframe_orfs"] = ";".join(str(n) for _, n in c) if c else ""
        r["excB_n_candidates"] = len(c) if scanned else ""
        r["excB_has_147aa"] = int(any(n == 147 for _, n in c)) if scanned else ""
        fam = [f for f in EEX_FAMS if f != "NF033891" and acc in eex.get(f, {})]
        r["excB_family_hit"] = ";".join(fam)
        if not ex:
            r["excB_call_mode"] = "no_excA"
        elif not scanned:
            # ExcA present but the genome was not in the nucleotide FASTA.
            # This is ABSENCE OF EVIDENCE, not evidence of absence, and must not
            # share a label with a genome that was searched and came back empty.
            r["excB_call_mode"] = "genome_not_retrieved"
        elif c:
            r["excB_call_mode"] = "sixframe_inframe"
        elif fam:
            r["excB_call_mode"] = "family_only"
        else:
            r["excB_call_mode"] = "excA_but_no_candidate"
        out.append(r)

    cols = list(out[0].keys())
    dest = os.path.join(PROJ, "data", "release", "v1.1", "catalogue_MPF_I_v1.tsv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           restval="")
        w.writeheader(); w.writerows(out)
    n = len(out)
    print("\nwrote %s (%d rows x %d cols)" % (dest, n, len(cols)))
    sc = collections.Counter(r["stratum"] for r in out)
    print("\n  stratum: %s" % ", ".join("%s %d" % kv for kv in sc.most_common()))
    print("  slot_resolved     %d (%.1f%%)"
          % (sum(r["slot_resolved"] for r in out), 100.0*sum(r["slot_resolved"] for r in out)/n))
    print("  slot_ready__traY  %d (%.1f%%)"
          % (sum(r["slot_ready__traY"] for r in out), 100.0*sum(r["slot_ready__traY"] for r in out)/n))
    print("  excA_hit          %d (%.1f%%)"
          % (sum(r["excA_hit"] for r in out), 100.0*sum(r["excA_hit"] for r in out)/n))
    print("  excB call mode: %s"
          % ", ".join("%s %d" % kv for kv in
                      collections.Counter(r["excB_call_mode"] for r in out).most_common()))
    sc = [r for r in out if r["excB_scanned"] == 1]
    print("  excB_has_147aa    %d of %d SCANNED (unscanned rows carry '', never 0)"
          % (sum(r["excB_has_147aa"] for r in sc), len(sc)))
    kd = sorted(r["slot_maxkd"] for r in out if r["slot_maxkd"] != "")
    if kd:
        print("  slot_maxkd: median %.2f  Q1 %.2f  Q3 %.2f  (COLUMN ONLY, no admission)"
              % (kd[len(kd)//2], kd[len(kd)//4], kd[3*len(kd)//4]))


if __name__ == "__main__":
    main()
