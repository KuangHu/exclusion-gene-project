#!/usr/bin/env python3
"""Measure the slot as it was actually defined: the gene immediately 5' of VirB6.

The GO on the core hypothesis came from one observation: all four MPF_T controls
carry a small lipoprotein IMMEDIATELY 5' OF VirB6. Everything since has measured
something else -- the interval BETWEEN VirB5 and VirB6 -- which coincides with it
only when the operon is in canonical order.

Removing the VirB5-before-VirB6 assumption recovered 554 plasmids in which VirB6
precedes VirB5. They are real: in six independent Mash clusters the two proteins
score 34-43 (VirB5, 213-215 aa) and 120-137 (VirB6, 328-331 aa), both above GA,
with no other Pfam family hitting either. But in those plasmids the between-anchor
interval is NOT the 5'-of-VirB6 position, so calling them "empty slots" measured
the wrong window and dragged record occupancy from 82.4% to 73.8%.

So measure the defined thing directly. Two further gains:

  * it does not depend on VirB5 at all, and VirB5 is the weakest anchor in the
    set (21.4% failure vs ~2% for the others). Every occupancy figure so far has
    been conditional on VirB5 being detectable; this one is not.
  * canonical and inverted plasmids get the same treatment, with no order rule.

An earlier draft of this file predicted occupancy would go to ~100% and
concluded that v2 6.1 therefore kills the presence/absence direction. That was
wrong in a dangerous way. There is ALWAYS a gene 5' of VirB6 unless VirB6 sits at
a record boundary, so 100% would not be a measurement at all -- it would be an
identity forced by the definition, and retiring a whole research direction on it
would repeat the \r bug exactly: a non-biological cause producing a number that
lands precisely on a decision threshold.

Requiring VirB5 was what made "empty" definable in the first place. Dropping it
does not remove the conditionality; it removes the ability to measure.

So "empty" is redefined by the neighbour's IDENTITY, not its existence:

    upstream neighbour hits...          verdict
    a canonical T4SS component          slot_empty      (VirB1-11, VirD4, TraN)
    a known Eex family                  eex_occupied
    nothing                             candidate

This still measures something real, and the denominator improves: it is now
"VirB6 detectable" rather than "VirB5 detectable". The requirement is spread over
twelve component families instead of resting on PF07996, whose 21.4% failure rate
is an order of magnitude worse than VirB6/VirB8/VirB9 at 4.9%. It also handles the
inverted plasmids with no order rule at all -- whichever side VirB5 is on, the
question is only what sits 5' of VirB6.

Reported alongside, never merged, because they measure different things:

    slot_occupancy | VirB5 detectable          = the existing between-anchor figure
    neighbour_is_unassigned | VirB6 detectable = this one

Red line (v2 1.3): no exclusion-gene family is used to find or bound anything
here. The upstream gene is located by position alone. Family labels, if any, are
attached afterwards and never feed the selection.
"""
import argparse
import collections
import csv
import os
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJ, "scripts", "lib"))

HMM = "/global/scratch/users/kh36969/exclusion_gene/hmm"
CACHE = "/global/scratch/users/kh36969/exclusion_gene/cds_cache"

# Positive controls: the four MPF_T systems the hypothesis was built on. Each must
# come back with its known occupant. A miss here invalidates every other row.
CONTROLS = {
    "BN000925.1":  ("RP4",     "trbK",  69),
    "NC_001735.4": ("R751",    "trbK",  73),
    "U09868.1":    ("pKM101",  "eex",   76),
    "NC_028464.1": ("R388",    "candidate", 76),
}
FORBIDDEN = {"TIGR04359", "NF033894", "NF041429", "NF033891", "PF10624", "PF14729"}
ANCHORS = ["anchor_SLT", "anchor_TrbC", "anchor_VirB3", "anchor_CagE_TrbE_VirB",
           "anchor_T4SS", "TrbL", "anchor_VirB8", "anchor_CagX", "anchor_TrbI",
           "anchor_T2SSE", "anchor_T4SS-DNA_transf", "anchor_P_T4SS_TraN",
           "anchor_F_T4SS_TraN"]

# POST-HOC CLASSIFICATION ONLY (v2 1.3 red line). These label a gene that has
# ALREADY been selected by position. They never bound a window, never choose an
# anchor pair, and never decide which gene is examined -- doing any of that would
# make the exclusion families their own evidence.
EEX_FAMILIES = ["TIGR04359", "NF033894", "NF041429", "NF033891",
                "TraS", "DUF4467"]


def load_cache(shard, nshards):
    """accession -> list of CDS dicts in coordinate order, from the frozen cache."""
    import hashlib
    _h = lambda x: int(hashlib.md5(x.encode()).hexdigest(), 16)
    genes = collections.defaultdict(list)
    for fn in sorted(os.listdir(CACHE)):
        if not fn.endswith(".faa"):
            continue
        name = None
        for line in open(os.path.join(CACHE, fn)):
            if line[0] == ">":
                name = line[1:].rstrip("\n")
            else:
                acc, i, s, e, st, ln = name.split("|")
                if acc in CONTROLS:
                    continue          # loaded from the seed record in shard 0
                if nshards > 1 and _h(acc) % nshards != shard:
                    continue
                genes[acc].append({"i": int(i), "start": int(s), "end": int(e),
                                   "strand": int(st), "aa_len": int(ln),
                                   "aa": line.rstrip("\n")})
    for v in genes.values():
        v.sort(key=lambda c: c["i"])
    return genes


CONTROL_GB = {
    "BN000925.1":  "RP4__BN000925.1.gb",
    "NC_001735.4": "R751__NC_001735.4.gb",
    "U09868.1":    "pKM101__U09868.1.gb",
    "NC_028464.1": "R388__NC_028464.1.gb",
}


def load_controls():
    """Call CDS on the four control records directly, bypassing the cache.

    pKM101 has no complete plasmid sequence -- U09868.1 is a tra-region fragment
    -- so it is absent from PLSDB, absent from the step-1a VirB4+ set, and would
    therefore be absent from the cache. Taking the baseline from whatever controls
    happened to survive the gating is exactly how a benchmark ends up scoring
    itself on the easy cases. All four are loaded here regardless.
    """
    import orf_caller
    from Bio import SeqIO
    out = {}
    for acc, fn in CONTROL_GB.items():
        p = os.path.join(PROJ, "data", "seed", "genbank", fn)
        if not os.path.exists(p):
            sys.stderr.write("control record missing: %s\n" % fn)
            continue
        rec = next(SeqIO.parse(p, "genbank"))
        out[acc] = [{"i": i, "start": c["start"], "end": c["end"],
                     "strand": c["strand"], "aa_len": c["aa_len"], "aa": c["aa"]}
                    for i, c in enumerate(orf_caller.call(str(rec.seq).upper()))]
    return out


def main():
    from assertions import require_compute_node
    require_compute_node()          # A11: no heavy scans on a login node
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--cpus", type=int, default=8)
    a = ap.parse_args()
    for x in ANCHORS:
        if x in FORBIDDEN:
            raise SystemExit("RED LINE: %s is an exclusion family" % x)
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()

    genes = load_cache(a.shard, a.nshards)
    if not genes:
        raise SystemExit("CDS cache empty or shard empty -- run 93_cds_cache.py")
    if a.shard == 0:
        genes.update(load_controls())
    sys.stderr.write("accessions in shard: %d\n" % len(genes))

    flat, owner = [], []
    for acc, cs in genes.items():
        for c in cs:
            owner.append((acc, c["i"]))
            flat.append(pyhmmer.easel.TextSequence(
                name=str(len(flat)).encode(), sequence=c["aa"]).digitize(alpha))

    def scan(fn, ev=None, minaa=0):
        with pyhmmer.plan7.HMMFile(os.path.join(HMM, fn + ".hmm")) as fh:
            model = next(iter(fh))
        kw = {"E": ev} if ev else {"bit_cutoffs": "gathering"}
        out = {}
        for top in pyhmmer.hmmsearch([model], flat, cpus=a.cpus, **kw):
            for h in top:
                nm = h.name
                j = int(nm.decode() if isinstance(nm, bytes) else nm)
                acc, i = owner[j]
                if minaa and genes[acc][i]["aa_len"] <= minaa:
                    continue
                if (acc, i) not in out or h.score > out[(acc, i)]:
                    out[(acc, i)] = h.score
        return out

    virb6 = scan("TrbL")
    virb5 = scan("anchor_T4SS", ev=1e-5, minaa=150)
    # family labels for the upstream gene -- attached AFTER selection, never before
    labels = collections.defaultdict(list)
    for fam in ANCHORS:
        for k in scan(fam):
            labels[k].append(fam.replace("anchor_", ""))
    eex = collections.defaultdict(list)
    for fam in EEX_FAMILIES:
        for k in scan(fam):
            eex[k].append(fam)

    rows = []
    for acc, cs in genes.items():
        sixes = sorted([i for (b, i) in virb6 if b == acc])
        if not sixes:
            continue
        for i6 in sixes:
            g6 = cs[i6]
            up = i6 - 1 if g6["strand"] == 1 else i6 + 1
            if up < 0 or up >= len(cs):
                continue
            u = cs[up]
            if labels.get((acc, up)) or (acc, up) in virb5:
                verdict = "slot_empty"
                cls = ("is_VirB5" if (acc, up) in virb5
                       else "is_" + labels[(acc, up)][0])
            elif eex.get((acc, up)):
                verdict = "eex_occupied"
                cls = "is_" + eex[(acc, up)][0]
            else:
                verdict = "candidate"
                cls = "unlabelled"
            if u["strand"] != g6["strand"]:
                cls += "|opposite_strand"
            gapbp = (g6["start"] - u["end"] - 1 if g6["strand"] == 1
                     else u["start"] - g6["end"] - 1)
            rows.append({
                "accession": acc, "virb6_index": i6,
                "virb6_score": round(virb6[(acc, i6)], 1),
                "virb6_aa": g6["aa_len"], "virb6_strand": g6["strand"],
                "virb6_coords": "%d..%d" % (g6["start"], g6["end"]),
                "upstream_aa": u["aa_len"],
                "upstream_coords": "%d..%d" % (u["start"], u["end"]),
                "intergenic_bp": gapbp,
                "slot_verdict": verdict,
                "upstream_class": cls,
                "upstream_eex_families": ";".join(eex.get((acc, up), [])),
                "upstream_families": ";".join(labels.get((acc, up), [])),
                "has_virb5_anywhere": int(any(b == acc for (b, i) in virb5)),
                "upstream_aa_seq": u["aa"],
            })

    out = os.path.join(PROJ, "data", "positional")
    os.makedirs(out, exist_ok=True)
    sfx = "" if a.nshards == 1 else "_%03d" % a.shard
    dest = os.path.join(out, "upstream_of_virb6%s.tsv" % sfx)
    cols = ["accession", "virb6_index", "virb6_score", "virb6_aa", "virb6_strand",
            "virb6_coords", "upstream_aa", "upstream_coords", "intergenic_bp",
            "slot_verdict", "upstream_class", "upstream_families",
            "upstream_eex_families", "has_virb5_anywhere",
            "upstream_aa_seq"]
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    # --- baseline: the four controls must each return their known occupant ----
    print("=" * 74)
    print("BASELINE -- the four MPF_T controls. Expected to score perfectly.")
    print("=" * 74)
    seen = 0
    for acc, (name, gene, aa) in CONTROLS.items():
        mine = [r for r in rows if r["accession"] == acc]
        if not mine:
            continue
        seen += 1
        ok = any(abs(r["upstream_aa"] - aa) <= 3 for r in mine)
        for r in mine:
            print("  %-8s %-10s expect %3d aa   got %3d aa (%s)  %s"
                  % (name, gene, aa, r["upstream_aa"], r["upstream_class"],
                     "OK" if abs(r["upstream_aa"] - aa) <= 3 else "MISMATCH"))
        if not ok:
            print("  *** %s BASELINE FAILED -- stop; do not read the "
                  "distribution below ***" % name)
    if seen < len(CONTROLS) and a.shard == 0:
        print("  *** only %d/%d controls returned a row -- the baseline is "
              "incomplete and nothing below may be read ***"
              % (seen, len(CONTROLS)))
    if not seen:
        print("  (no controls in this shard)")

    print("\n  rows: %d over %d accessions" % (len(rows), len(genes)))
    vv = collections.Counter(r["slot_verdict"] for r in rows)
    print("\n  neighbour_is_unassigned | VirB6 detectable")
    for k in ("slot_empty", "eex_occupied", "candidate"):
        print("    %-16s %6d  (%.1f%%)" % (k, vv.get(k, 0),
                                           100.0 * vv.get(k, 0) / len(rows)))
    print("    -> slot NOT empty (eex + candidate): %.1f%%"
          % (100.0 * (vv.get("eex_occupied", 0) + vv.get("candidate", 0)) / len(rows)))
    cc = collections.Counter(r["upstream_class"] for r in rows)
    print("\n  WHAT SITS 5' OF VirB6")
    for k, v in cc.most_common(15):
        print("    %-22s %6d  (%.1f%%)" % (k, v, 100.0 * v / len(rows)))

    unl = [r["upstream_aa"] for r in rows if r["slot_verdict"] == "candidate"]
    if unl:
        s = sorted(unl)
        import statistics
        print("\n  LENGTH of the CANDIDATE upstream gene (n=%d)" % len(s))
        h = collections.Counter(min(16, l // 25) for l in s)
        mx = max(h.values())
        for b in range(17):
            lab = "%d-%d" % (b * 25, b * 25 + 24) if b < 16 else ">=400"
            star = "   <- controls 69-76 aa" if b == 2 else ""
            print("    %-10s %6d  %s%s" % (lab, h.get(b, 0),
                                           "#" * int(50 * h.get(b, 0) / mx), star))
        print("    median %d  Q1 %d  Q3 %d" % (statistics.median(s),
                                               s[len(s)//4], s[3*len(s)//4]))
        band = sum(1 for l in s if 60 <= l <= 100)
        print("    in the 60-100 aa control band: %d (%.1f%%)"
              % (band, 100.0 * band / len(s)))
    print("\nwrote %s" % dest)


if __name__ == "__main__":
    main()
