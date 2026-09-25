#!/usr/bin/env python3
"""Stage 1b -- triage accession candidates into actionable classes.

Key insight: for most seed elements the top two hits are the SAME SEQUENCE
surfaced twice, once as RefSeq (NC_/NZ_) and once as the original GenBank
submission. That is the dual-accession pair the schema wants, not an ambiguity.
Identical length + one of each namespace => RESOLVED_PAIR.

Biological filter: a conjugative element carries a full T4SS and cannot be 2 kb.
Candidates below min_expected_bp are fragments (tra-region submissions, single
resistance genes) and are reported as such -- that is the pKM101 case, where no
complete-plasmid record exists and the element must be assembled from parts.
"""
import csv, collections, os, re, sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFSEQ = re.compile(r"^(NC_|NZ_|NG_)")
TPA    = re.compile(r"^(BR|BK)\d")          # third-party annotation, not a submission

def load(path, key=None):
    with open(path) as f:
        rows = [l.rstrip("\n").split("\t") for l in f
                if not l.startswith("#") and l.strip()]
    hdr, data = rows[0], rows[1:]
    d = [dict(zip(hdr, r)) for r in data]
    return {r[key]: r for r in d} if key else d

def main():
    els   = load(os.path.join(PROJ, "config", "seed_elements.tsv"), key="element_id")
    cands = load(os.path.join(PROJ, "data", "seed", "accession_candidates.tsv"))
    by = collections.OrderedDict()
    for c in cands:
        by.setdefault(c["element_id"], []).append(c)

    out_rows = []
    for eid, cs in by.items():
        el = els[eid]
        minbp = int(el["min_expected_bp"])
        exp_org = el["expected_organism"]

        real = [c for c in cs if c["source"] != "NONE"]
        # organism assertion is a HARD filter when asserted
        if exp_org not in ("", "NA"):
            real = [c for c in real if c["organism_ok"] in ("YES", "unasserted")]
        full = [c for c in real
                if c["length"].isdigit() and int(c["length"]) >= minbp
                and not c["is_derivative"]]

        if not real:
            verdict, refseq, gb, note = "NO_CANDIDATE", "", "", "no hit survives organism filter"
        elif not full:
            best = max((int(c['length']) for c in real if c['length'].isdigit()), default=0)
            verdict, refseq, gb = "FRAGMENTS_ONLY", "", ""
            note = f"largest hit {best}bp < {minbp}bp minimum; assemble from parts"
        else:
            groups = collections.defaultdict(list)
            for c in full:
                groups[c["length"]].append(c)
            # prefer the length group with the most independent records
            best_len = max(groups, key=lambda L: (len(groups[L]), int(L)))
            grp = groups[best_len]
            rs = [c["accver"] for c in grp if REFSEQ.match(c["accver"])]
            gbs = [c["accver"] for c in grp
                   if not REFSEQ.match(c["accver"]) and not TPA.match(c["accver"])]
            refseq = rs[0] if rs else ""
            gb     = gbs[0] if gbs else ""
            others = [L for L in groups if L != best_len]
            if refseq and gb:
                verdict = "RESOLVED_PAIR"
            elif refseq or gb:
                verdict = "SINGLE_ACCESSION"
            else:
                verdict = "TPA_ONLY"
            note = f"{best_len}bp"
            if others:
                note += f"; {len(others)} other length group(s) rejected: {','.join(sorted(others)[:3])}"
            if len(rs) > 1:
                note += f"; MULTIPLE RefSeq at this length: {','.join(rs)}"
        out_rows.append([eid, el["role"], el["resolution_tier"], verdict,
                         refseq, gb, note])

    outp = os.path.join(PROJ, "data", "seed", "accession_triage.tsv")
    with open(outp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["element_id","role","tier","verdict",
                    "refseq_accession","genbank_accession","note"])
        w.writerows(out_rows)

    order = ["RESOLVED_PAIR","SINGLE_ACCESSION","TPA_ONLY","FRAGMENTS_ONLY","NO_CANDIDATE"]
    for v in order:
        hits = [r for r in out_rows if r[3] == v]
        if not hits: continue
        print(f"\n=== {v}  ({len(hits)}) ===")
        for r in hits:
            print(f"  {r[0]:<12} {r[4]:<16} {r[5]:<14} {r[6][:78]}")
    print(f"\nwrote {outp}")

if __name__ == "__main__":
    main()
