#!/usr/bin/env python3
"""Regression check: re-extract every target and compare against the frozen SHA.

This is the payoff of the inverted source-of-truth design. The extractor is not
an oracle; the frozen (coords + aa_sha256) is. Re-running this after any NCBI
refresh answers one question: does the current annotation still yield the same
protein?

Fails LOUD on mismatch. A silent pass is the only acceptable success.

Coordinates are validated ONLY against accession_authoritative. Never against
accession_related -- RP4's two records differ by a 2 bp frame offset while
carrying identical proteins, so a related record would silently resolve to the
wrong span.
"""
import csv, glob, os, sys
from Bio import SeqIO
sys.path.insert(0,os.path.join(os.path.dirname(__file__),"lib"))
import ncbi

PROJ=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GB=os.path.join(PROJ,"data","seed","genbank")

def load_targets():
    rows=[]
    with open(os.path.join(PROJ,"config","extraction_targets.tsv")) as f:
        for line in f:
            if line.startswith("#"): continue
            rows.append(line.rstrip("\n").split("\t"))
    return [dict(zip(rows[0],r)) for r in rows[1:]]

HOSTS="/global/scratch/users/kh36969/exclusion_gene/hosts"

def find_record(acc):
    hits=glob.glob(os.path.join(GB,f"*__{acc}.gb"))
    if hits: return hits[0]
    # CHROMOSOMAL hosts live on scratch (bulk). This is a DURABILITY GAP: the
    # element's authoritative record is then outside the repo, contrary to S7's
    # frozen-accession-set requirement. See docs/durability_gap.md.
    hits=glob.glob(os.path.join(HOSTS,f"*__{acc}.gb"))
    return hits[0] if hits else None

def protein_sha(acc, gene, coords=""):
    """Pull the translation from the AUTHORITATIVE record.

    COORDINATES ARE THE ANCHOR, not the gene name. The freeze stores
    (coords + aa_sha256); re-extraction asks whether the protein at those
    coordinates still hashes the same. Looking up by name instead would miss
    ICEBs1 yddJ, which is annotated only as DUF4467 with no yddJ in the product.
    Name matching remains as a fallback for rows whose coords are not yet frozen.
    """
    p=find_record(acc)
    if not p: return None,"RECORD_MISSING"
    rec=next(SeqIO.parse(p,"genbank"))
    import re
    if coords and ":" in coords:
        want=coords.split(":")[1]
        try: lo,hi=(int(x) for x in want.split(".."))
        except ValueError: lo=hi=None
        if lo:
            for f in rec.features:
                if f.type!="CDS": continue
                if int(f.location.start)+1==lo and int(f.location.end)==hi:
                    tr=(f.qualifiers.get("translation") or [""])[0]
                    if tr: return ncbi.aa_sha256(tr),f"{lo}..{hi} (by coords)"
            return None,f"NO_CDS_AT_{lo}..{hi}"
    for f in rec.features:
        if f.type!="CDS": continue
        g=((f.qualifiers.get("gene") or [""])[0]).strip()
        prod=(f.qualifiers.get("product") or [""])[0]
        base=gene.split("-")[0].split("_")[0]
        if g.lower() in (gene.lower(),base.lower()) or \
           re.search(rf"\b{re.escape(base)}\b",prod,re.I) or \
           (gene=="eex" and g.lower().startswith("eex")):
            tr=(f.qualifiers.get("translation") or [""])[0]
            if tr: return ncbi.aa_sha256(tr),f"{int(f.location.start)+1}..{int(f.location.end)}"
    return None,"GENE_NOT_FOUND"

def main():
    tg=load_targets()
    checked=ok=fail=skipped=0
    print(f"{'element':<12} {'gene':<11} {'result':<10} detail")
    print("-"*84)
    for t in tg:
        auth=t["accession_authoritative"]
        for side in ("exclusion","partner"):
            gene=t[f"{side}_gene"]; exp=t[f"{side}_expected_sha256"]
            if not exp:
                skipped+=1
                print(f"{t['element_id']:<12} {gene:<11} {'no-baseline':<10} not yet frozen")
                continue
            checked+=1
            got,detail=protein_sha(auth,gene,t.get(f"{side}_coords",""))
            if got==exp:
                ok+=1; print(f"{t['element_id']:<12} {gene:<11} {'OK':<10} {auth}:{detail}")
            else:
                fail+=1
                print(f"{t['element_id']:<12} {gene:<11} {'MISMATCH':<10} {detail}")
                print(f"{'':<12} {'':<11} {'':<10} expected {exp[:16]}  got {(got or 'None')[:16]}")
    print("-"*84)
    print(f"INTERNAL CONSISTENCY: checked {checked}  ok {ok}  MISMATCH {fail}  no-baseline {skipped}")

    # COVERAGE is a separate question from consistency. A green consistency run
    # says nothing about how much of the seed is actually populated -- 16/16 PASS
    # with three systems unrecovered would otherwise read as "Level 1 complete".
    import csv as _csv
    refs=os.path.join(PROJ,"data","level1_literature","references.csv")
    n_def=0
    if os.path.exists(refs):
        elems=set()
        for r in _csv.DictReader(open(refs)):
            if r["ref_role"]=="defining_sequence":
                for e in r["element_id"].split(";"): elems.add(e.strip())
        n_def=len(elems)
    got=sum(1 for t in tg if t["exclusion_expected_sha256"])
    print(f"COVERAGE            : exclusion gene recovered for {got}/{len(tg)} elements")
    if n_def: print(f"                      defining_sequence papers cover {n_def} elements")
    missing=[t["element_id"] for t in tg if not t["exclusion_expected_sha256"]]
    if missing: print(f"                      NOT RECOVERED: {', '.join(missing)}")
    if fail:
        print("\nFAIL: a frozen protein no longer matches its record. Do NOT proceed.")
        return 1
    print("\nPASS" if checked else "\nNo baselines to check yet.")
    return 0

if __name__=="__main__": sys.exit(main())
