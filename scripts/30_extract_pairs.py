#!/usr/bin/env python3
"""Stage 3 -- ONE-TIME candidate extraction. Not the source of truth.

Route order per target gene:  /gene (incl. synonyms) -> /product regex -> miss.
Every hit carries its route, coords, protein_id, locus_tag and aa_sha256 so the
human pass can verify it and FREEZE it. After the freeze, 50_regression_check.py
re-runs this and compares SHAs; a RefSeq re-annotation that renames a product
then fails loudly instead of silently dropping the gene.

Negative arm: for role=negative elements every exclusion-side pattern is run and
ANY hit is a FALSE POSITIVE. A pipeline that finds all positives and also
hallucinates hits in RSF1010/pAM373 has a precision problem invisible until S4.
"""
import csv, os, re, sys
from Bio import SeqIO
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
import ncbi

PROJ   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GB_DIR = os.path.join(PROJ,"data","seed","genbank")
OUT    = os.path.join(PROJ,"data","seed","extraction_candidates.tsv")

def load(p):
    with open(p) as f:
        rows=[l.rstrip("\n").split("\t") for l in f if not l.startswith("#") and l.strip()]
    return [dict(zip(rows[0],r)) for r in rows[1:]]

def main():
    els   = {e["element_id"]:e for e in load(os.path.join(PROJ,"config","extraction_targets.tsv"))}
    pins  = {p["element_id"]:p for p in load(os.path.join(PROJ,"config","pinned_accessions.tsv"))}
    tgts  = {t["gene_canonical"]:t for t in load(os.path.join(PROJ,"config","gene_targets.tsv"))}
    syns  = load(os.path.join(PROJ,"config","synonyms.tsv"))
    man   = load(os.path.join(PROJ,"data","seed","manifest.tsv"))

    # historical -> canonical, for the /gene route
    syn_map={}
    for s in syns:
        if s["relation"]=="synonym" and s["canonical"]!="NA":
            syn_map[s["historical_name"].lower()]=s["canonical"]
    deprecated={s["historical_name"].lower() for s in syns if s["relation"]=="deprecated"}

    excl_genes=[g for g,t in tgts.items() if t["side"]=="exclusion"]

    # records per element, primary first
    recs={}
    order={"primary":0,"aux":1,"crosscheck":2}
    for m in man:
        if m["status"]!="OK" or not m["path"]: continue
        p=os.path.join(PROJ,m["path"]) if not m["path"].startswith("/") else m["path"]
        if not os.path.exists(p): continue
        recs.setdefault(m["element_id"],[]).append((order.get(m["record_role"],9),m["record_role"],m["accession"],p))
    for k in recs: recs[k].sort()

    # an element row may consolidate several sequenced representatives; the
    # manifest files records under the REPRESENTATIVE name, not the row id
    RECORD_ALIAS={"IncC":["pVCR94","R16a"],"SXT":["SXT","R391"]}
    for row_id,alias in RECORD_ALIAS.items():
        if row_id not in recs:
            merged=[]
            for a in alias: merged+=recs.get(a,[])
            if merged: recs[row_id]=sorted(merged)

    rows=[]
    undefined=set()
    for eid,el in els.items():
        if eid not in recs:
            # NEVER skip silently: a missing record set must appear in the output
            # or it vanishes from every recall count
            rows.append([eid,el.get("role_exclusion","positive"),"*","*","","","NO_RECORDS",
                         "","","","","","","","NO_RECORDS_FOR_ELEMENT"])
            continue
        pin=pins.get(eid,{})
        inc=el.get("inc_group","NA")
        role=el.get("role_exclusion","positive")
        if role=="negative":
            wanted=[(g,"exclusion") for g in excl_genes]
        else:
            wanted=[(g,"exclusion") for g in el["exclusion_gene"].split(",") if g and g!="NA"]
            wanted+=[(g,"partner")  for g in el["partner_gene"].split(",")  if g and g not in("NA","self")]
        found=set()
        for _,rrole,acc,path in recs[eid]:
            try: rec=next(SeqIO.parse(path,"genbank"))
            except Exception: continue
            for f in rec.features:
                if f.type!="CDS": continue
                gene=((f.qualifiers.get("gene") or [""])[0]).strip()
                prod=(f.qualifiers.get("product") or [""])[0]
                if gene.lower() in deprecated: continue
                tr=(f.qualifiers.get("translation") or [""])[0]
                if not tr: continue
                canon_gene=syn_map.get(gene.lower(),gene)
                for g,side in wanted:
                    if (g,acc) in found: continue
                    t=tgts.get(g)
                    if not t:
                        # an undefined gene used to vanish here, surfacing as
                        # "no-baseline" rather than as a real problem
                        undefined.add((eid,g)); continue
                    # sfx and friends may be restricted to specific Inc groups
                    restr=t.get("restrict_inc_group","NA")
                    if restr and restr!="NA" and inc not in restr.split(","):
                        continue
                    route=None
                    if canon_gene.lower()==g.lower(): route="gene"
                    elif gene.lower()==g.lower():     route="gene"
                    # name_regex is tested against BOTH the gene qualifier and the
                    # product. Testing only the product made locus-tag aliases
                    # (R27's R0018/R0128) unmatchable, since their product is just
                    # "hypothetical protein".
                    elif re.search(t["name_regex"],gene,re.I): route="gene_regex"
                    elif re.search(t["name_regex"],prod,re.I): route="product"
                    if not route: continue
                    # a phenotype phrase never assigns identity -- it only flags
                    ph=t.get("phenotype_regex","NA")
                    phen_only = (route=="product" and ph not in("","NA")
                                 and re.search(ph,prod,re.I)
                                 and not re.search(t["name_regex"],prod,re.I))
                    lo,hi=int(f.location.start)+1,int(f.location.end)
                    amin,amax=t["aa_min"],t["aa_max"]
                    flag=""
                    if amin!="NA" and amax!="NA" and not (int(amin)<=len(tr)<=int(amax)):
                        flag="LENGTH_OUT_OF_WINDOW"
                    if phen_only: flag=(flag+";" if flag else "")+"PHENOTYPE_ONLY_REVIEW"
                    rows.append([eid,role,g,side,rrole,acc,route,
                                 f"{acc}:{lo}..{hi}",("+" if f.location.strand==1 else "-"),
                                 (f.qualifiers.get("protein_id") or [""])[0],
                                 (f.qualifiers.get("locus_tag") or [""])[0],
                                 len(tr),ncbi.aa_sha256(tr),prod[:70],flag])
                    found.add((g,acc))
        for g,side in wanted:
            if not any(r[0]==eid and r[2]==g for r in rows) and role!="negative":
                rows.append([eid,role,g,side,"","","MISS","","","","","","","","NOT_FOUND"])

    cols=["element_id","role","gene","side","record_role","accession","route",
          "coords","strand","protein_id","locus_tag","aa_len","aa_sha256","product","flag"]
    with open(OUT,"w",newline="") as fh:
        w=csv.writer(fh,delimiter="\t", lineterminator="\n"); w.writerow(cols); w.writerows(rows)
    if undefined:
        print("UNDEFINED GENES (add to config/gene_targets.tsv):")
        for eid,g in sorted(undefined): print(f"   {eid}: {g}")
    print(f"wrote {OUT}  ({len(rows)} rows)")

if __name__=="__main__": main()
