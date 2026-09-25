#!/usr/bin/env python3
"""Mine the level-2 raw JSON for everything beyond the feature table. No requests.

Four payloads the first pass left on the floor:

  1. CAUTION / SEQUENCE CAUTION comments -- a DIFFERENT field from the
     "Sequence conflict" feature. Subtypes include Erroneous initiation,
     Frameshift, Erroneous termination. This is where start-codon disputes live.
  2. Pfam / InterPro / NCBIfam / CDD / PROSITE / PIRSF domain accessions -- the
     scanning strategy is anchored on DOMAINS, not gene names, so these ids are
     the actual anchors. Their ABSENCE on an exclusion protein is itself a
     usable feature: "small domainless ORF inside a T4SS operon".
  3. SUBCELLULAR LOCATION comments -- curated IM vs OM, more reliable than a
     SignalP prediction, and a primary discriminator for entry vs surface.
  4. AlphaFoldDB / PDB cross-references -- existing structures, so the roadmap's
     AF2/AF3 step does not re-predict what is already deposited.
"""
import csv, json, os
PROJ=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L2=os.path.join(PROJ,"data","level2_uniprot")
DOMAIN_DBS={"Pfam","InterPro","NCBIfam","CDD","PROSITE","PIRSF","SUPFAM","Gene3D"}
STRUCT_DBS={"AlphaFoldDB","PDB","SMR"}

def main():
    raw=json.load(open(os.path.join(L2,"uniprot_raw.json")))
    idx={}
    for r in csv.DictReader(open(os.path.join(L2,"uniprot_entries.tsv")),delimiter="\t"):
        if r["uniprot"]: idx[r["uniprot"]]=r
    rows=[]
    for acc,e in raw.items():
        meta=idx.get(acc,{})
        cautions=[]; subloc=[]; func=[]
        for c in e.get("comments",[]):
            t=c.get("commentType")
            if t in ("CAUTION","SEQUENCE CAUTION"):
                txt=(c.get("texts") or [{}])[0].get("value","")
                sub=c.get("sequenceCautionType","")
                cautions.append(f"[{t}{'/'+sub if sub else ''}] {txt}")
            elif t=="SUBCELLULAR LOCATION":
                for l in c.get("subcellularLocations",[]):
                    v=(l.get("location") or {}).get("value","")
                    top=(l.get("topology") or {}).get("value","")
                    subloc.append(f"{v}{' ('+top+')' if top else ''}")
            elif t=="FUNCTION":
                func.append((c.get("texts") or [{}])[0].get("value",""))
        doms=[]; structs=[]
        for x in e.get("uniProtKBCrossReferences",[]):
            db=x.get("database")
            if db in DOMAIN_DBS:
                props={p["key"]:p["value"] for p in x.get("properties",[])}
                nm=props.get("EntryName") or props.get("MatchStatus") or ""
                doms.append(f"{db}:{x['id']}{'('+nm+')' if nm and nm!='1' else ''}")
            elif db in STRUCT_DBS:
                structs.append(f"{db}:{x['id']}")
        rows.append({"uniprot":acc,
                     "element_id":meta.get("element_id",""),"gene":meta.get("gene",""),
                     "side":meta.get("side",""),"reviewed":meta.get("reviewed",""),
                     "length":e["sequence"]["length"],
                     "n_domains":len(doms),"domains":";".join(sorted(doms)),
                     "subcellular_location":" | ".join(subloc),
                     "cautions":" | ".join(cautions),
                     "structures":";".join(sorted(structs)),
                     "function":(" ".join(func))[:300]})
    cols=["element_id","gene","side","uniprot","reviewed","length","n_domains",
          "domains","subcellular_location","cautions","structures","function"]
    out=os.path.join(L2,"uniprot_mined.tsv")
    with open(out,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t", lineterminator="\n",extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows,key=lambda x:(x["element_id"],x["side"])): w.writerow(r)
    print(f"wrote {os.path.relpath(out,PROJ)}\n")
    print("=== CAUTIONS ===")
    n=0
    for r in rows:
        if r["cautions"]:
            n+=1; print(f"  {r['element_id']} {r['gene']} ({r['uniprot']}):\n     {r['cautions'][:260]}")
    if not n: print("  none")
    print("\n=== SUBCELLULAR LOCATION (curated) ===")
    for r in sorted(rows,key=lambda x:x["element_id"]):
        if r["subcellular_location"]:
            print(f"  {r['element_id']:<11} {r['gene']:<10} {r['subcellular_location']}")
    print("\n=== DOMAIN ANCHORS ===")
    for r in sorted(rows,key=lambda x:(x["side"],x["element_id"])):
        tag="EXCL" if r["side"]=="exclusion" else "TGT "
        print(f"  {tag} {r['element_id']:<11} {r['gene']:<10} n={r['n_domains']}  {r['domains'][:96]}")
    print("\n=== STRUCTURES ===")
    for r in sorted(rows,key=lambda x:x["element_id"]):
        print(f"  {r['element_id']:<11} {r['gene']:<10} {r['structures'] or 'none'}")

if __name__=="__main__": main()
