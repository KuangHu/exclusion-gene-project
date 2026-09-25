#!/usr/bin/env python3
"""Flag UniProt Signal features that likely ignore a lipobox.

Generalises the RP4 TrbK finding. Q79AR9 gives Signal 1..21 / Chain 22..69, but
residue 22 is Gly and a lipoprotein's mature +1 MUST be the lipidated cysteine at
23. haase1996 is right; the UniProt call looks like SignalP applied without
lipoprotein awareness.

That failure mode applies to every lipoprotein in the set, and there are several:
pKM101 Eex (pohlman1994's title names a lipid attachment motif), F TraT, IncC Sfx,
R27 EexB. This checks all of them mechanically rather than one at a time.

Rule: if UniProt annotates a Signal ending at position p, and a canonical lipobox
[LVIAM][ASTVI][GAS]C places a cysteine at p+2 rather than p+1, the boundary is
off by one and gets uniprot_signal_suspect=true.
"""
import csv, json, os, re
PROJ=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L2=os.path.join(PROJ,"data","level2_uniprot")
LIPOBOX=re.compile(r"[LVIAM][ASTVI][GAS]C")

def main():
    raw=json.load(open(os.path.join(L2,"uniprot_raw.json")))
    meta={}
    for r in csv.DictReader(open(os.path.join(L2,"uniprot_entries.tsv")),delimiter="\t"):
        if r["uniprot"]: meta[r["uniprot"]]=r
    rows=[]
    for acc,e in raw.items():
        seq=e["sequence"]["value"]
        m=meta.get(acc,{})
        sig=[f for f in e.get("features",[]) if f["type"]=="Signal"]
        # every cysteine reachable by a lipobox in the N-terminal region
        lip=[mm.end() for mm in LIPOBOX.finditer(seq[:45])]   # 1-based Cys position
        if not sig and not lip: continue
        sig_end=sig[0]["location"]["end"].get("value") if sig else None
        verdict="ok"; detail=""
        if sig_end and lip:
            cys=lip[0]
            if cys==sig_end+1: verdict="ok"; detail=f"mature +1 = Cys{cys}, consistent"
            elif cys==sig_end+2:
                verdict="SUSPECT_OFF_BY_ONE"
                detail=(f"Signal ends {sig_end} so mature +1 = {seq[sig_end-1+1-1]}{sig_end+1}, "
                        f"but the lipobox Cys is at {cys}. A lipoprotein's mature +1 must be the Cys.")
            else:
                verdict="SUSPECT_MISMATCH"
                detail=f"Signal ends {sig_end}, lipobox Cys at {cys} (offset {cys-sig_end})"
        elif lip and not sig:
            verdict="LIPOBOX_UNANNOTATED"
            detail=f"lipobox Cys at {lip[0]} but UniProt annotates no Signal"
        elif sig and not lip:
            verdict="signal_no_lipobox"
            detail=f"Signal 1..{sig_end}, no canonical lipobox -- SPase I substrate, not a lipoprotein"
        rows.append({"element_id":m.get("element_id",""),"gene":m.get("gene",""),
                     "uniprot":acc,"reviewed":m.get("reviewed",""),
                     "length":e["sequence"]["length"],
                     "uniprot_signal_end":sig_end or "","lipobox_cys":lip[0] if lip else "",
                     "uniprot_signal_suspect":str(verdict.startswith("SUSPECT")).lower(),
                     "verdict":verdict,"detail":detail,"nterm":seq[:32]})
    out=os.path.join(L2,"signal_audit.tsv")
    cols=["element_id","gene","uniprot","reviewed","length","uniprot_signal_end",
          "lipobox_cys","uniprot_signal_suspect","verdict","detail","nterm"]
    with open(out,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t", lineterminator="\n",extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows,key=lambda x:(x["verdict"],x["element_id"])): w.writerow(r)
    print(f"{'element':<11} {'gene':<10} {'uniprot':<11} {'sigend':>6} {'cys':>4}  verdict")
    print("-"*82)
    for r in sorted(rows,key=lambda x:(0 if x['uniprot_signal_suspect']=='true' else 1,x['element_id'])):
        print(f"{r['element_id']:<11} {r['gene']:<10} {r['uniprot']:<11} "
              f"{str(r['uniprot_signal_end']):>6} {str(r['lipobox_cys']):>4}  {r['verdict']}")
        if r['detail'] and r['verdict']!='ok': print(f"      {r['detail'][:96]}")
    n=sum(1 for r in rows if r['uniprot_signal_suspect']=='true')
    print(f"\n{len(rows)} proteins with a Signal feature or an N-terminal lipobox")
    print(f"{n} flagged uniprot_signal_suspect")
    print(f"\nwrote {os.path.relpath(out,PROJ)}")

if __name__=="__main__": main()
