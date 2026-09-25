#!/usr/bin/env python3
"""Precision of the S4 discovery heuristic, measured against in-operon decoys.

The question this answers: inside the +/-6 kb window around the VirB6-family
target, how many ORFs pass the S4 step-3 length filter (40-260 aa) but are NOT
the exclusion gene? Those are the false positives the heuristic must reject, and
until they are counted, `supports_discovery_heuristic` carries no precision
information at all.

F is the worked case that prompted this: trbH and trbF sit in the SAME operon,
inside the window, and are not exclusion genes.
"""
import csv, glob, os, sys
from Bio import SeqIO
PROJ=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GB=os.path.join(PROJ,"data","seed","genbank")
HOSTS="/global/scratch/users/kh36969/exclusion_gene/hosts"
WINDOW=6000; AA_MIN,AA_MAX=40,260

def load(p):
    rows=[l.rstrip("\n").split("\t") for l in open(p) if not l.startswith("#") and l.strip()]
    return [dict(zip(rows[0],r)) for r in rows[1:]]

def find(acc):
    h=glob.glob(os.path.join(GB,f"*__{acc}.gb")) or glob.glob(os.path.join(HOSTS,f"*__{acc}.gb"))
    return h[0] if h else None

def main():
    out=[]
    print(f"{'element':<12} {'target':<10} {'in-window':>9} {'pass-len':>9} {'decoys':>7} "
          f"{'true d(bp)':>11} {'rank':>6}")
    print("-"*76)
    for t in load(os.path.join(PROJ,"config","extraction_targets.tsv")):
        pc,ec = t.get("partner_coords",""), t.get("exclusion_coords","")
        if not (pc and ec): continue
        acc=t["accession_authoritative"]; p=find(acc)
        if not p: continue
        rec=next(SeqIO.parse(p,"genbank"))
        pl,ph=(int(x) for x in pc.split(":")[1].split(".."))
        el,eh=(int(x) for x in ec.split(":")[1].split(".."))
        inwin=[];
        for f in rec.features:
            if f.type!="CDS": continue
            s,e=int(f.location.start)+1,int(f.location.end)
            if s==pl and e==ph: continue                       # the target itself
            d = 0 if (s<=ph and e>=pl) else (s-ph if s>ph else pl-e)
            if d>WINDOW: continue
            tr=(f.qualifiers.get("translation") or [""])[0]
            if not tr: continue
            inwin.append({"gene":(f.qualifiers.get("gene") or [""])[0],"aa":len(tr),
                          "d":d,"is_true":(s==el and e==eh)})
        passlen=[x for x in inwin if AA_MIN<=x["aa"]<=AA_MAX]
        decoys=[x for x in passlen if not x["is_true"]]
        true=[x for x in passlen if x["is_true"]]
        td=true[0]["d"] if true else None
        rank=(sorted(passlen,key=lambda x:x["d"]).index(true[0])+1) if true else None
        print(f"{t['element_id']:<12} {t['partner_gene']:<10} {len(inwin):>9} {len(passlen):>9} "
              f"{len(decoys):>7} {(td if td is not None else '-'):>11} "
              f"{(f'{rank}/{len(passlen)}' if rank else '-'):>6}")
        out.append({"element_id":t["element_id"],"n_in_window":len(inwin),
                    "n_pass_length":len(passlen),"n_decoys":len(decoys),
                    "true_distance_bp":td,"true_rank_by_distance":rank,
                    "decoy_genes":";".join(x["gene"] or "?" for x in sorted(decoys,key=lambda y:y["d"]))})
    f=os.path.join(PROJ,"data","exports","decoy_analysis.tsv")
    with open(f,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(out[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(out)
    tot_d=sum(o["n_decoys"] for o in out); tot_p=sum(o["n_pass_length"] for o in out)
    r1=sum(1 for o in out if o["true_rank_by_distance"]==1)
    print("-"*76)
    print(f"length filter alone: {tot_p-tot_d}/{tot_p} = {100*(tot_p-tot_d)/tot_p:.0f}% precision")
    print(f"nearest-ORF rule   : {r1}/{len(out)} elements rank the true gene FIRST by distance")
    print(f"\nwrote {os.path.relpath(f,PROJ)}")

if __name__=="__main__": main()
