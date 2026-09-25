#!/usr/bin/env python3
"""EMPIRICAL NULL for the TraG 610-673 enrichment claim.

The binomial p=4.2e-19 assumes substitutions are iid along the protein. They are
not: real proteins have conserved cores and variable loops, so SOME 64-aa window
will look enriched in almost any ortholog pair. The binomial therefore massively
overstates significance.

Correct test: build the distribution of MAXIMUM 64-aa window enrichment across
every F/R100 tra ortholog pair, then ask where TraG's PRE-SPECIFIED window falls
in it. Conservative on purpose -- TraG's window is fixed a priori while every
comparator gets its own best window.
"""
import os, sys, statistics
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

PROJ=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
W=64                      # the reported window width, 610-673
F_GB =os.path.join(PROJ,"data","seed","genbank","F__AP001918.1.gb")
R_GB =os.path.join(PROJ,"data","seed","genbank","R100__AP000342.1.gb")

def cds_by_gene(path):
    out={}
    rec=next(SeqIO.parse(path,"genbank"))
    for f in rec.features:
        if f.type!="CDS": continue
        g=((f.qualifiers.get("gene") or [""])[0]).strip()
        tr=(f.qualifiers.get("translation") or [""])[0]
        if g and tr: out.setdefault(g,tr)
    return out

def sub_profile(a,b,al):
    """Aligned substitution mask in seq-a coordinates. Gapped columns excluded."""
    aln=al.align(a,b)[0]; x,y=aln[0],aln[1]
    pos=0; mask=[]
    for p,q in zip(x,y):
        if p=="-": continue
        pos+=1
        if q=="-": mask.append(None)          # gap: not scoreable
        else: mask.append(p!=q)
    return mask

def max_window_enrichment(mask,w=W):
    """Best w-residue window: (enrichment, start, subs_in_window)."""
    sc=[m for m in mask if m is not None]
    n=len(sc); tot=sum(sc)
    if n<w or tot==0: return None
    rate=tot/n
    best=None
    run=sum(1 for m in mask[:w] if m)
    valid=sum(1 for m in mask[:w] if m is not None)
    for i in range(0,len(mask)-w+1):
        if i>0:
            if mask[i-1] is True: run-=1
            if mask[i-1] is not None: valid-=1
            if mask[i+w-1] is True: run+=1
            if mask[i+w-1] is not None: valid+=1
        if valid<w*0.8: continue
        exp=rate*valid
        if exp<=0: continue
        e=run/exp
        if best is None or e>best[0]: best=(e,i+1,run)
    return best

def window_enrichment(mask,lo,hi):
    sc=[m for m in mask if m is not None]
    n=len(sc); tot=sum(sc)
    seg=[m for m in mask[lo-1:hi] if m is not None]
    subs=sum(1 for m in mask[lo-1:hi] if m)
    exp=(tot/n)*len(seg)
    return (subs/exp if exp>0 else None), subs, exp, len(seg)

def main():
    al=PairwiseAligner(mode="global",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-11,extend_gap_score=-1)
    Fg,Rg=cds_by_gene(F_GB),cds_by_gene(R_GB)
    shared=sorted(set(Fg)&set(Rg))
    print(f"F/R100 orthologous pairs by gene name: {len(shared)}")
    print(f"of those with >= {W} aa and >0 substitutions:\n")
    results=[]
    for g in shared:
        a,b=Fg[g],Rg[g]
        if min(len(a),len(b))<W: continue
        mask=sub_profile(a,b,al)
        mw=max_window_enrichment(mask)
        if not mw: continue
        sc=[m for m in mask if m is not None]
        ident=100*(len(sc)-sum(sc))/len(sc)
        results.append({"gene":g,"len":len(a),"identity":ident,
                        "max_enrich":mw[0],"at":mw[1],"subs":mw[2],
                        "total_subs":sum(sc)})
    results.sort(key=lambda r:-r["max_enrich"])
    print(f"  {'gene':<8} {'len':>5} {'ident%':>7} {'max 64aa':>9} {'at':>5} {'subs':>5}")
    for r in results:
        star=" <== TraG (max-window, NOT the reported one)" if r["gene"]=="traG" else ""
        print(f"  {r['gene']:<8} {r['len']:>5} {r['identity']:>7.1f} {r['max_enrich']:>9.2f} "
              f"{r['at']:>5} {r['subs']:>5}{star}")

    # the pre-specified window
    mask=sub_profile(Fg["traG"],Rg["traG"],al)
    e,subs,exp,seglen=window_enrichment(mask,610,673)
    others=[r["max_enrich"] for r in results if r["gene"]!="traG"]
    print("\n" + "="*72)
    print("PRE-SPECIFIED WINDOW, F TraG 610-673")
    print(f"  substitutions {subs}  expected {exp:.2f}  enrichment {e:.2f}x")
    print("\nEMPIRICAL NULL -- max 64-aa enrichment over the other tra ortholog pairs")
    print(f"  n = {len(others)}")
    print(f"  median {statistics.median(others):.2f}x   mean {statistics.mean(others):.2f}x")
    print(f"  min {min(others):.2f}x   max {max(others):.2f}x")
    beat=sum(1 for o in others if o>=e)
    print(f"\n  comparator MAX windows >= TraG's pre-specified {e:.2f}x : {beat}/{len(others)}")
    print(f"  empirical p (conservative)                       : {(beat+1)/(len(others)+1):.4f}")
    print("\n  Note: every comparator is credited with its OWN BEST window while")
    print("  TraG's window is fixed a priori, so this is an upper bound on p.")

if __name__=="__main__": main()
