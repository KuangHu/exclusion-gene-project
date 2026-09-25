#!/usr/bin/env python3
"""LEVEL 2 -- UniProt (Swiss-Prot preferred) for every frozen seed protein.

Why this level exists: the curated entries carry three things the GenBank
records do not.

  1. SEQUENCE CONFLICT features -- a curated map of where different papers or
     accessions reported DIFFERENT sequences. This is the historical dispute
     record. F TraS (P09129) carries "Sequence conflict 144..173 in Ref. 1 and 2";
     the 1987 submission CAA30011.1 is 149 aa and diverges from residue 144,
     24 residues short of the modern 173 aa. Everything before 144 is identical.
  2. SIGNAL / LIPID / TRANSMEM features -- curated lipobox and topology calls
     rather than my regex guesses.
  3. Reference lists with PMIDs, which feed back into level 1.

Cross-check: UniProt's sequence SHA is compared against the frozen level-1 SHA.
A mismatch means the two databases disagree about the protein, which must be
resolved before either is cited.

LOOKUP ROUTE -- exact, not name-based. The frozen protein's CRC64 checksum is
resolved through UniParc to a UniParc identifier, then UniProtKB is queried for
the ACTIVE entries on that identifier, preferring Swiss-Prot. Searching by my own
EMBL protein accession fails for 14 of 18 proteins, because UniProt cross-
references whichever EMBL record it curated from, not necessarily the one I
pinned. The checksum route is identity-by-sequence, so it cannot miss for that
reason -- and a hit is guaranteed to be the same protein.
"""
import csv, json, os, sys, time, urllib.parse, urllib.request
sys.path.insert(0,os.path.join(os.path.dirname(__file__),"lib"))
import ncbi

PROJ=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT=os.path.join(PROJ,"data","level2_uniprot")
API="https://rest.uniprot.org/uniprotkb"
_LAST=[0.0]

def _get(url,tries=4):
    for a in range(tries):
        dt=time.time()-_LAST[0]
        if dt<0.4: time.sleep(0.4-dt)
        _LAST[0]=time.time()
        try:
            with urllib.request.urlopen(url,timeout=90) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if a==tries-1: raise
            time.sleep(2**a)

UNIPARC="https://rest.uniprot.org/uniparc"

def search(term,size=5):
    q=urllib.parse.urlencode({"query":term,"size":size})
    return _get(f"{API}/search?{q}").get("results",[])

def upi_for(seq):
    """Exact sequence -> UniParc identifier, via CRC64 checksum."""
    from Bio.SeqUtils.CheckSum import crc64
    ck=crc64(seq).replace("CRC-","")
    q=urllib.parse.urlencode({"query":f"checksum:{ck}","size":1,"fields":"upi"})
    res=_get(f"{UNIPARC}/search?{q}").get("results",[])
    return (res[0]["uniParcId"] if res else None), ck

def kb_entries_for_upi(upi):
    """Active UniProtKB entries on a UniParc id, Swiss-Prot first."""
    hits=search(f"uniparc:{upi}",size=25)
    hits.sort(key=lambda h: 0 if "reviewed (Swiss-Prot)" in h.get("entryType","") else 1)
    return hits

def entry(acc):
    return _get(f"{API}/{acc}.json")

def load(p):
    rows=[l.rstrip("\n").split("\t") for l in open(p) if not l.startswith("#") and l.strip()]
    return [dict(zip(rows[0],r)) for r in rows[1:]]

FEATURES_WANTED={"Sequence conflict","Signal","Lipid binding","Lipidation",
                 "Transmembrane","Chain","Mutagenesis","Propeptide","Domain",
                 "Region","Modified residue","Topological domain"}

def _frozen_sequences():
    """protein_accession -> translation, read from the frozen GenBank records."""
    import glob as _g
    from Bio import SeqIO
    out={}
    for pat in (os.path.join(PROJ,"data","seed","genbank","*.gb"),
                "/global/scratch/users/kh36969/exclusion_gene/hosts/*.gb"):
        for f in _g.glob(pat):
            try: rec=next(SeqIO.parse(f,"genbank"))
            except Exception: continue
            for ft in rec.features:
                if ft.type!="CDS": continue
                pid=(ft.qualifiers.get("protein_id") or [""])[0]
                tr=(ft.qualifiers.get("translation") or [""])[0]
                if pid and tr: out.setdefault(pid,tr)
    return out

def main():
    os.makedirs(OUT,exist_ok=True)
    tg=load(os.path.join(PROJ,"config","extraction_targets.tsv"))
    # every frozen protein: (element, gene, side, protein_accession, frozen_sha)
    want=[]
    import glob, yaml
    for f in sorted(glob.glob(os.path.join(PROJ,"data","elements","*.yaml"))):
        e=yaml.safe_load(open(f))
        for side in ("exclusion","partner"):
            sq=e[side].get("sequence") or {}
            if sq.get("protein_accession"):
                want.append((e["element_id"],e[side]["gene"],side,
                             sq["protein_accession"],sq.get("aa_sha256"),sq.get("length_aa")))
    rows=[]; entries={}
    seqs=_frozen_sequences()
    for eid,gene,side,pacc,sha,alen in want:
        seq=seqs.get(pacc)
        upi,ck=(None,None)
        hits=[]
        if seq:
            upi,ck=upi_for(seq)
            if upi: hits=kb_entries_for_upi(upi)
        route="checksum" if hits else ""
        if not hits:
            hits=search(pacc.split(".")[0],size=5)
            hits.sort(key=lambda h: 0 if "reviewed (Swiss-Prot)" in h.get("entryType","") else 1)
            if hits: route="accession_fallback"
        if not hits:
            route="none"
            rows.append({"element_id":eid,"gene":gene,"side":side,
                         "embl_protein":pacc,"uniprot":"","reviewed":"NOT_FOUND",
                         "uniprot_len":"","frozen_len":alen,"sha_match":"",
                         "n_conflicts":"","conflicts":"","features":"","pmids":"",
                         "uniparc":upi or "","crc64":ck or "","lookup_route":route})
            print(f"  {eid:<11} {gene:<10} {pacc:<14} NOT FOUND (uniparc={upi})")
            continue
        h=hits[0]; acc=h["primaryAccession"]
        d=entry(acc); entries[acc]=d
        useq=d["sequence"]["value"]; usha=ncbi.aa_sha256(useq)
        reviewed="reviewed" if "reviewed (Swiss-Prot)" in d.get("entryType","") else "unreviewed"
        feats=[]; confl=[]
        for f in d.get("features",[]):
            t=f["type"]
            if t not in FEATURES_WANTED: continue
            s=f["location"]["start"].get("value"); e2=f["location"]["end"].get("value")
            rec=f"{t}:{s}-{e2}"
            if f.get("description"): rec+=f"({f['description'][:40]})"
            feats.append(rec)
            if t=="Sequence conflict": confl.append(rec)
        pmids=[]
        for r in d.get("references",[]):
            for x in r.get("citation",{}).get("citationCrossReferences",[]):
                if x["database"]=="PubMed": pmids.append(x["id"])
        m = "MATCH" if usha==sha else "DIFFER"
        rows.append({"element_id":eid,"gene":gene,"side":side,
                     "embl_protein":pacc,"uniprot":acc,"reviewed":reviewed,
                     "uniprot_len":d["sequence"]["length"],"frozen_len":alen,
                     "sha_match":m,"n_conflicts":len(confl),
                     "conflicts":";".join(confl),"features":";".join(feats),
                     "pmids":";".join(sorted(set(pmids))),
                     "uniparc":upi or "","crc64":ck or "",
                     "lookup_route":route})
        flag=" *** CONFLICT ***" if confl else ""
        print(f"  {eid:<11} {gene:<10} {pacc:<14} {acc:<8} {reviewed:<10} "
              f"{d['sequence']['length']:>4}aa {m:<7}{flag}")
    cols=["element_id","gene","side","embl_protein","crc64","uniparc","uniprot",
          "reviewed","uniprot_len","frozen_len","sha_match","n_conflicts",
          "conflicts","features","pmids","lookup_route"]
    with open(os.path.join(OUT,"uniprot_entries.tsv"),"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t", lineterminator="\n",extrasaction="ignore")
        w.writeheader()
        for r in rows: w.writerow({c:r.get(c,"") for c in cols})
    json.dump(entries,open(os.path.join(OUT,"uniprot_raw.json"),"w"),indent=1)
    print(f"\n{len(rows)} proteins queried")
    print(f"reviewed (Swiss-Prot) : {sum(1 for r in rows if r['reviewed']=='reviewed')}")
    print(f"unreviewed (TrEMBL)   : {sum(1 for r in rows if r['reviewed']=='unreviewed')}")
    print(f"not in UniProt        : {sum(1 for r in rows if r['reviewed']=='NOT_FOUND')}")
    print(f"SHA mismatch vs level1: {sum(1 for r in rows if r['sha_match']=='DIFFER')}")
    print(f"entries WITH sequence conflicts: {sum(1 for r in rows if r['n_conflicts'] not in ('',0))}")

if __name__=="__main__": main()
