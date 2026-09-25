#!/usr/bin/env python3
"""Stage 2 -- fetch and FREEZE the seed records.

Records land in the repo (data/seed/genbank/), not scratch: S7 requires a stable
accession set per release. The one exception is CHROMOSOMAL hosts (e.g. the 4.2 Mb
B. subtilis 168 genome carrying ICEBs1) -- those are bulk and go to scratch; only
the extracted element is frozen in the repo.

Idempotent: an existing file whose md5 matches the manifest is not refetched.
"""
import csv, hashlib, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
import ncbi

PROJ    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GB_DIR  = os.path.join(PROJ, "data", "seed", "genbank")
SCRATCH = "/global/scratch/users/kh36969/exclusion_gene/hosts"
MANIFEST= os.path.join(PROJ, "data", "seed", "manifest.tsv")

def load(path):
    with open(path) as f:
        rows=[l.rstrip("\n").split("\t") for l in f if not l.startswith("#") and l.strip()]
    return [dict(zip(rows[0], r)) for r in rows[1:]]

def targets():
    """(element_id, accession, role_of_record, destination_dir)"""
    out=[]
    for el in load(os.path.join(PROJ,"config","pinned_accessions.tsv")):
        eid, st = el["element_id"], el["status"]
        if st=="UNRESOLVED": continue
        dest = SCRATCH if st=="CHROMOSOMAL" else GB_DIR
        for field, kind in (("genbank_accession","primary"),
                            ("refseq_accession","crosscheck")):
            acc=el[field]
            if acc and acc!="NA":
                out.append((eid,acc,kind,dest))
        aux=el["aux_accessions"]
        if aux and aux!="NA":
            for acc in aux.split(";"):
                if acc.strip(): out.append((eid,acc.strip(),"aux",dest))
    return out

def main():
    os.makedirs(GB_DIR,exist_ok=True); os.makedirs(SCRATCH,exist_ok=True)
    prev={}
    if os.path.exists(MANIFEST):
        for r in load(MANIFEST): prev[(r["element_id"],r["accession"])]=r
    rows=[]
    for eid,acc,kind,dest in targets():
        path=os.path.join(dest,f"{eid}__{acc}.gb")
        old=prev.get((eid,acc))
        if os.path.exists(path) and old:
            md5=hashlib.md5(open(path,'rb').read()).hexdigest()
            if md5==old["md5"]:
                sys.stderr.write(f"  skip  {eid:<12} {acc}\n")
                rows.append(old); continue
        sys.stderr.write(f"  FETCH {eid:<12} {acc:<16} ({kind})\n")
        try:
            txt=ncbi.efetch_gb(acc)
        except Exception as e:
            sys.stderr.write(f"    !! {e}\n")
            rows.append({"element_id":eid,"accession":acc,"record_role":kind,
                         "path":"","bp":"","md5":"","fetched":"","status":f"FETCH_FAILED: {e}"})
            continue
        if "LOCUS" not in txt[:200]:
            sys.stderr.write("    !! not a GenBank record\n"); continue
        open(path,"w").write(txt)
        bp=txt.split()[2] if len(txt.split())>2 else ""
        rows.append({"element_id":eid,"accession":acc,"record_role":kind,
                     "path":os.path.relpath(path,PROJ) if dest==GB_DIR else path,
                     "bp":bp,"md5":hashlib.md5(txt.encode()).hexdigest(),
                     "fetched":time.strftime("%Y-%m-%d"),"status":"OK"})
    cols=["element_id","accession","record_role","path","bp","md5","fetched","status"]
    with open(MANIFEST,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t", lineterminator="\n"); w.writeheader()
        for r in rows: w.writerow({c:r.get(c,"") for c in cols})
    ok=sum(1 for r in rows if r.get("status")=="OK")
    sys.stderr.write(f"\n{ok}/{len(rows)} records OK -> {MANIFEST}\n")

if __name__=="__main__": main()
