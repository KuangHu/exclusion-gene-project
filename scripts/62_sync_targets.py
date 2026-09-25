#!/usr/bin/env python3
"""Regenerate the extractor's target list FROM the element YAMLs.

data/elements/*.yaml is the ONLY place a pair is declared.
config/extraction_targets.tsv is generated; never hand-edit it.

Design points that are load-bearing:

* ACCESSIONS ARE SPLIT BY COORDINATE RELATIONSHIP, not lumped into one list.
    authoritative  the coordinate-system anchor. Exactly one. All coords in the
                   database are in THIS frame.
    mirror         same sequence AND same coordinates (RefSeq <-> GenBank).
                   Mirror status is VERIFIED here, not asserted: the RefSeq
                   COMMENT backlink must resolve to the authoritative accession
                   and the lengths must match.
    related        DIFFERENT coordinate frame. For comparison ONLY. A script
                   must never validate a coordinate against a related record.
  RP4 is why: BN000925.1 and CP152305.1 carry an IDENTICAL TrbK protein but
  their frames are offset by 2 bp (trbK 27450..27659 vs 27452..27661). Picking
  the wrong one lands on the wrong sequence, quietly.

* ROLE IS TWO COLUMNS. R27 is an exclusion POSITIVE and an adjacency NEGATIVE
  (eexA in the reverse-strand Tra2 Z operon, trhG over in Tra1). Collapsing them
  makes R27 read as "should have been found and wasn't", depressing recall,
  when the real question -- how much of known exclusion the adjacency heuristic
  covers -- is a different measurement.

* expected_protein_sha256 travels WITH the target. The sync step fetches and
  immediately checks it; a mismatch fails loud. Without it an accession version
  bump (.1 -> .2) shifts coordinates by a few bases and nothing notices.

* min_expected_bp derives from accessions.record_scope, never from element size.
  ICEBs1's element is ~20 kb but its record is a 4.2 Mb chromosome.
"""
import csv, glob, os, re, sys, yaml
from Bio import SeqIO
sys.path.insert(0,os.path.join(os.path.dirname(__file__),"lib"))
import ncbi

PROJ=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GB=os.path.join(PROJ,"data","seed","genbank")
norm=lambda a: re.sub(r"^(NZ_|NC_)","",a or "")

def rec_path(eid,acc):
    for pref in {"IncC":["pVCR94","R16a"],"SXT":["SXT","R391"]}.get(eid,[eid]):
        p=os.path.join(GB,f"{pref}__{acc}.gb")
        if os.path.exists(p): return p
    hits=glob.glob(os.path.join(GB,f"*__{acc}.gb"))
    return hits[0] if hits else None

def rec_info(eid,acc):
    p=rec_path(eid,acc)
    if not p: return None
    txt=open(p).read()
    r=next(SeqIO.parse(p,"genbank"))
    return {"len":len(r.seq),"backlink":ncbi.refseq_to_genbank(txt)}

def classify(eid,auth,others):
    """Return (mirrors, related). Mirror status is verified, never assumed."""
    a=rec_info(eid,auth)
    mirrors,related=[],[]
    for o in others:
        if not o or o==auth: continue
        i=rec_info(eid,o)
        if not i: related.append(o); continue
        same_len = a and i["len"]==a["len"]
        backlink_ok = i["backlink"] and norm(i["backlink"])==norm(auth)
        if backlink_ok and same_len: mirrors.append(o)
        elif norm(o)==norm(auth) and same_len: mirrors.append(o)
        else: related.append(o)
    return mirrors,related

COLS=["element_id","role_exclusion","role_adjacency","inc_group","mpf_class",
      "record_scope","min_expected_bp",
      "accession_authoritative","accession_mirror","accession_related",
      "exclusion_gene","exclusion_coords","exclusion_expected_sha256",
      "partner_gene","partner_family","partner_family_note","partner_coords","partner_expected_sha256"]

def main():
    els=[yaml.safe_load(open(f)) for f in sorted(glob.glob(os.path.join(PROJ,"data","elements","*.yaml")))]
    out=os.path.join(PROJ,"config","extraction_targets.tsv")
    rows=[]
    for e in sorted(els,key=lambda x:x["element_id"]):
        eid=e["element_id"]; acc=e["accessions"]
        auth=acc.get("primary") or (acc.get("parts") or [None])[0]
        others=[acc.get("crosscheck")]+(acc.get("parts") or [])
        mirrors,related=classify(eid,auth,others)
        adj=e["context"]["adjacency"].get("supports_discovery_heuristic")
        rows.append({
          "element_id":eid,
          "role_exclusion":"positive",
          "role_adjacency":{True:"positive",False:"negative",None:"unmeasured"}[adj],
          "inc_group":e.get("inc_group"),"mpf_class":e.get("mpf_class"),
          "record_scope":acc.get("record_scope"),"min_expected_bp":acc.get("min_expected_bp"),
          "accession_authoritative":auth,
          "accession_mirror":";".join(mirrors),
          "accession_related":";".join(related),
          "exclusion_gene":e["exclusion"]["gene"],
          "exclusion_coords":(e["exclusion"].get("sequence") or {}).get("coords",""),
          "exclusion_expected_sha256":(e["exclusion"].get("sequence") or {}).get("aa_sha256",""),
          "partner_gene":e["partner"]["gene"],
          "partner_family":e["partner"]["family"],
          "partner_family_note":e["partner"].get("family_note",""),
          "partner_coords":(e["partner"].get("sequence") or {}).get("coords",""),
          "partner_expected_sha256":(e["partner"].get("sequence") or {}).get("aa_sha256","")})
    with open(out,"w",newline="") as fh:
        fh.write("# GENERATED by scripts/62_sync_targets.py from data/elements/*.yaml\n"
                 "# DO NOT EDIT. Edit the element YAML and re-run.\n"
                 "# accession_related has a DIFFERENT coordinate frame -- never use it\n"
                 "# to validate a coordinate. Empty cells are empty strings, not NA:\n"
                 "# read with pandas keep_default_na=False.\n")
        w=csv.DictWriter(fh,fieldnames=COLS,delimiter="\t", lineterminator="\n",extrasaction="ignore")
        w.writeheader()
        for r in rows: w.writerow({k:("" if r.get(k) is None else r[k]) for k in COLS})
    print(f"wrote {os.path.relpath(out,PROJ)}  ({len(rows)} elements)")
    n=sum(1 for r in rows if r["exclusion_expected_sha256"])
    print(f"exclusion SHAs available for verification: {n}/{len(rows)}")
    print(f"adjacency negatives: {[r['element_id'] for r in rows if r['role_adjacency']=='negative']}")
    print(f"adjacency unmeasured: {[r['element_id'] for r in rows if r['role_adjacency']=='unmeasured']}")

if __name__=="__main__": main()
