#!/usr/bin/env python3
"""Generate the machine-readable exports S7 requires, from the YAML source.

  data/exports/entry_exclusion_pairs.tsv   flat, one row per pair -- greppable,
                                           joinable, diffable
  data/exports/entry_exclusion_db.json     nested, full context

Never edit the exports. Edit data/elements/*.yaml and re-run.
"""
import glob, json, os, yaml, csv
PROJ=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT=os.path.join(PROJ,"data","exports"); os.makedirs(OUT,exist_ok=True)

els=[yaml.safe_load(open(f)) for f in sorted(glob.glob(os.path.join(PROJ,"data","elements","*.yaml")))]
els.sort(key=lambda e:e["element_id"])

def g(d,*ks,default=None):
    for k in ks:
        if not isinstance(d,dict): return default
        d=d.get(k)
        if d is None: return default
    return d

COLS=["element_id","element_type","inc_group","mpf_class","host",
      "accession_primary","accession_crosscheck","accession_status",
      "exclusion_gene","exclusion_type","evidence_tier","exclusion_protein",
      "exclusion_aa_len","exclusion_aa_sha256","exclusion_coords","exclusion_strand",
      "partner_gene","partner_family","partner_virb","partner_status",
      "partner_specificity_residues",
      "adjacency_bp","within_discovery_window","same_strand","same_operon",
      "genes_between","supports_discovery_heuristic","supports_synteny_claim"]

rows=[]
for e in els:
    ex,pa,ad = e["exclusion"], e["partner"], g(e,"context","adjacency",default={})
    rows.append({
     "element_id":e["element_id"],"element_type":e["element_type"],
     "inc_group":e.get("inc_group"),"mpf_class":e.get("mpf_class"),"host":e.get("host"),
     "accession_primary":g(e,"accessions","primary"),
     "accession_crosscheck":g(e,"accessions","crosscheck"),
     "accession_status":g(e,"accessions","status"),
     "exclusion_gene":ex.get("gene"),"exclusion_type":ex.get("type"),
     "evidence_tier":ex.get("evidence_tier"),
     "exclusion_protein":g(ex,"sequence","protein_accession"),
     "exclusion_aa_len":g(ex,"sequence","length_aa"),
     "exclusion_aa_sha256":g(ex,"sequence","aa_sha256"),
     "exclusion_coords":g(ex,"sequence","coords"),
     "exclusion_strand":g(ex,"sequence","strand"),
     "partner_gene":pa.get("gene"),"partner_family":pa.get("family"),
     "partner_virb":pa.get("virb_equivalent"),"partner_status":pa.get("status"),
     "partner_specificity_residues":pa.get("specificity_residues"),
     "adjacency_bp":g(ad,"proximity","measured_bp"),
     "within_discovery_window":g(ad,"proximity","within_discovery_window"),
     "same_strand":g(ad,"synteny","same_strand"),
     "same_operon":g(ad,"synteny","same_operon"),
     "genes_between":g(ad,"synteny","genes_between"),
     "supports_discovery_heuristic":ad.get("supports_discovery_heuristic"),
     "supports_synteny_claim":ad.get("supports_synteny_claim")})

tsv=os.path.join(OUT,"entry_exclusion_pairs.tsv")
with open(tsv,"w",newline="") as fh:
    fh.write("# supports_discovery_heuristic is a RECALL statement: the exclusion gene\n"
             "# lies inside the +/-6 kb window. It is NOT precision -- 53 non-exclusion\n"
             "# ORFs also satisfy it (11%% precision). See docs/heuristic_precision.md.\n"
             "# Empty cells are empty strings, not NA: read with keep_default_na=False.\n")
    w=csv.DictWriter(fh,fieldnames=COLS,delimiter="\t", lineterminator="\n",extrasaction="ignore")
    w.writeheader()
    for r in rows: w.writerow({k:("" if r.get(k) is None else r.get(k)) for k in COLS})

js=os.path.join(OUT,"entry_exclusion_db.json")
json.dump({"schema_version":"1.0.0","scope":"entry_exclusion_only",
           "n_elements":len(els),"elements":els}, open(js,"w"), indent=2, sort_keys=False)

print(f"{len(rows)} pairs -> {os.path.relpath(tsv,PROJ)}")
print(f"{len(els)} elements -> {os.path.relpath(js,PROJ)}")
d=[r["supports_discovery_heuristic"] for r in rows]
print(f"\ndiscovery heuristic  holds {d.count(True)}  fails {d.count(False)}  unmeasured {d.count(None)}")
print(f"synteny claim        holds {sum(1 for r in rows if r['supports_synteny_claim'])}/{len(rows)}")
print(f"target is VirB6-family     {sum(1 for r in rows if r['partner_family']=='VirB6')}/{len(rows)}")
print(f"target experimentally verified  {sum(1 for r in rows if r['partner_status']=='verified')}/{len(rows)}")
