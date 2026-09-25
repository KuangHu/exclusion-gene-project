#!/usr/bin/env python3
"""Run the hard-fail guards over every config and element file. Exit 1 on failure."""
import glob, os, sys, yaml, csv
sys.path.insert(0,os.path.join(os.path.dirname(__file__),"lib"))
from assertions import (check_enum, check_inc_group, AssertionFailure, ENUMS,
                        assert_no_cr, is_blank, assert_curated_unchanged,
                        assert_no_generated_collision, assert_accessions_resolve)

PROJ=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fails=[]

def guard(fn,*a,**k):
    try: fn(*a,**k)
    except AssertionFailure as e: fails.append(str(e))

# element YAMLs
for f in sorted(glob.glob(os.path.join(PROJ,"data","elements","*.yaml"))):
    e=yaml.safe_load(open(f)); w=os.path.basename(f)
    guard(check_inc_group, e.get("inc_group"), w)
    guard(check_enum,"record_scope",(e.get("accessions") or {}).get("record_scope"),w)
    guard(check_enum,"status",(e.get("accessions") or {}).get("status"),w)
    guard(check_enum,"type",e["exclusion"].get("type"),w)
    guard(check_enum,"evidence_tier",e["exclusion"].get("evidence_tier"),w)
    guard(check_enum,"hmm_hit",e["exclusion"].get("hmm_hit"),w)
    guard(check_enum,"partner_status",e["partner"].get("status"),w)

# gene_targets: restrict_inc_group must be NA or real Inc groups, never prose
gt=os.path.join(PROJ,"config","gene_targets.tsv")
rows=[l.rstrip("\n").split("\t") for l in open(gt) if not l.startswith("#") and l.strip()]
hdr=rows[0]
for r in rows[1:]:
    d=dict(zip(hdr,r)); w=f"gene_targets.tsv[{d['gene_canonical']}]"
    guard(check_enum,"side",d.get("side"),w)
    guard(check_inc_group, d.get("restrict_inc_group","NA"), w)
    for col in ("aa_min","aa_max"):
        v=d.get(col,"NA")
        if v not in ("NA","") and not v.isdigit():
            fails.append(f"{w}: {col}={v!r} is not numeric -- column shift?")

# extraction_targets
et=os.path.join(PROJ,"config","extraction_targets.tsv")
rows=[l.rstrip("\n").split("\t") for l in open(et) if not l.startswith("#") and l.strip()]
hdr=rows[0]
for r in rows[1:]:
    d=dict(zip(hdr,r)); w=f"extraction_targets.tsv[{d['element_id']}]"
    guard(check_enum,"role_exclusion",d.get("role_exclusion"),w)
    guard(check_enum,"role_adjacency",d.get("role_adjacency"),w)
    guard(check_enum,"record_scope",d.get("record_scope"),w)
    guard(check_inc_group,d.get("inc_group"),w)
    if not d.get("accession_authoritative"):
        fails.append(f"{w}: no accession_authoritative -- coordinates have no anchor")

# every gene referenced must be DEFINED (the L1 silent skip)
defined={r[0] for r in [l.rstrip("\n").split("\t") for l in open(gt)
                        if not l.startswith("#") and l.strip()][1:]}
for r in rows[1:]:
    d=dict(zip(hdr,r))
    for g in (d.get("exclusion_gene"),d.get("partner_gene")):
        if g and g not in defined:
            fails.append(f"extraction_targets.tsv[{d['element_id']}]: gene {g!r} is "
                         f"NOT defined in gene_targets.tsv -- it would be silently skipped")

# A7: no generated table may contain carriage returns
import glob as _g
for _p in (_g.glob(os.path.join(PROJ, "data", "**", "*.tsv"), recursive=True)
           + _g.glob(os.path.join(PROJ, "config", "*.tsv"))):
    try:
        assert_no_cr(_p)
    except AssertionFailure as e:
        fails.append(str(e))

# A9: curated files must still hash to their sealed values, and no script may
# open one for writing. This is the anchor_set.tsv clobber guard.
fails.extend(assert_curated_unchanged(PROJ))
fails.extend(assert_no_generated_collision(PROJ))
# A10: a cited accession must resolve to a model we actually hold.
fails.extend(assert_accessions_resolve(PROJ))

if fails:
    print(f"VALIDATION FAILED — {len(fails)} problem(s):\n")
    for f in fails: print("  * "+f+"\n")
    sys.exit(1)
print(f"validation PASSED — {len(glob.glob(os.path.join(PROJ,'data','elements','*.yaml')))} elements, "
      f"{len(rows)-1} target rows, all enums and references check out")
