#!/usr/bin/env python3
"""Stage 1 -- propose an accession for each seed element, with evidence.

This script COMMITS NOTHING. It emits a review table; a human pins the accession
into config/pinned_accessions.tsv. That is deliberate: protocol S3 requires the
accession to come from the primary paper's data-availability statement, and the
tier-B traps (pCF10 -> Citrobacter, pVCR94 -> deltaX deletion derivative,
pLS20 -> Xylella pPLS206) are exactly the cases an automatic pick gets wrong.

Evidence per candidate:
  source        plsdb | ncbi
  organism_ok   does /organism match the asserted expected_organism?
  is_derivative does the title look like a deletion/derivative construct?
"""
import csv, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
import ncbi

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLSDB_NUCCORE = "/global/scratch/users/kh36969/plsdb/nuccore.csv"
OUT = os.path.join(PROJ, "data", "seed", "accession_candidates.tsv")

DERIVATIVE = re.compile(r"delta|Delta|\bmini\b|::|derivative|deletion|knockout|"
                        r"transposon|insertion mutant", re.I)

def load_elements():
    rows = []
    with open(os.path.join(PROJ, "config", "seed_elements.tsv")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.rstrip("\n").split("\t"))
    hdr, data = rows[0], rows[1:]
    return [dict(zip(hdr, r)) for r in data]

def search_patterns(el):
    """Terms to search on. Falls back to element_id; overridden via search_alias
    for ids that are bad search strings (bare 'F' matches 147 PLSDB rows) or
    that have a second historical name (RP4/RK2)."""
    alias = el.get("search_alias", "NA")
    if alias and alias != "NA":
        return alias.split("|")
    return [el["element_id"]]

def plsdb_hits(terms):
    """Search the local PLSDB description field. Word-boundary anchored so
    'R100' does not match 'R1000' and 'pCF10' does not match 'pCF10-H'."""
    pat = re.compile("|".join(rf"(?<![A-Za-z0-9]){re.escape(t)}(?![A-Za-z0-9-])"
                              for t in terms), re.I)
    out = []
    with open(PLSDB_NUCCORE, newline="") as f:
        for row in csv.DictReader(f):
            if pat.search(row["NUCCORE_Description"]):
                out.append({"source": "plsdb",
                            "accver": row["NUCCORE_ACC"],
                            "title": row["NUCCORE_Description"],
                            "slen": row["NUCCORE_Length"],
                            "organism": ""})
    return out

def ncbi_hits(terms, expected_organism, retmax=12):
    org = (f' AND "{expected_organism}"[Organism]'
           if expected_organism and expected_organism != "NA" else "")
    out = []
    for t in terms:
        term = f'"{t}"[Title] AND complete[Title]{org}'
        try:
            out += ncbi.esummary(ncbi.esearch(term, retmax=retmax))
        except Exception as e:
            sys.stderr.write(f"  !! esearch failed [{term}]: {e}\n")
    return out

def main():
    elements = load_elements()
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["element_id", "role", "tier", "expected_organism", "source",
                    "accver", "length", "organism", "organism_ok",
                    "is_derivative", "title"])
        for el in elements:
            eid, exp_org = el["element_id"], el["expected_organism"]
            sys.stderr.write(f"[{eid}] tier {el['resolution_tier']} ...\n")
            terms = search_patterns(el)
            cands = plsdb_hits(terms) + ncbi_hits(terms, exp_org)
            if not cands:
                w.writerow([eid, el["role"], el["resolution_tier"], exp_org,
                            "NONE", "", "", "", "", "", "NO CANDIDATE FOUND"])
                continue
            seen = set()
            for c in cands:
                if c["accver"] in seen:
                    continue
                seen.add(c["accver"])
                blob = f"{c['organism']} {c['title']}"
                if exp_org in ("", "NA"):
                    org_ok = "unasserted"
                else:
                    org_ok = "YES" if re.search(re.escape(exp_org), blob, re.I) else "NO"
                deriv = "YES" if DERIVATIVE.search(c["title"]) else ""
                w.writerow([eid, el["role"], el["resolution_tier"], exp_org,
                            c["source"], c["accver"], c["slen"], c["organism"],
                            org_ok, deriv, c["title"][:110]])
    sys.stderr.write(f"\nwrote {OUT}\n")

if __name__ == "__main__":
    main()
