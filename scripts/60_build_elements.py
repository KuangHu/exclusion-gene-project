#!/usr/bin/env python3
"""Emit one YAML per element: curated operon context + verified sequence facts.

Curated fields come from the domain expert. Sequence fields (protein_accession,
length_aa, aa_sha256, coords, strand) are looked up from
data/seed/extraction_candidates.tsv so nothing is transcribed by hand, and
anything not yet recovered is written as null with a `needs` note rather than
guessed.

The emitted YAML is the SOURCE OF TRUTH from here on; this script exists to
create it once. Curators edit the YAML, not this file.
"""
import csv, os, re, sys, yaml
PROJ=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT=os.path.join(PROJ,"data","elements")

seq={}
with open(os.path.join(PROJ,"data","seed","extraction_candidates.tsv")) as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["route"] in ("MISS",""): continue
        k=(r["element_id"],r["gene"])
        pref={"primary":0,"aux":1,"crosscheck":2}
        if k not in seq or pref.get(r["record_role"],9)<pref.get(seq[k]["record_role"],9):
            seq[k]=r

def S(eid,gene,alias=None):
    r=seq.get((eid,gene)) or (seq.get((alias,gene)) if alias else None)
    if not r: return None
    return {"protein_accession":r["protein_id"] or None,
            "length_aa":int(r["aa_len"]), "aa_sha256":r["aa_sha256"],
            "coords":r["coords"], "strand":r["strand"], "route":r["route"]}

# ============================ CURATED CONTEXT ============================
# operon gene lists are ORDERED as they appear on the element.
# adjacency.supports_rule records whether the Garcillan-Barcia & de la Cruz
# observation (exclusion gene abuts the VirB6-family target) HOLDS for this
# element. It is a measured property with real FALSE values -- the database must
# be able to falsify the S4 heuristic, not only confirm it.
E=[
{"element_id":"F","element_type":"plasmid","inc_group":"IncFI","mpf_class":"MPF_F",
 "aa_len_status":"DISPUTED",
 "aa_len_dispute":"TraS length is NOT settled. UniProt P09129 carries 'Sequence conflict 144..173 in Ref. 1 and 2' -- the last 30 residues, 17% of the protein. The 1987 submission CAA30011.1 is 149 aa and diverges from residue 144; AP001918.1 gives 173 aa. audette2007's abstract reports no TraS length, so the C-terminus is unverified by literature, by UniProt, and between the two GenBank records. The frozen aa_sha256 therefore locks a C-terminus of unknown provenance. RP4 TrbK loses all activity to an 8-aa C-terminal truncation, so this conflict sits in the worst possible place.",
 "host":"Escherichia coli K-12","size_bp":99159,
 "accessions":{"primary":"AP001918.1","crosscheck":None,"parts":[],"status":"PINNED"},
 "exclusion":{"gene":"traS","type":"entry","evidence_tier":"E1",
   "specificity_determinant":"not mapped"},
 "partner":{"gene":"traG","family":"VirB6","virb_equivalent":"VirB6","status":"verified",
   "role":"N-term assembles pilus; C-term periplasmic domain TraG* handles Mps and docks TraN",
   "specificity_residues":"TraG 610-673",
   "specificity_window_provenance":"PRE-SPECIFIED, NOT DERIVED. The window came from the reference protocol S5.2 and was written into this record BEFORE any divergence scan was run, so the enrichment test below is NOT circular. HOWEVER the protocol's own source for 610-673 is unverified: audette2007's abstract gives no residue range. Confirm from Audette 2007 body text or the TraG* domain boundaries in bragagnolo2022 before publishing.",
   "specificity_enrichment":"6.88x. 31 substitutions vs 4.50 expected. Window identity 51.6% vs 93.0% overall.",
   "specificity_enrichment_pvalue":"EMPIRICAL p = 0.022 (scripts/46_window_null.py). SUPERSEDES the binomial p=4.2e-19, which was wrong by ~17 orders of magnitude: the binomial assumes substitutions are iid along the protein, but conserved cores and variable loops mean SOME 64-aa window looks enriched in almost any ortholog pair.",
   "specificity_enrichment_null":"Null = max 64-aa window enrichment across the other 45 F/R100 tra ortholog pairs: median 1.60x, mean 2.02x, range 0.73-6.75x. 0/45 comparators reach TraG's 6.88x, but traB reaches 6.75x -- essentially tied. TraG's window is the top by a hair, NOT an outlier.",
   "specificity_enrichment_caveat":"TraG's pre-specified 610-673 IS also its empirical maximum window (max-enrichment scan returns start=610). Reassuring for the window being real, but it removes the a-priori advantage IF the protocol's number was itself scan-derived -- which cannot be ruled out, since the same protocol sentence's other number (the 17%) was already shown to be imported from TraS.",
   "specificity_contrast":"TraS itself shows max window enrichment of only 1.06x despite 27.9% F/R100 identity -- diffusely divergent, not locally. So TraG has a discrete variable segment and TraS does not, which is the substantive result here rather than the p-value.",
   "specificity_measured":"MEASURED from pinned accessions (F BAA97969.1 vs R100 BAA78880.1): whole protein 92.8% identical; window 610-673 only 51.6% (33/64). The 601-640 and 641-680 windows are the two lowest across all 938 aa; every other 40-aa window is >=85%. The window claim holds.",
   "specificity_correction":"NOT ~17%. Per humbert2019, that figure is TraS(F) vs TraS(R100) identity -- measured here as 17.7% (36/203). The reference protocol S5.2 welded a TraS number onto a TraG window"},
 "context":{"operon":{"name":"tra operon (P_Y), single polycistron ~33 kb",
   "genes":["traY","traA","traL","traE","traK","traB","traP","trbD","trbG","traV","traR","traC","trbI","traW","traU","trbC","traN","trbE","traF","trbA","artA","traQ","trbB","trbJ","trbF","traH","traG","traS","traT","traD","trbH","traI","traX"],
   "orientation":"+","promoter":"P_Y; TraJ + ArcA co-activated, H-NS silenced",
   "outside_operon":["traM (own promoter)","traJ (own promoter, FinP/FinO antisense)","oriT","finO (inactivated in F itself)"],
   "note":"traS/traT/traD each have internal promoters. trbH lies between traD and traI, trbF between trbJ and traH -- both routinely omitted from simplified gene maps but present in the scan window"},
  "system":{"t4cp":"TraD","relaxase":"TraI","adhesin":"TraN"}},
 "notes":["Longest known E. coli transcript by some reports (~40 kb)"]},

{"element_id":"R64","element_type":"plasmid","inc_group":"IncI1","mpf_class":"MPF_I",
 "host":"Salmonella enterica ser. Typhimurium","size_bp":120826,
 "accessions":{"primary":"AP005147.1","crosscheck":"NC_005014.1","parts":[],"status":"PROVISIONAL"},
 "exclusion":{"gene":"excA","type":"entry","evidence_tier":"E1","second_gene":"excB",
   "specificity_determinant":None,
   "specificity_determinant_note":"sakuma2013 maps the discriminating segment on the DONOR TraY, not on ExcA. ExcA is the recogniser; the recognised element differs. Do not mirror the TraY segment onto this field",
   "note":"excA and excB are overlapping ORFs translated by in-frame reinitiation. furuya1994 is titled \"surface exclusion\" but the gene is entry exclusion -- another instance of title/product phenotype words not being reliable classifiers"},
 "partner":{"gene":"traY","family":"MPF_I inner-membrane Mpf protein","virb_equivalent":None,
   "status":"verified","specificity_residues":"TraY internal variable segment (R64); contrast R621a C-terminal",
   "specificity_source":"sakuma2013 -- reciprocal TraY swap",
   "note":"NOT VirB6. The R64/R621a Tra/Trb system is homologous to Legionella Dot/Icm and Ti VirB/D4, not to IncF. No TraF or TraG homolog exists in Dot/Icm or Ti -- there is no 'TraG' concept in this system at all"},
 "context":{"operon":{"name":"excA-traY cluster (one of four transfer-region clusters)",
   "genes":["excA","excB","traY"],"orientation":"+","promoter":"not characterised in the literature"},
  "transfer_region":{"size_kb":54,"n_genes":">=49","n_proteins":58,
   "clusters":{"regulation":["traA","traB","traC","traD"],
     "relaxosome":["oriT","nikA","nikB"],
     "type_IVb_thin_pilus":["pilI","pilJ","pilK","pilL","pilM","pilN","pilO","pilP","pilQ","pilR","pilS","pilT","pilU","pilV"],
     "conjugation_apparatus":["traE..traY (22 genes, Komano 2000)","trbA","trbB","trbC","nuc"]},
   "note":"58 proteins from >=49 genes via shufflon rearrangement plus overlapping genes. Thin pilus is required for liquid mating"}},
 "notes":["trbABC is adjacent to but antiparallel with the oriT operon","trbA and trbC essential for transfer; trbB deletion retains residual activity"]},

{"element_id":"R621a","element_type":"plasmid","inc_group":"IncI-gamma","mpf_class":"MPF_I",
 "host":"Salmonella enterica ser. Typhimurium","size_bp":93185,
 "accessions":{"primary":"AP011954.1","crosscheck":"NC_015965.1","parts":[],"status":"PROVISIONAL"},
 "exclusion":{"gene":"excA","type":"entry","evidence_tier":"E1",
   "specificity_determinant":None,
   "specificity_determinant_note":"sakuma2013 maps the discriminating segment on the DONOR TraY, not on ExcA"},
 "partner":{"gene":"traY","family":"MPF_I inner-membrane Mpf protein","virb_equivalent":None,
   "status":"verified","specificity_residues":"TraY C-terminal segment (R621a); contrast R64 internal",
   "specificity_source":"sakuma2013 -- reciprocal TraY swap"},
 "context":{"operon":{"name":"excA-traY cluster, as R64","genes":["excA","traY"],"orientation":"+"}},
 "notes":["CONTROL VALUE: every Tra/Trb protein of R64 and R621a is >95% identical EXCEPT ExcA and TraY. The cleanest available instance of exclusion specificity co-varying with its target and nothing else"]},

{"element_id":"RP4","element_type":"plasmid","inc_group":"IncP-alpha","mpf_class":"MPF_T",
 "host":"Escherichia coli","size_bp":60096,
 "accessions":{"primary":"BN000925.1","crosscheck":"CP152305.1","parts":[],"status":"PINNED"},
 "exclusion":{"gene":"trbK","type":"entry","evidence_tier":"E2",
   "specificity_determinant":"not mapped",
   "note":"47-aa mature lipoprotein after signal cleavage (69-aa precursor)"},
 "partner":{"gene":"trbL","family":"VirB6","virb_equivalent":"VirB6","status":"predicted",
   "specificity_residues":None,
   "note":"E4 hypothesis. Adjacency is confirmed but the interaction is untested"},
 "context":{"operon":{"name":"Tra2 trb operon (Mpf)",
   "genes":["trbA","trbB","trbC","trbD","trbE","trbF","trbG","trbH","trbI","trbJ","trbK","trbL","trbM","trbN","trbO","trbP"],
   "orientation":"+","promoter":"trbBp, main promoter of the mpf genes",
   "coords_kb":"Tra2 core 18.825-29.276 kb (18.03-29.26 kb whole region)"},
  "virb_mapping":{"trbB":"VirB11","trbC":"VirB2 (pilin)","trbD":"VirB3","trbE":"VirB4","trbF":"VirB8","trbG":"VirB9","trbJ":"VirB5","trbL":"VirB6","trbN":"VirB1"},
  "tra1":{"note":"DNA processing + coupling, 13 genes traA-traM in three operons",
   "operons":{"leading":["traK","traL","traM"],"relaxase":["traH","traI","traJ"],"primase":["traA","traB","traC","traD","traE","traF","traG"]},
   "system":{"t4cp":"TraG","relaxase":"TraI","primase":"TraC"}},
  "regulation":"TrbA (first gene of the trb operon) trans-represses traJp, traKp, traGp and itself; KorA/KorB global control"},
 "notes":["SIGNATURE: of the 12 core ORFs trbB-trbM, the ten trbB-trbL are all required for transfer EXCEPT trbK. A non-essential gene embedded in an otherwise essential mpf operon is the exclusion-gene signature"]},

{"element_id":"pKM101","element_type":"plasmid","inc_group":"IncN","mpf_class":"MPF_T",
 "host":"Escherichia coli","size_bp":None,
 "accessions":{"primary":None,"crosscheck":None,"parts":["U09868.1","U72482.2"],"status":"ASSEMBLE",
   "note":"no complete-plasmid record exists; U09868.1 carries the whole virB-colinear cluster plus eex"},
 "exclusion":{"gene":"eex","type":"entry","evidence_tier":"E2",
   "specificity_determinant":"not mapped"},
 "partner":{"gene":"traD","family":"VirB6","virb_equivalent":"VirB6","status":"predicted",
   "specificity_residues":None,
   "note":"CORRECTION: pKM101 TraG is VirB11 (331 aa ATPase, 3' end of the cluster), NOT the VirB6 homolog. TraD (346 aa) is VirB6 and is the gene immediately downstream of eex"},
 "context":{"operon":{"name":"tra mating-bridge cluster (virB1-virB11 fully colinear)",
   "genes":["traL","traM","traA","traB","traC","eex","traD","traN","traE","traO","traF","traG"],
   "orientation":"+","promoter":"two promoters, upstream of traL and upstream of traN"},
  "virb_mapping":{"traL":"VirB1","traM":"VirB2 (pilin)","traA":"VirB3","traB":"VirB4","traC":"VirB5","traD":"VirB6","traN":"VirB7","traE":"VirB8","traO":"VirB9","traF":"VirB10","traG":"VirB11"},
  "system":{"relaxase":"TraI (43% identical to R388 TrwC)","t4cp":"TraJ","accessory":"TraK (TrwA-like)",
   "omcc":"TraN/TraO/TraF core, 1.05 MDa, 14 copies (PDB 2YPW / EMD-2232)","second_adhesin":"pep"}},
 "notes":["VERIFIED 2026-09-01 from the U09868.1 feature table: eex 6533..6760 (+), traC ends 6525 (7 bp gap), traD starts 6776 (15 bp gap). The GenBank DEFINITION line places eex after traG, but that line contains a duplicated '(traB), (traC)' which shifts every subsequent name by two. The feature table is the authority and the original curation was correct"]},
]

E += [
{"element_id":"R27","element_type":"plasmid","inc_group":"IncHI1","mpf_class":"MPF_F chimera",
 "host":"Salmonella enterica ser. Typhi","size_bp":180461,
 "accessions":{"primary":"AF250878.1","crosscheck":"NC_002305.1","parts":[],"status":"PROVISIONAL"},
 "exclusion":{"gene":"eexA","type":"entry","evidence_tier":"E3","synonym":"trhZ",
   "evidence_tier_note":"E3 not E1/E2: the coordinates are HOMOLOGY-INFERRED plus positional corroboration, NOT extracted from gunton2008. The phenotype is E2 in the literature but this ROW's sequence provenance is inference.",
   "homolog_source":"A0A0H2XLF7 (TrhZ, 129 aa) from pAPEC-O1-R = NC_009838.1 / DQ517526.1, 241,387 bp. PlasmidFinder (local PLSDB): IncHI2A + IncHI2. R27 is IncHI1A + IncHI1B(R27) -- so this is a CROSS-Inc-group comparison inside the IncH family, and the 78.1% identity must be read as such. The positional corroboration (R0015-R0018 matching the curated Z-operon ORF numbering) is independent of the homology and is the stronger half of the evidence.",
   "localisation":"inner membrane","specificity_determinant":"not mapped",
   "resolved_orf":"R0018 / AAF69857.1, 130 aa, AF250878.1:22505..22897 (+)",
   "resolution_method":"UniProt homolog anchoring, not the paper. AF250878.1 carries NO trhZ/eexA name (locus tags R0001-R0210 only). Searched the record with UniProt TrhZ homolog A0A0H2XLF7 (129 aa, pAPEC-O1-R): best hit R0018 at 78.1% local identity, score 533 vs 52 for the runner-up. Corroborated positionally: R0015/R0016/R0017/R0018 correspond exactly to the curated Z operon trhO(015)-orf016(eexB)-orf017-trhZ(eexA), including the ORF numbering.",
   "note":"eexB = R0016 / AAF69855.1, 279 aa (orf016, outer membrane, SURFACE exclusion) -- out of scope under the entry-only cut, recorded here because the homolog search located it too"},
 "partner":{"gene":"trhG","family":"VirB6","virb_equivalent":"VirB6 (F-TraG homolog)","status":"predicted",
   "evidence_tier":"E3",
   "homolog_source":"A0A0H2XKG9 (TrhG, 1317 aa) from pAPEC-O1-R (IncHI2), same cross-group caveat as eexA. Corroborated by PF07916 TraG_N hit at score 72.1 and by sitting 9 bp after trhH.",
   "located_in":"Tra1, F operon","specificity_residues":None,
   "resolved_orf":"R0128 / AAF69966.1, 1329 aa, AF250878.1:113294..117283 (+)",
   "resolution_method":"UniProt homolog anchoring with A0A0H2XKG9 (TrhG, 1317 aa): 76.1% local identity, score 5240 vs 69 for the runner-up. Corroborated positionally -- it sits 9 bp downstream of trhH (AAF69965.1), matching the curated Tra1 F operon trhF-trhH-trhG."},
 "context":{"operon":{"name":"Tra2 Z operon","genes":["trhO (orf015)","orf016 (eexB)","orf017","trhZ (eexA)"],
   "orientation":"-","note":"REVERSE relative to the other two Tra2 operons"},
  "tra1":{"size_orfs":14,"operons":{"H":["traJ","traG","traI","traH"],"R":["trhR","trhY","trhX"],"F":["trhF","trhH","trhG"]},
   "note":"oriT lies between traH and trhR. TraG here is the T4CP (VirD4 class), interacts with TrhB -- NOT the VirB6 homolog; that is TrhG"},
  "tra2":{"size_kb":36,"n_orfs":28,
   "operons":{"AC":["trhA","trhL","trhE","trhK","orf030","trhB","orf028","orf027","trhV","trhC"],
     "Z":["trhO","orf016","orf017","trhZ"],
     "AN":["htdA","htdF","htdK","orf009","trhP","trhW","trhU","trhN"]},
   "orientations":{"AC":"+","Z":"-","AN":"+"}},
  "components":{"trhA":"pilin","trhK":"secretin-like (VirB9/TraK)","trhC":"VirB4 ATPase","trhN":"TraN homolog (Mps)","trhW":"N-terminus is a fusion of F TrbC"},
  "regulation":"TrhR/TrhY activate, HtdA represses; Hha/H-NS thermoregulation (optimal 30 C, suppressed at 37 C); cAMP/CRP and growth-phase dependent"},
 "notes":["THE COUNTEREXAMPLE. eexA sits in Tra2 on the REVERSE strand; its predicted target TrhG is in Tra1; TrhN (the eexB target) is in the Tra2 AN operon. None of the three are adjacent. This is the most complete refutation available of the adjacency-plus-same-strand rule",
  "Temperature MUST be recorded for any IncHI1 exclusion index -- conjugation is thermosensitive"]},

{"element_id":"IncC","element_type":"plasmid","inc_group":"IncC (also IncA)","mpf_class":"MPF_F",
 "host":"Vibrio cholerae / Escherichia coli","size_bp":171829,
 "accessions":{"primary":"CP033514.1","crosscheck":"NZ_CP033514.1","parts":[],"status":"PROVISIONAL",
   "represents":["pVCR94","R16a"]},
 "exclusion":{"gene":"eexC","type":"entry","evidence_tier":"E1",
   "specificity_determinant":"C-terminus","expression":"AcaCD-regulated, low"},
 "partner":{"gene":"traG_C","family":"VirB6","virb_equivalent":"VirB6","status":"verified",
   "specificity_residues":"not resolved to residue level"},
 "context":{"operon":{"name":"traF-traH-traG operon (AcaCD dependent)","genes":["traF","traH","traG_C"],
   "orientation":"+","note":"eexC lies adjacent to traG_C"},
  "tra_operons":[["traL","traE","traK","traB"],["traV","traA"],
    ["dsbC","traC","trhF","traW","traU","traN"],["traF","traH","traG_C"]],
  "dna_processing":{"t4cp":"traD","relaxase":"traI","other":"mobI"},
  "regulation":"acaC-acaD is the master activator (AcaCD binds and activates 18 promoters); acaB activates the acaDC promoter (2020 TraDIS)",
  "note":"traN has its own additional AcaCD-dependent promoter"},
 "notes":["SGI1 encodes traNS/traHS/traGS which replace IncC TraNC/TraHC/TraGC wholesale. TraN-TraH-TraG is therefore an interchangeable mating-pore core, which is also why exclusion targets both TraN and TraG",
  "TraDIS identified 27 conjugation genes (Hancock et al. 2020 Nat Microbiol)",
  "sfx (surface exclusion, targets TraN_C, constitutive at ~150x eexC) is OUT OF SCOPE under the entry-only cut"]},

{"element_id":"SXT","element_type":"ICE","inc_group":None,"mpf_class":"MPF_F",
 "host":"Vibrio cholerae","size_bp":104025,
 "accessions":{"primary":"KJ817376.1","crosscheck":None,"parts":[],"status":"PROVISIONAL",
   "represents":["SXT","R391"],
   "note":"ICDC-1307. Phenotype and the 606-608 mapping derive from SXT-MO10, for which only the ARG cluster (AY034138.1) exists"},
 "exclusion":{"gene":"eex","type":"entry","evidence_tier":"E1","gene_name_in_record":"EexR1",
   "specificity_determinant":"C-terminal 56 aa","compartment":"cytoplasmic"},
 "partner":{"gene":"traG","family":"VirB6","virb_equivalent":"VirB6","status":"verified",
   "specificity_residues":"aa 606-607-608. SXT (S group) P-G-E; R391 (R group) T-G-D",
   "compartment":"cytoplasmic",
   "note":"DISCREPANCY: the reference protocol records the R-group triplet as T-D-D, the curated table as T-G-D. Position 607 is unresolved and must be settled from R391 sequence before any S/R group assignment is published"},
 "context":{"operon":{"name":"P_L regulatory transcript (eex is NOT in a tra operon)",
   "genes":["croS","setD","setC","...","eex"],"orientation":"convergent with traG",
   "note":"in R391 the eex gene is separated from the P_L module by the mer operon. eex is not SetCD-regulated"},
  "tra_operons":[["traL","traE","traK","traB"],["traV","traA"],
    ["s054 (dsbC)","traC","trhF","traW","traU","traN"],["traF","traH","traG"]],
  "note":"tra operon synteny is identical to IncC",
  "other_required":["traI (relaxase)","traD (coupling)","traJ = s043 (coupling)","s063","mobI"],
  "integration":["int","xis"],
  "regulation":["setR (repressor, lambda-CI-like)","setC","setD (FlhD2C2-like activator)","croS"],
  "hotspots":{"HS1":"s043-traL","HS2":"traA-s054","HS3":"s073-traF"}},
 "notes":["VERIFIED 2026-09-01: TraG residues 606-608 read P-G-E in KJ817376.1, confirming ICDC-1307 is an S-group element despite the EexR1 gene name",
  "Both specificity determinants are CYTOPLASMIC. A naive periplasmic-interface model will be wrong -- store predicted topology and interface compartment separately"]},

{"element_id":"ICEBs1","element_type":"ICE","inc_group":None,"mpf_class":"MPF_FA",
 "host":"Bacillus subtilis subsp. subtilis 168","size_bp":20000,
 "accessions":{"primary":"CP171645.1","crosscheck":None,"parts":[],"status":"CHROMOSOMAL",
   "note":"integrated in the 4.2 Mb host chromosome; element extracted by coordinates"},
 "exclusion":{"gene":"yddJ","type":"entry","evidence_tier":"E1",
   "fold":"cystatin-like, DUF4467","specificity_determinant":"not mapped"},
 "partner":{"gene":"conG","family":"VirB6","virb_equivalent":"VirB6-like, 7 TM","status":"verified",
   "synonym":"yddG",
   "specificity_residues":"extracellular loop aa 276-295 is sufficient to confer specificity; resistance mutation E288K falls inside it"},
 "context":{"operon":{"name":"P_xis large operon",
   "genes":["conB (yddB)","conC (yddC)","conD (yddD)","conE (yddE)","yddF","conG (yddG)","cwlT (yddH)","yddI","yddJ","yddK","rapI","phrI","yddM"],
   "orientation":"+"},
  "upstream_module":{"genes":["int","xis","immR","immA","ydzL","ydcO","helP","nicK","ydcS","ydcT","conQ (ydcQ)","yddA"],
    "note":"internal order to be taken from the ICEBs1 annotation, not from this list"},
  "virb_mapping":{"conB":"VirB8-like","conD":"VirB3-like","conE":"VirB4 ATPase","conG":"VirB6-like (7 TM)","cwlT":"VirB1-like (bifunctional cell-wall hydrolase)","conQ":"VirD4 T4CP"},
  "essentiality":{"required":["conB","conC","conD","conE","conG"],"not_required":["yddF"]},
  "interactions":"BACTH: ConE interacts with itself, ConB and ConQ only -- NOT with ConG",
  "regulation":"P_xis operon; RapI-PhrI quorum sensing + ImmR/ImmA + RecA/SOS"},
 "notes":["Gram-positive. Only element in the set with no periplasm, so the compartment logic differs from every Gram-negative row"]},

{"element_id":"pKPC_UVA01","element_type":"plasmid","inc_group":"IncP/IncN-like","mpf_class":"MPF_T",
 "host":"Klebsiella pneumoniae CAV1016","size_bp":43621,
 "accessions":{"primary":"CP009465.1","crosscheck":"NZ_CP009465.1","parts":["CP017937.1"],
   "status":"PROVISIONAL",
   "note":"two live duplicate records at identical length; the tie cannot be broken by aa_sha256 because the trbK-like gene is absent from BOTH annotations"},
 "exclusion":{"gene":"trbK-like","type":"entry","evidence_tier":"E2",
   "exclusion_index":"EI >300, ~400, 200-400 fold, max >2800 (kamruzzaman2022 as retrieved)",
   "protein_length_reported_aa":151,
   "resolved_orf":"AJE44581.1 (152 aa, CP009465.1:21029..21487, 47 bp upstream of trbL-like)",
   "resolution_evidence":"kamruzzaman2022 states 151 aa, explicitly TWICE RP4 TrbK's 69 aa. That excludes the 63-aa lipobox candidate AJE44587.1 outright. AJE44581.1 is 152 aa (off by one, most likely initiator-Met counting) at the operon position the paper's trbJ-trbK-trbL model requires.",
   "resolution_caveat":"UNCONFIRMED. The paper gives NO locus tag and NO coordinates. It describes a lipoprotein signal, but AJE44581.1 has NO lipobox cysteine -- a comprehensive scan of BOTH pKPC records found no ORF with a lipobox AND Asp at +2. Cleaving AJE44581.1 after ...SANA| does put Asp at +2 of the mature peptide, consistent with the paper, but via signal peptidase I rather than lipoprotein processing. Do not promote above E4 without the locus tag.",
   "requirement":"required in BOTH donor and recipient -- the only dual-requirement system in the set",
   "specificity_determinant":"not mapped"},
 "partner":{"gene":"trbL-like","family":"VirB6","virb_equivalent":"VirB6","status":"predicted",
   "specificity_residues":None},
 "context":{"operon":{"name":"IncP-type trb cluster","genes":["trbJ-like","trbK-like","trbL-like"],
   "orientation":"+","note":"structurally isomorphic to RP4 Tra2. Full gene list and coordinates are not given in the literature and must be computed from the accession, ideally by aligning against the RP4 trbA-trbP frame"}},
 "notes":["Highest exclusion index in the set (EI ~1000). Dual donor+recipient requirement is mechanistically distinct from every other row"]},
]

E += [
{"element_id":"R100","element_type":"plasmid","inc_group":"IncFII","mpf_class":"MPF_F",
 "host":"Shigella flexneri 2b","size_bp":94281,
 "accessions":{"primary":"AP000342.1","crosscheck":"NC_002134.1","parts":[],"status":"PINNED"},
 "exclusion":{"gene":"traS","type":"entry","evidence_tier":"E1",
   "specificity_determinant":"not mapped"},
 "partner":{"gene":"traG","family":"VirB6","virb_equivalent":"VirB6","status":"verified",
   "specificity_residues":"TraG 610-673 in F numbering; R100 equivalent must be projected through the alignment, NOT copied as a bare range",
   "specificity_measured":"F vs R100 TraG: 92.8% identical overall, 51.6% over the window, 6.9x substitution enrichment (p=4.2e-19)"},
 "context":{"operon":{"name":"tra operon (P_Y) -- same architecture as F",
   "genes":["traY","traA","traL","traE","traK","traB","traP","traV","traC","traW","traU","trbC","traN","traF","trbB","traH","traG","traS","traT","traD","traI","traX"],
   "orientation":"+","note":"gene order around traG-traS-traT is identical to F"}},
 "notes":["THE POSITION-CONSERVED / SEQUENCE-DIVERGED CONTROL, argued from HMM detection rather than percent identity. The Pfam family PF10624 ('Plasmid conjugative transfer entry exclusion protein TraS'), built on the IncFII type, scores R100 traS at 312.5 at the GA threshold and DOES NOT DETECT F traS at all. A purpose-built family for this exact gene cannot cross from IncFII to IncFI. Yet both sit in the same operon slot immediately downstream of traG. That is a threshold-anchored, reproducible statement and it is the cleanest single argument for the position-first strategy.",
  "DO NOT use the F/R100 TraS percent identity for this argument. The alignment is marginal -- 74 of 203 columns gapped (36%) -- and the value is unstable to method: 17.7% (BLOSUM62 -11/-1, id/alignment length), 27.9% (same alignment, id/ungapped columns), 25.5% (BLOSUM45), 40.0% (local, over a 15-residue block). Both original figures were correct arithmetic under different denominators (see scripts/lib/metrics.py, A6), but no single number here is a stable claim",
  "PGAP re-annotation of NC_002134.1 DROPPED the /gene=\"traS\" qualifier that AP000342.1 carries -- the origin of the GenBank-primary policy in scripts/lib/ncbi.py"]},
]

# Adjacency, split into the TWO claims that were previously conflated.
#
#   proximity  -- is the exclusion gene within S4 step 2's discovery window
#                 (+/-5 genes / +/-6 kb, BOTH STRANDS)? This is the operational
#                 rule the pipeline actually implements.
#   synteny    -- same strand AND same operon? This is the stronger structural
#                 claim in the Garcillan-Barcia & de la Cruz observation.
#
# SXT is why the split matters: its eex is 32 bp from traG and so is trivially
# recovered by the discovery window, but it is transcribed CONVERGENTLY from a
# regulatory transcript and is not in the tra operon at all. Scoring it as a
# single boolean either wrongly fails the discovery heuristic or wrongly claims
# operon membership.
ADJ={
 "F":         {"same_operon":True, "same_strand":True,  "genes_between":0},
 "R100":      {"same_operon":True, "same_strand":True,  "genes_between":0},
 "R64":       {"same_operon":True, "same_strand":True,  "genes_between":1,
               "note":"excA-excB-traY consecutive; excB is the overlapping reinitiation ORF"},
 "R621a":     {"same_operon":True, "same_strand":True,  "genes_between":0},
 "RP4":       {"same_operon":True, "same_strand":True,  "genes_between":0},
 "pKM101":    {"same_operon":True, "same_strand":True,  "genes_between":0,
               "note":"eex sits between the VirB5 (traC) and VirB6 (traD) homologs"},
 "R27":       {"same_operon":False,"same_strand":True,  "genes_between":None,
               "proximity_override":False,
               "note":"eexA in the Tra2 Z operon (reverse strand); target TrhG in the Tra1 F operon -- a different transfer region entirely. NOW MEASURED (2026-09-02): eexA 22505..22897 and trhG 113294..117283 are 90,396 bp apart on the SAME strand -- 15x outside the +/-6 kb discovery window. Note the correction: they ARE co-oriented, so the rule fails on DISTANCE and operon membership, not on strand. The quantitative counterexample the set needed."},
 "IncC":      {"same_operon":True, "same_strand":True,  "genes_between":0,
               "note":"adjacency is asserted by the literature (eexC abuts traG_C) but is UNMEASURED here: eexC is not annotated in either representative. An earlier automated measurement of 128,521 bp was a false positive -- the generic phrase \"entry exclusion\" had matched R16a /gene=\"trbK\" at 35.7 kb against a TraG at 164 kb"},
 "SXT":       {"same_operon":False,"same_strand":False, "genes_between":0,
               "note":"eex is in the P_L regulatory transcript, transcribed CONVERGENTLY with traG. Not in the tra operon, not same-strand -- but only 32 bp away, so the discovery window finds it regardless"},
 "ICEBs1":    {"same_operon":True, "same_strand":True,  "genes_between":2,
               "note":"conG and yddJ separated by cwlT and yddI inside the same P_xis operon"},
 "pKPC_UVA01":{"same_operon":True, "same_strand":True,  "genes_between":0,
               "proximity_override":True,
               "note":"by homology to RP4 Tra2; gene absent from both annotations so unmeasured"},
}
DISCOVERY_WINDOW_BP=6000
DOMAIN_ID={'F': 'PF07916 (TraG_N)', 'R100': 'PF07916 (TraG_N)', 'IncC': 'PF07916 (TraG_N)', 'SXT': 'PF07916 (TraG_N)', 'R27': 'PF07916 (TraG_N)', 'RP4': 'PF04610 (TrbL)', 'pKM101': 'PF04610 (TrbL)', 'pKPC_UVA01': 'PF04610 (TrbL)', 'R64': 'IPR027628 (DotA_TraY) -- no Pfam', 'R621a': 'IPR027628 (DotA_TraY) -- no Pfam', 'ICEBs1': 'CDD:cd06261 (TM_PBP2) / IPR000515 (MetI-like)'}
HMM_HIT={'F': ('no_hit', ''), 'R100': ('HIT', 'TraS PF10624.15 score=312.5'), 'ICEBs1': ('HIT', 'DUF4467 PF14729.13 score=98.9'), 'R64': ('HIT', 'surf_exc_IncI1 NF033891.1 score=289.7'), 'R621a': ('HIT', 'surf_exc_IncI1 NF033891.1 score=216.4'), 'RP4': ('HIT', 'TrbK_RP4 TIGR04359.1 score=75.6'), 'SXT': ('HIT', 'EexR NF041429.1 score=274.5'), 'pKM101': ('HIT', 'Eex_IncN NF033894.1 score=70.5'), 'R27': ('no_hit', ''), 'IncC': ('no_seq', ''), 'pKPC_UVA01': ('no_seq', '')}

# authoritative = the coordinate-system anchor (exactly one per element)
# mirror        = same sequence, same coordinates (RefSeq <-> GenBank)
# related       = DIFFERENT coordinates; comparison only, NEVER coordinate validation
RECORD_SCOPE={"F":"plasmid","R100":"plasmid","R64":"plasmid","R621a":"plasmid","RP4":"plasmid",
 "pKM101":"fragment","R27":"plasmid","IncC":"plasmid","SXT":"ICE",
 "ICEBs1":"chromosome","pKPC_UVA01":"plasmid"}
MIN_BP_BY_SCOPE={"chromosome":1000000,"plasmid":30000,"ICE":15000,"fragment":1000}


from Bio import SeqIO
import glob as _glob

# an element row may consolidate several sequenced representatives, and the
# frozen records are filed under the representative's name, not the row's
RECORD_ALIAS={"IncC":["pVCR94","R16a"],"SXT":["SXT","R391"]}
# gene names as they actually appear in the records, where they differ
GENE_ALIAS={("SXT","eex"):["EexR1","eex"],
            # R27 carries only systematic locus tags; these were resolved by
            # UniProt-homolog anchoring, see the element notes
            ("R27","eexA"):["R0018"],("R27","trhG"):["R0128"],("ICEBs1","conG"):["conG","yddG"],
            ("ICEBs1","yddJ"):["yddJ"],("IncC","traG_C"):["traG"],("IncC","eexC"):["eexC"]}

def _locate(eid, gene, only_prefix=None):
    """Find a gene's span by scanning the element's frozen GenBank records.
    Matches /gene exactly, else the bare name inside /product. Returns
    (start, end, strand, accession) or None."""
    names=GENE_ALIAS.get((eid,gene)) or [gene]
    bases=[n.split("-")[0].split("_")[0] for n in names]
    paths=[]
    for pref in ([only_prefix] if only_prefix else RECORD_ALIAS.get(eid,[eid])):
        paths+=sorted(_glob.glob(os.path.join(PROJ,"data","seed","genbank",f"{pref}__*.gb")))
    for path in paths:
        try: rec=next(SeqIO.parse(path,"genbank"))
        except Exception: continue
        for f in rec.features:
            if f.type!="CDS": continue
            g=((f.qualifiers.get("gene") or [""])[0]).strip()
            prod=(f.qualifiers.get("product") or [""])[0]
            if any(g.lower()==n.lower() for n in names+bases) or \
               any(re.search(rf"\b{re.escape(b)}\b",prod,re.I) for b in bases):
                return (int(f.location.start)+1,int(f.location.end),f.location.strand,rec.id)
    return None

def measure(eid, egene, pgene):
    """bp between the exclusion gene and its target.

    Resolution is per REPRESENTATIVE and both genes must resolve within the same
    one: an IncC row consolidates pVCR94 and R16a, and coordinates from two
    different plasmids are not comparable. Within a representative, NZ_X and X
    mirror the same sequence and ARE comparable.

    Uses the corrected partner assignment, so it is independent of the stage-3
    extraction table (built before pKM101's partner was fixed from traG/VirB11
    to traD/VirB6).
    """
    norm=lambda acc: re.sub(r"^(NZ_|NC_)","",acc)

    def from_table(rep,g):
        r=seq.get((rep,g))
        if not r: return None
        acc,rng=r["coords"].split(":"); lo,hi=rng.split("..")
        return (int(lo),int(hi),1 if r["strand"]=="+" else -1,acc)

    for rep in RECORD_ALIAS.get(eid,[eid]):
        a=_locate(eid,egene,rep) or from_table(rep,egene)
        b=_locate(eid,pgene,rep) or from_table(rep,pgene)
        if a and b and norm(a[3])==norm(b[3]):
            gap=(b[0]-a[1]-1) if b[0]>a[1] else (a[0]-b[1]-1)
            return gap, a[3]
    return None, None

def main():
    os.makedirs(OUT,exist_ok=True)
    for el in E:
        eid=el["element_id"]
        egene=el["exclusion"]["gene"]; pgene=el["partner"]["gene"]
        el["exclusion"]["sequence"]=S(eid,egene) or {"status":"NOT_RECOVERED",
            "needs":"manual coordinates from the primary paper"}
        el["partner"]["sequence"]=S(eid,pgene) or {"status":"NOT_RECOVERED",
            "needs":"manual coordinates from the primary paper"}
        adj=dict(ADJ.get(eid,{}))
        bp,bp_acc=measure(eid,egene,pgene)
        ov=adj.pop("proximity_override",None)
        # tri-state: True / False / None. None means UNMEASURED -- distinct from
        # a measured failure. Collapsing them would let a gap in the annotation
        # masquerade as evidence against the heuristic.
        if ov is not None:        within=ov
        elif bp is None:          within=None
        else:                     within=abs(bp)<=DISCOVERY_WINDOW_BP
        el["context"]["adjacency"]={
            "proximity":{"measured_bp":bp,
                         "measured_from":bp_acc,
                         "within_discovery_window":within,
                         "window_bp":DISCOVERY_WINDOW_BP},
            "synteny":{"same_strand":adj.get("same_strand"),
                       "same_operon":adj.get("same_operon"),
                       "genes_between":adj.get("genes_between")},
            "supports_discovery_heuristic":within,
            "supports_synteny_claim":bool(adj.get("same_strand") and adj.get("same_operon")),
            "note":adj.get("note")}
        scope=RECORD_SCOPE.get(eid,"plasmid")
        el["accessions"]["record_scope"]=scope
        el["accessions"]["min_expected_bp"]=MIN_BP_BY_SCOPE[scope]
        if el.get("inc_group") in (None,"NA",""):
            el["inc_group"]="not_applicable"      # ICEs have no Inc group; not missing
        # partner_family is the MECHANISTIC label, sourced from the literature.
        # partner_domain_id is the SEARCHABLE anchor, sourced from InterPro/Pfam.
        # These are different things and both are kept: an automatic CDS-level
        # domain call must not be allowed to overturn an experimental assignment
        # (the same error avoided on RP4 TrbK's signal peptide).
        fam=el["partner"].get("family")
        if fam and fam!="VirB6":
            el["partner"]["family_note"]=fam
            el["partner"]["family"]="unassigned"
        el["partner"]["domain_id"]=DOMAIN_ID.get(eid,"")
        hh,hd=HMM_HIT.get(eid,("unknown",""))
        el["exclusion"]["hmm_hit"]=hh
        el["exclusion"]["hmm_model"]=hd
        el["exclusion"]["benchmark_role"]=(
            "positive_control -- already recoverable by an existing public HMM"
            if hh=="HIT" else
            "DE NOVO TARGET -- no existing HMM family finds this gene"
            if hh=="no_hit" else "pending -- no sequence yet")
        el["schema_version"]="1.1.0"
        el["scope"]="entry_exclusion_only"
        hdr=(f"# {eid} -- entry exclusion pair record\n"
             "# Curated context + sequence facts verified from the pinned accession.\n"
             "# THIS FILE IS THE SOURCE OF TRUTH. data/exports/*.tsv|json are generated.\n")
        body=yaml.safe_dump(el,sort_keys=False,allow_unicode=True,default_flow_style=False,width=100)
        open(os.path.join(OUT,f"{eid}.yaml"),"w").write(hdr+body)
    print(f"wrote {len(E)} element files -> {OUT}")
    print(f"\n{'element':<11} {'exclusion':<10} {'target':<10} {'VirB6?':<7} {'bp':>7} {'discov':<7} {'synteny':<7}")
    print("-"*64)
    for el in E:
        a=el["context"]["adjacency"]; bp=a["proximity"]["measured_bp"]
        print(f"{el['element_id']:<11} {el['exclusion']['gene']:<10} {el['partner']['gene']:<10} "
              f"{('yes' if el['partner']['family']=='VirB6' else 'NO'):<7} "
              f"{(str(bp) if bp is not None else '-'):>7} "
              f"{str(a['supports_discovery_heuristic']):<7} {str(a['supports_synteny_claim']):<7}")

if __name__=="__main__": main()
