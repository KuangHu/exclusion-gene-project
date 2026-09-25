# Gold-Standard Collection Protocol: A Curated Database of Conjugation Exclusion Genes and Their Cognate Partners

## 0. Core design decision

**The unit of record is the pair, not the gene.**

Exclusion has no meaning as a single-gene annotation. What is biologically real — and what defines an exclusion group — is a two-sided recognition event:

| Side | Location | Function | Examples |
|---|---|---|---|
| **Exclusion factor** | Recipient | Blocks transfer | TraS, TraT, TrbK, ExcA, Eex, EexA/EexB, EexC, Sfx, YddJ, MbeD, Ses, PrgA, Sea1 |
| **Cognate target** | Donor | Recognized element of the T4SS | TraG, TraY, ConG, TraG_C, TraN, **TrbL** (hypothesized) |

So every row in the database is `element × exclusion_gene × partner_gene × specificity_evidence`. A database of exclusion genes alone cannot express exclusion groups, which is the whole point.

The single most important structural insight for collection: **the donor-side target is, in almost every solved case, a VirB6/TraG-family protein**, and the exclusion gene is very often its immediate genomic neighbor. TrbK/TrbL in RP4 is the cleanest instance of this and is the basis of the main search heuristic (§4).

---

## 1. Schema

Minimum viable table. One row per pair; use `NULL` + evidence tier rather than guessing.

```
element_id              # F, R100, RP4, R64, pKM101, R27, pVCR94, SXT, R391, ICEBs1, pLS20, pCF10, pAD1, ColE1, pKPC_UVA01
element_type            # plasmid | ICE | mobilizable
inc_group               # IncFI, IncFII, IncI1, IncN, IncP-alpha, IncW, IncHI1, IncA, IncC, ...
mob_class               # MOBF, MOBP, MOBH, MOBQ, MOBC (Garcillán-Barcia MOB typing)
t4ss_type               # T4SS-typeF, -typeT, -typeG, -typeI, -typeB (CONJscan)

exclusion_gene          # traS, traT, trbK, excA, excB, eex, eexA, eexB, eexC, sfx, yddJ, mbeD, ses, prgA, sea1
exclusion_type          # entry | surface | both | undetermined
protein_length_aa       # precursor / mature if cleaved
localization_pred       # DeepTMHMM + LipoP call
localization_exp        # fractionation result, if published
lipoprotein             # Y/N + lipobox sequence
protein_accession       # RefSeq WP_ / GenBank, WITH version suffix
nuc_accession_coords    # accession:start..end, strand

partner_gene            # traG, traY, conG, traN, trbL(?), NULL
partner_family          # VirB6 | VirB5 | MPS-adhesin | unknown
partner_accession
adjacency_bp            # distance between exclusion gene and nearest VirB6-family gene
specificity_residues    # e.g. TraG 606-608 (PGE/TDD); TraT 116-120; ConG 276-295 (E288); Eex C-term
exclusion_group_label   # S/R (SXT), C/D/E (IncC), Sfx I-IV (F-like), R64/R621a (IncI1)

exclusion_index         # fold reduction, numeric
index_assay             # plate | liquid | single-cell fluorescence
index_conditions        # temp, mating time, D:R ratio, in-trans vs native expression
evidence_tier           # E1..E4 (below)
pmid / doi
curator, curation_date, schema_version
```

### Evidence tiers (non-negotiable — this is what makes it "gold standard")

- **E1** — Exclusion demonstrated genetically *and* cognate partner identified *and* specificity residues mapped. → F *traS*/TraG; SXT *eex*/TraG; ICEBs1 *yddJ*/ConG; IncI1 *excA*/TraY; IncC *eexC*/TraG_C; IncC *sfx*/TraN.
- **E2** — Exclusion gene identified and validated, partner unknown. → RP4 *trbK*; pKM101 *eex*; R27 *eexA*/*eexB*; ColE1 *mbeD*; pLS20 *ses*; pAD1 *sea1*; F *traT* (target disputed).
- **E3** — Exclusion phenotype mapped to a region, no gene isolated. → IncW R388 (PILW region only).
- **E4** — Homology/context prediction only, no phenotype. → everything your pipeline generates in §4.

Never mix E4 rows into analyses that report exclusion-group structure. Keep them in a separate `candidates` table.

---

## 2. Seed set (the hand-curated positives)

Build these ~18 rows manually from the primary papers before writing any code. This is your training/benchmark set, and every automated method gets scored against it.

| Element | Exclusion gene(s) | Partner | Tier |
|---|---|---|---|
| F (IncFI) | *traS* / *traT* | TraG / disputed (TraN–OmpA?) | E1 / E2 |
| R100, R1, R6-5, ColB2, pED208 (IncFII) | *traS* / *traT* | TraG | E1 / E2 |
| R64, R144 (IncI1); R621a (IncIγ) | *excA*, *excB* | TraY | E1 |
| pKM101 (IncN) | *eex* | unknown (adjacent VirB6-family gene) | E2 |
| RP4 (IncPα) | *trbK* | unknown — **adjacent gene is *trbL* (VirB6)** | E2 |
| R388 (IncW) | none named | unknown (TrwI = VirB6) | E3 |
| R27 (IncHI1) | *eexA* (entry, IM), *eexB* (surface, OM) | each other | E2 |
| pVCR94 / R16a (IncC, IncA) | *eexC*, *sfx* | TraG_C, TraN | E1 |
| ColE1 | *mbeD* | unknown | E2 |
| pKPC_UVA01 | *trbK*-like (dual donor+recipient) | self | E2 |
| SXT, R391 (ICE) | *eex* | TraG | E1 |
| ICEBs1 | *yddJ* | ConG (VirB6) | E1 |
| pLS20 (B. subtilis) | *ses* | unknown | E2 |
| pCF10 (E. faecalis) | *prgA* | steric, no partner | E2 |
| pAD1 | *sea1* | unknown | E2 |
| pAM373 | *sea1* absent | — | negative control |

**Explicit negatives matter as much as positives.** Reserve rows for: pAM373 (*sea1*-less pheromone plasmid), RSF1010/IncQ (mobilizable, no self-exclusion), ColE1 *exc1*/*exc2* (annotated as exclusion, experimentally shown inactive), and virulence T4SSs (Ti VirB/D4, *Helicobacter* Cag, *Legionella* Dot/Icm — systematically exclusion-free). Without negatives you cannot compute precision in §5.

---

## 3. Sequence space and retrieval

**Do not** pull sequences from memory or by editing accession numbers. For each seed element, get the accession from the primary paper's data-availability statement, then fetch the versioned record from NCBI Nucleotide, and record `accession.version`.

Search space for expansion, in order of preference:

1. **PLSDB** — curated, deduplicated complete bacterial plasmids with metadata. Primary search space.
2. **COMPASS** — complete plasmid set with taxonomic/curation annotation; useful cross-check.
3. **ICEberg 3.0** — curated ICEs and IMEs; the only sane starting point for SXT/R391, ICEBs1, ICESt1/3, ICEclc.
4. **RefSeq complete genomes** — for chromosomally integrated ICEs and for chromosomal exclusion-gene homologs (recent TraT surveys found broad divergent chromosomal lineages, so do not restrict to plasmids).

Deduplicate at ≥99% ANI before any frequency or phylogenetic claim, or plasmid-sequencing bias will manufacture "exclusion group" structure that does not exist.

---

## 4. The pipeline — genomic context is the discovery engine

Sequence homology alone fails on exclusion genes. They are small (47–250 aa), fast-evolving, and most have no Pfam family. The method that actually works is **anchor on the machinery, then mine the neighborhood.**

### Step 1 — Type the conjugation machinery
- **CONJscan / MacSyFinder** — this is the reference standard. It gives you T4SS type, MOB class, and, critically, **coordinates of the VirB6-family gene** (*trbL*, *traG*, *traY*, *conG*, *trwI*).
- **oriTfinder / oriTDB** — *oriT*, relaxase, T4CP; confirms the element is genuinely conjugative rather than a decayed remnant.
- **MOB-suite / MOB-typer**, **PlasmidFinder + pMLST** — replicon and Inc typing. Non-optional, because every exclusion group in the literature is indexed against an Inc group.

Reject elements failing CONJscan's conjugative-system criteria before neighborhood mining. Exclusion genes in non-conjugative contexts are almost always false positives.

### Step 2 — Extract the neighborhood
For each VirB6-family hit, extract ORFs within ±5 genes / ±6 kb, both strands.

### Step 3 — Filter for exclusion-gene phenotype signature
Flag ORFs meeting:
- length 40–260 aa (covers TrbK's 47-aa mature form through TraT's ~245);
- **AND** either a lipobox (**LipoP / PRED-LIPO / SignalP 6.0**) or 1–2 predicted TM helices (**DeepTMHMM**);
- **AND** no confident functional assignment from **InterProScan** (Pfam, TIGRFAM, CDD);
- **AND** short intergenic distance to the VirB6 gene, often translationally coupled (pKM101 *eex* is coupled to its neighbors — check for overlapping ORFs and Shine-Dalgarno within the upstream CDS, and for translational reinitiation as in R64 *excAB*).

This heuristic is not invented here — it is the operational form of the Garcillán-Barcia & de la Cruz (2008) observation that in IncP, IncW, and IncN the exclusion gene abuts a VirB6 homolog. RP4 *trbK*–*trbL* is the validating case: it is recovered by this rule with no prior knowledge.

### Step 4 — Structure-based homology (the modern addition)
Sequence-level HMM search saturates fast. Run:
- **AlphaFold2/3** (or ESMFold for throughput) on all flagged ORFs;
- **Foldseek** against AFDB + PDB and against your own seed-set models.

This recovers remote homologs that BLAST/HMMER miss, and it is how you can plausibly assign families to the newly described **Sfx** (IncC) and **Ses** (pLS20) proteins, which currently have essentially no sequence-family representation. YddJ's cystatin-like DUF4467 fold is the proof of concept that fold assignment carries information here.

### Step 5 — Predict the pair, don't assume it
- **AlphaFold-Multimer / Boltz / AF3** on candidate `exclusion protein + VirB6 homolog` dyads. Score with ipTM/pDOCKQ and, importantly, check whether the predicted interface **contains the known specificity segments** (TraG 606–608, ConG 276–295). An interface prediction that ignores the experimentally mapped specificity residues is a red flag, not a result.
- Caveat to record in the metadata: in SXT/R391 the functionally critical Eex and TraG determinants are **cytoplasmic**, so a naive periplasmic-interface model will be wrong. Store predicted topology and interface compartment as separate fields so this is auditable.

### Step 6 — Confounder subtraction (routinely skipped, routinely fatal)
Run **DefenseFinder** and **PADLOC** on every element and host genome. Restriction-modification, CRISPR-Cas, BREX, and Wadjet all reduce measured transfer frequency and will be misattributed to exclusion if you only have a phenotype. Record co-occurring defense systems as a covariate column. Same for **FinO/FinP** fertility-inhibition loci — these reduce donation, not reception, and must not be scored as exclusion.

---

## 5. Pair assignment, exclusion groups, and specificity

1. **Reciprocal co-phylogeny.** Build independent trees for the exclusion protein and its partner across the element family. Congruent clades = exclusion groups. This is exactly how the SXT/R391 S and R groups and the IncC C/D/E groups were assigned, and it is the only defensible way to predict a group without wet-lab work.
2. **Variable-segment extraction.** Within each family, align and locate the low-identity window. This is where specificity lives, every single time:
   - TraG_F vs TraG_R100: 93% identical overall but only ~17% over residues **610–673** → TraS recognition.
   - SXT/R391 TraG: residues **606–608**, P-G-E (S group) vs T-D-D (R group).
   - F-like TraT: residues **116–120** (five-residue window).
   - ICEBs1 ConG: residues **276–295**, with E288 pivotal; plus two YddJ regions.
   - IncI1 TraY: internal variable segment (R64) vs C-terminal segment (R621a).
   - IncC Eex: C-terminus.
3. Store specificity windows as coordinate ranges **against a named reference accession**, never as bare numbers. Residue numbering drifts between papers and between precursor/mature forms.

---

## 6. Experimental validation loop — reporting standard

Any row promoted above E4 needs a mating experiment. Minimum reportable set:

- **Exclusion index** = transfer frequency into exclusion-negative recipient ÷ frequency into exclusion-expressing recipient. Report as fold-reduction with confidence interval and *n* biological replicates.
- **Both denominators.** Transconjugants per donor *and* per recipient. Surface exclusion changes aggregate formation, so the two normalizations diverge diagnostically.
- **Plate and liquid matings in parallel.** Liquid mating additionally requires thin/flexible pili — indispensable for IncI1, and the difference between the two formats is itself evidence for surface vs entry exclusion.
- **Isogenic controls:** Δexclusion recipient, exclusion gene supplied *in trans* on a compatible vector, and a vector-only recipient. Complement in the recipient only, then in the donor only — this is how the pKPC_UVA01 dual-requirement system and the recipient-only requirement of TrbK and YddJ were established.
- **Incompatibility and lethality controls.** The ColE1 case is the cautionary tale: resident-plasmid incompatibility and colicin killing both mimic exclusion. Use an exclusion gene cloned away from its replicon.
- **Temperature.** Record it explicitly. IncHI conjugation is thermosensitive (optimal ~30 °C, suppressed at 37 °C); an unrecorded temperature makes an IncHI1 index uninterpretable.
- **Aggregate counts.** Microscopic quantification of stable mating pairs separates "pairs don't form" (surface) from "pairs form, DNA doesn't move" (entry). Without this the `exclusion_type` field is a guess.
- **Discovery-mode methods**, for filling E3 gaps: **TraDIS**/TnSeq on the donor or recipient (this is how IncC *sfx* was found after *eexC* failed to explain residual exclusion), and single-cell fluorescent DNA reporters (e.g. TetR-mNeonGreen on a *tetO*-arrayed incoming plasmid) for directed mutagenesis at scale.

---

## 7. Curation QC

- **Dual independent extraction.** Two curators read each primary paper into the schema; reconcile discrepancies against the figures, not the abstract.
- **Sequence-provenance auditing.** The F vs R100 TraT episode — where a reported F sequence variant was later found to be an error, and where *traN* swaps failed to flip specificity — means historical single-residue specificity claims must be re-derived from current accessions, not copied from the 1990s text. Add a `sequence_verified_against` field with accession.version and date.
- **Nomenclature reconciliation.** Maintain a synonym table: *trhZ*→*eexA*; ORF016→*eexB*; *sec10*→*prgA*; *exc*→*excAB*; ColE1 *exc1*/*exc2* → deprecated, not exclusion. Store the historical name so old literature remains searchable.
- **Versioning.** Semantic-version the release, freeze the accession set per release, and keep a changelog with reason-for-change per row. Deposit in Zenodo with a DOI so the exclusion-group assignments you publish are citable and reproducible.
- **Machine-readable output.** TSV + JSON, plus a GFF3 of exclusion-gene coordinates per element so the pairs can be viewed in a genome browser alongside the *tra* region.

---

## 8. Known gaps worth designing the database to accommodate

- **Missing partners:** RP4 *trbK* (test *trbL* first), pKM101 *eex*, ColE1 *mbeD*, R27, pLS20 *ses*, pAD1 *sea1*.
- **Missing genes entirely:** IncW (R388 region-level only), IncL/M, IncX, IncU, IncHI2.
- **Structurally unresolved:** no high-resolution structure for TraS, any Eex, TrbK, ExcA, EexC, or YddJ. TraG has solution/SAXS characterization only.
- **Mechanistic split:** the TraT cryo-EM work argues against a specific partner interaction for surface exclusion, while entry exclusion is squarely interaction-based. Do not force one mechanistic model onto both `exclusion_type` values.
- **Evasion:** SGI1 remodels the IncC T4SS to escape IncC exclusion. Add a boolean `evasion_documented` field — exclusion escape is its own emerging category.
- **Non-canonical systems:** Streptomyces/actinomycete conjugation uses dsDNA TraB/FtsK translocation with *kil/kor* and *spd* genes, and has no established Eex/Sfx analog. Schema should tolerate `partner_family = N/A` rather than forcing a VirB6 slot.

---

## Suggested build order

1. Hand-curate the 18 seed rows + 4 negative rows from primary papers (§2), with verified accessions.
2. Run CONJscan + PlasmidFinder + oriTfinder over PLSDB and ICEberg; store machinery coordinates.
3. Apply the neighborhood rule (§4 steps 2–3); score recovery of the seed set. Target: recover *trbK*, pKM101 *eex*, *yddJ*, *eexC* without prior knowledge. If it misses these, the filter is wrong, not the biology.
4. Add Foldseek/AF structure layer; re-score.
5. Co-phylogeny → predicted exclusion groups (E4 candidates table).
6. Prioritize wet-lab: the IncW R388 gene and the RP4 TrbK–TrbL hypothesis are the two highest-value, most tractable open questions.
