# Entry-exclusion gene discovery: pipeline and results

Four conjugation classes catalogued from PLSDB 2024_05_31_v2 (72,556 plasmids;
72,438 with called CDS). All figures below are recomputed from the sealed release
files, not quoted from working notes.

---

# 1. Pipeline

## 1.1 Inputs, frozen

| | |
|---|---|
| plasmids | PLSDB 2024_05_31_v2, RefSeq `NZ_`/`NC_` |
| gene calling | Pyrodigal anon mode: `ANON=True, MIN_GENE=90, CLOSED=True, TRANSL_TABLE=11` |
| protein cache | `cds_cache_full` 7,725,419 proteins / 72,438 plasmids |
| profiles | Pfam-A (30,134), CONJScan (404) |
| search | pyhmmer `hmmsearch`, **gathering thresholds only** (A13) |

`cds_cache` (1,099,911 proteins over the 7,436 MPF_T plasmids) was verified
byte-identical to `cds_cache_full` on those accessions — 0 disagreements in
per-accession protein counts.

## 1.2 The shape every class follows

```
entry criterion  ->  architecture  ->  slot  ->  catalogue  ->  post-hoc exclusion families
```

**Entry criterion.** One profile (or conjunction) defines class membership. Chosen
on: seed recall, admission count, and cross-class contamination against the other
catalogues. *Not* chosen on MOB-suite purity — that was measured to be
class-dependent.

**Architecture.** Gene order and layout from anchors only. Never uses the slot
occupant, so the slot cannot be defined circularly.

**Slot.** A position relative to an anchor, validated on seeds and then at class
scale. Offsets are taken in the anchor's **transcription frame**, and same-strand
fraction is reported separately.

**Catalogue.** One row per plasmid, anchors with coordinates, `slot_ready` =
anchor present AND slot occupant resolvable AND architecture callable.

**Exclusion families are POST-HOC.** `TIGR04359`, `NF033894`, `NF041429`,
`NF033891`, `PF10624`, `PF14729` are scanned *after* the catalogue is built and
feed only annotation columns. They never enter an entry criterion, an anchor set,
a length prior, or a slot definition.

## 1.3 Standing rules, enforced in code

| | |
|---|---|
| A9 | curated and generated files never share a path; curated files are SHA256-sealed |
| A10 | a cited accession must resolve to a model actually held |
| A11 | heavy work must not run on a login node |
| A12 | a lookup key must exist in the database it is looked up in |
| A13 | **an anchor threshold may never be an E-value** |

Plus rules carried as practice: anchors defined by VirB number / Pfam family and
never by IncP gene name; `assign_method = positional` never feeds `slot_status`;
`(element_id, gene_name)` as join key, never gene name alone; every benchmark
carries a baseline expected to score perfectly, and nothing else is read if it
fails.

## 1.4 Scripts

| script | role |
|---|---|
| `106_catalogue_v1.py` | MPF_T catalogue |
| `115_catalogue_mpff_v1.py` | MPF_F catalogue |
| `130_catalogue_v1_1.py` | MPF_T v1.1, corrected VirB5 |
| `131_merge_annotation_v1_1.py` | unified MPF_T/F annotation |
| `132_catalogue_mpfi_v1.py` | MPF_I catalogue |
| `133`–`138` | MPF_FA: seed readout, entry criterion, residual identity, slot on seeds, PF14729 census, catalogue |

---

# 2. Results

## 2.1 The four catalogues

| class | file | rows × cols | entry criterion | slot |
|---|---|---|---|---|
| MPF_T | `catalogue_MPF_T_v1.1.tsv` | 7,436 × 36 | `PF03135` (VirB4) | VirB5±1 |
| MPF_F | `catalogue_MPF_F_v1.tsv` | 11,120 × 30 | `PF11130` (TraC_F_IV) | TraG_N+1 |
| MPF_I | `catalogue_MPF_I_v1.tsv` | 3,732 × 46 | `T4SS_I_traU` | TraY+1 |
| MPF_FA | `catalogue_MPF_FA_v1.tsv` | 2,391 × 29 | `virb4 AND NOT(T/F/I) AND ≥2 FA/FATA` | **none** |

**24,470 distinct plasmids**, 33.8% of PLSDB.

## 2.2 Per class

### MPF_T (IncP / trb) — 7,436
* tier A 6,717 (90.3%)
* `slot_ready__virb5_virb6` 6,548 (**88.1%**)
* VirB5 present 6,881 — failure **7.5%**
* VirB5 call source: `pfam_only` 3,468, `both` 2,077, `conjscan_only` **1,336**, none 555
* architecture: canonical 6,020, rearranged 777, uncallable 639
* named family at slot: TrbK / `NF033894`, **224**

### MPF_F (F / tra) — 11,120
* tier A 9,994 (89.9%)
* `slot_ready__traG` 7,803 (70.2%)
* named family at slot: TraS / `PF10624`, **2,534** = 32.5% of `slot_ready`, 22.8% of the class
* layout: contiguous 7,209, split_2 2,463, split_3 1,304, split_4+ 144

### MPF_I (IncI) — 3,732
* `slot_ready__traY` 3,607 (**96.7%**)
* ExcA: **59.3%** of admitted (2,213/3,732) — and **96.5%** of all PLSDB ExcA carriers (2,213/2,293). Different denominators; neither implies the other.
* strata: order#1 1,544, split_2 1,195, order#3 439, order#4+ 305, order#2 132
* `slot_admission = undetermined` on every row; `slot_maxkd` carried as a column (median 2.26), never as a filter
* excB in seven columns; unscanned rows blank, never 0

### MPF_FA (Gram-positive) — 2,391
* subclass: FATA 1,170, unassigned 726, FA 495
* `profile_support`: high 1,041, medium 593, minimal 757 (median 4 profiles, max 16)
* relaxase 1,963 (82.1%)
* **`slot_status = no_slot_defined` on all 2,391 rows**

## 2.3 Unified annotation

`annotation_MPF_T_F_v1.1.tsv` — 18,381 × 53. One row per plasmid, both class
blocks coexisting; `not_admitted` rather than dropped.

**175 dual-system plasmids, 144 tier A in both.** MOB-suite's single-valued
`mpf_type` cannot represent a dual system — it calls 143 of them MPF_F and 31
MPF_T, reporting one and missing the other. This is why the classes are kept as
separate blocks rather than resolved to one label.

## 2.4 Findings that generalised

**The VirB4-family ATPase is the entry criterion in every class.** `T4SS_virb4`
is mandatory in every CONJScan plasmid class definition and fires on all five
Gram-positive seeds (91.0–279.8). MPF_T `PF03135`, MPF_F `PF11130`, MPF_I via
`traU`, MPF_FA `conE`.

**Criteria do not transfer across classes — four measured instances.**

| instance | measurement |
|---|---|
| MPF_T lipobox | 56.1% in-slot vs 2.8% background — 20× |
| MPF_F lipobox | 1.8% |
| MPF_I hydrophobicity | 2.3× against a 40.5% background |
| MPF_FA tiering | no mandatory core exists to tier on |

**A generic profile can fail on an entire class.** `T4SS_virb4` admits 99.99% of
MPF_T and 99.6% of MPF_F but **3.9% of MPF_I** (145/3,732).

**The named-family control does not scale equally.**

Two different quantities; do not mix them.

**Family present anywhere on a plasmid of that class** (census over all 72,437):

| class | family | on-plasmid |
|---|---|---|
| MPF_F | TraS `PF10624` | 2,725 (of 3,011 db-wide) |
| MPF_I | ExcA `NF033891` | 2,213 (of 2,293 db-wide) |
| MPF_T | TrbK `TIGR04359` | 224 (of 263 db-wide) |
| MPF_FA | DUF4467 `PF14729` | **11** (of 13 db-wide) |

**Family at the slot position specifically** — a strictly smaller number:

| class | at-slot | share of `slot_ready` |
|---|---|---|
| MPF_F | TraS 2,534 | 32.5% |

---

# 3. Known limits

## 3.1 MPF_FA has no slot

`conG+3` is confirmed on ICEBs1 alone. Seed check: 1 of 4 locatable seeds — two
occupants **anti-oriented** to the anchor (not plausibly co-transcribed), one an
identified type II toxin-antitoxin system, one seed with no anchor. Class census:
DUF4467 on 11 of 2,391 (0.5%), and within those 11 the modal offset from
`FA_orf15` is **−1**, not ICEBs1's **+3** — the one known case is not the modal
case.

A slot column validated on a single element would be empty, not caveated. It is
omitted.

## 3.2 Withdrawn during this work

**The VirB5 `E≤1e-5` threshold.** Its stated rationale was that GA loses IncP
TrbJ while E-values rescue it. Measured under the pipeline's own search, RP4 TrbJ
is 17.7 bits **E=0.22** and R751 22.2 bits **E=0.0093** — both FAIL `E≤1e-5`. The
rescuing E-values came from a protein-vs-Pfam scan (N≈30,134 models), not a
model-vs-proteome search (N=1,099,911). The rule never rescued its own controls,
and at ~31.7 bits was *stricter* than the GA 24.7 it replaced, costing ~177
plasmids while documented as a loosening. Replaced by
`PF07996 @ GA OR T4SS_T_virB5 @ GA`, both `aa>150`. A13 now forbids E-value
thresholds outright.

**`VirB5 27.1% → 3.1%` and `slot_ready 70.9% → 87.2%`.** Numerator counted over
72,556 plasmids, denominator over 7,436. Correct: **25.4% → 7.5%** and
**70.9% → 88.1%**.

**"The two-tier anchor structure collapses."** VirB2 and VirB3 appear to drop to
5.2% and 5.0% under CONJScan disjunctions, but 33.1% and 57.4% of their new
admissions are MPF_F. Both disjunctions rejected; GA figures stand. VirB1 passed
both controls and stayed peripheral at 13.2% — the negative control that made the
rejections interpretable.

**"MPF_FA slot supported" (census verdict).** A ≥30% share rule calibrated on
thousands of observations was applied to n=11. Minimum-n guard (n≥100) added.

## 3.3 Open items

* **no relaxase anchor** in MPF_T/F/I — mandatory in every CONJScan class
  definition; only MPF_FA records it (82.1%)
* the **98.1% architecture fallback** was measured under the withdrawn VirB5 rule
  and needs re-measurement before reuse
* seven analysis scripts still hardcode the withdrawn rule; left as they ran to
  preserve reproduction, but their VirB5-conditional outputs are undercounts
* seed representativeness unmeasured for MPF_FA
* 1,239 of the 4,747 MPF_FA residual (26.1%) carry no FA/FATA profile
* MPF_FA seed checks are on seed files, not PLSDB (pCF10 absent, ICEBs1
  chromosomal) — validated on seeds, applied to PLSDB

## 3.4 What these catalogues are not

HMM hits and gene positions. **Not operon annotations** — no promoter, terminator
or co-transcription evidence. Slot occupancy is positional adjacency plus, where
available, a named-family hit.
