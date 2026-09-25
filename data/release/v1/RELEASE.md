# Release v1 -- MPF_T and MPF_F tra gene cluster annotation over PLSDB

Frozen. Every file here is sealed by SHA256 in `config/curated_files.tsv` (A9) and
verified by `scripts/74_validate.py`. Generators write `*_regen.tsv` in `data/` and
never this directory; to change a release, regenerate, diff, and promote by hand.

NOT an operon annotation. HMM hits and gene positions only -- no promoter,
terminator or co-transcription evidence was gathered, and calling it an operon
would foreclose the moron line (RP4's trbJ-trbK terminator, IncC sfx's own
constitutive promoter).

## Files

| file | rows x cols | entry criterion | slot |
|---|---|---|---|
| `catalogue_MPF_T_v1.tsv` | 7,436 x 34 | PF03135 VirB4 | VirB5 +/-1 |
| `catalogue_MPF_F_v1.tsv` | 11,120 x 30 | PF11130 TraC_F_IV | TraG_N +1 |
| `annotation_MPF_T_F_v1.tsv` | **18,381 x 51** | both, as parallel blocks | both |
| `catalogue_v1_caveats.md` | | MPF_T caveats | |
| `mpf_f_entry.md` | | MPF_F entry derivation | |

18,381 = 7,436 + 11,120 - 175 admitted by both.

## Headline numbers

| | MPF_T | MPF_F |
|---|---|---|
| tier A | 6,717 (90.3%) | 9,994 (89.9%) |
| slot-ready | 5,272 (70.9%) | 7,803 (70.2%) |
| slot occupied by a named Eex family | 46.2% | 32.5% (all TraS) |
| occupant length | ~70-85 aa | median 173 |
| occupant lipobox-positive | 56.1% | **1.8%** |

Two independently built pipelines -- different entry criteria, anchor sets and
slot definitions -- give nearly identical tier A and slot-ready proportions. That
was not designed; it is two measurements agreeing.

## What does NOT transfer between classes

* **The lipobox criterion.** 56.1% in MPF_T against 1.8% in MPF_F. F traS is
  UniProt-curated inner membrane with no lipobox, and the database says that is
  the class norm. `docs/admission_criterion.md` scoped the claim to MPF_T in
  advance; this measures it.
* **Which side of VirB6 the slot sits on.** MPF_T's slot is 5' of VirB6, MPF_F's
  is 3' of the VirB6 analogue. Opposite sides. MPF_F also has no VirB5 analogue in
  its anchor set. Recorded, not explained.
* **Slot adjacency itself.** MPF_FA's exclusion gene (ICEBs1 yddJ/conJ, DUF4467)
  sits at conG+3, not adjacent.

## Known limits

* **VirB7 is unscoreable.** GA detects it on 206/7,436 (2.8%); phmmer recall
  recovered zero; six-frame in VirB6|VirB8 gaps returns ORFs at 2.3% lipobox,
  indistinguishable from background, while pKM101 traN itself has one.
* **Second-pass phmmer recall is an ANNOTATION, not a measure.** Its precision
  control is blind where it claims most. See `docs/ga_threshold_failure.md`.
* **`n_unresolved_gaps` needs its 10 kb ceiling** (median 3, Q3 8 with it; median
  11, Q3 36, max 4,400 without). Anchor-list adjacency is not genome adjacency.
* **MOB-suite `mpf_type` is a reference, not ground truth** -- MPF_T used
  CONJscan. Its single-valued field cannot represent the 144 plasmids that are
  tier A in BOTH blocks.

## Named exclusion families we currently cannot see

Census over all 72,437 cached accessions (GA):

| family | plasmids | in MPF_F | in MPF_T |
|---|---|---|---|
| NF033894 Eex_IncN | 4,740 | 150 | 4,559 |
| TraS PF10624 | 3,011 | 2,725 | 28 |
| **NF033891 ExcA** | **2,293** | **18** | **10** |
| TIGR04359 TrbK | 263 | 6 | 224 |
| DUF4467 YddJ | 13 | 0 | 0 |
| NF041429 EexR | 4 | 4 | 0 |

**~2,278 plasmids carry a named exclusion family and fall outside both entry
criteria** (ExcA 2,265 + YddJ 13). These are IncI1/MPF_I and Gram-positive
MPF_FA. That is a cost of the entry criteria, not a failure of the slot rules,
and it sizes what MPF_I and MPF_FA would add: a named-positive population about a
third of MPF_T's.

The candidate pool is only lightly diluted: at most 191 of the 5,269 unnamed MPF_F
slot occupants (3.6%) carry a TraS that is not at TraG_N+1.

## Open, recorded, not worked

* MPF_F's admission criterion needs recalibration on TraS-class proteins.
* MPF_I (IPR027628 must be resolved to member signatures; excA/excB share a stop
  codon so six-frame is mandatory) and MPF_FA (ConG has only the generic MetI-like
  fold -- uninformative, NOT evidence against Avello 2019).
* The 84 aa IncI2 slot lipoprotein: 25 independent Mash clusters, no named family,
  passes the MPF_T criterion. First candidate through.
