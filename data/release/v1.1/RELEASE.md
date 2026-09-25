# Release v1.1 (MPF_T only)

Corrects the VirB5 detection rule. v1 remains frozen, sealed and unedited.

| file | rows x cols | what |
|---|---|---|
| `catalogue_MPF_T_v1.1.tsv` | 7,436 x 36 | MPF_T catalogue, corrected VirB5 |
| `catalogue_v1_1_caveats.md` | -- | the v1 -> v1.1 delta and what it corrects |
| `anchor_tiers_v1_1.tsv` | 12 anchors | failure rates, continuous, no binary tier |
| `annotation_MPF_T_F_v1.1.tsv` | 18,381 x 53 | unified annotation, MPF_T block regenerated |
| `catalogue_MPF_I_v1.tsv` | 3,732 x 46 | MPF_I (IncI) catalogue |
| `catalogue_MPF_I_caveats.md` | -- | MPF_I caveats |
| `catalogue_MPF_FA_v1.tsv` | 2,391 x 29 | MPF_FA (Gram-positive) catalogue |
| `catalogue_MPF_FA_caveats.md` | -- | MPF_FA caveats |
| `slot_occupant_families.tsv` | 20,349 x 8 | slot-occupant exclusion family, all 6 families x all 4 classes |

## One change

`VirB5`: `PF07996 @ E<=1e-5 AND aa>150` -> `(PF07996 @ GA AND aa>150) OR
(T4SS_T_virB5 @ GA AND aa>150)`. The withdrawn rule failed the two controls it
was written for, was stricter than the GA it replaced, and was database-size
dependent. See `docs/virb5_evalue_threshold_invalid.md`.

## Headline

| | v1 | v1.1 |
|---|---|---|
| VirB5 present | 5,424 (27.1% fail) | 6,881 (7.5% fail) |
| `slot_ready__virb5_virb6` | 5,272 (70.9%) | 6,548 (88.1%) |
| `tier`, `architecture` | -- | unchanged (asserted, build aborts on drift) |

## Superseded figures

`VirB5 27.1% -> 3.1%` and `slot_ready 70.9% -> 87.2%` were reported in session and
are **wrong** -- numerator over 72,556 plasmids, denominator over 7,436. Correct
values are in the table above.

## Unified annotation regenerated

`annotation_MPF_T_F_v1.1.tsv` (18,381 x 53). Two columns added:
`mpft_VirB5_arm` (which arm produced the bitscore) and `mpft_virB5_source`.

VirB5 bitscores are RESCANNED from both arms rather than read from
`anchor_matrix.tsv`, which holds PF07996 scores only -- the 1,336 CONJScan-only
plasmids would otherwise carry an empty score in the one column whose definition
changed.

Verified against v1: rows 18,381 unchanged, accessions identical, `dual_system`
175 unchanged, **MPF_F block identical on all 18,381 rows** (0 differing), VirB5
presence flipped on 1,457, `mpft_slot_ready__virb5_virb6` 5,272 -> 6,548.

The build ABORTS if the VirB5 rescan disagrees with `catalogue_MPF_T_v1.1` on any
plasmid. It agreed on all 7,436 -- and since the catalogue was built on the 1.1M
cache and the merge scanned the 7.7M one, that agreement is also an empirical
confirmation that GA is database-independent, which is the property the withdrawn
E-value rule lacked.

## NOT reissued

`catalogue_MPF_F_v1.tsv` stays at v1: its entry criterion (PF11130) is GA-based
and the MPF_F block above is byte-identical to v1.

## MPF_I added

Entry `T4SS_I_traU` (5/5 seeds, 0.46% cross-class). Slot TraY+1,
`slot_ready__traY` 3,607 / 3,732 = 96.7%, same three-part definition as
`slot_ready__virb5_virb6`.

`slot_admission` is `undetermined` on every row; `slot_maxkd` is carried as a
column (median 2.26) and never used as a filter -- hydrophobicity separates at
only 2.3x against a 40.5% background.

Two ExcA figures are recorded together with their denominators, because reading
one as the other is the error that produced the withdrawn `VirB5 3.1%`:
**59.3%** (2,213/3,732) of admitted plasmids carry ExcA; **96.5%** (2,213/2,293)
of PLSDB's ExcA carriers are recovered by the entry criterion.

`excB` is reported in seven columns rather than one call, and unscanned rows carry
blank rather than 0 -- `NC_011077.1` is absent from the nucleotide FASTA and is
labelled `genome_not_retrieved`, not `no candidate`.

## MPF_FA added -- four classes complete

Entry: `T4SS_virb4 AND NOT (T or F or I) AND >=2 FA/FATA profiles` -> 2,391.
Cross-class is excluded BY CONSTRUCTION, not measured; no 0% figure is quoted.

**MPF_FA has NO SLOT COLUMN.** `conG+3` is confirmed on ICEBs1 alone: 1 of 4
locatable seeds (two occupants anti-oriented, one an identified toxin-antitoxin
system), and DUF4467 occurs on just 11 of 2,391 admitted plasmids (0.5%) against
2,725 / 2,213 / 224 for the other three classes. `slot_status = no_slot_defined`
on every row. A one-element slot column would be empty, not caveated.

`tier` is replaced by `profile_support` (minimal/medium/high) -- the specified
"entry hit + slot-flanking anchors" is impossible with no slot, and A/B/C/D would
import MPF_T's 7/7-core meaning, which does not exist for a class whose every
class-specific profile is `accessory`.

## Status against the minimum bar

| class | entry | slot verified | named-family control | caveats |
|---|---|---|---|---|
| MPF_T | PF03135 | VirB5+-1 | TrbK 224 | yes |
| MPF_F | PF11130 | TraG_N+1 | TraS 2,725 | yes |
| MPF_I | T4SS_I_traU | TraY+1 | ExcA 2,213 | yes |
| MPF_FA | virb4 conjunction | **NONE** | DUF4467 **11** | yes |

MPF_FA meets the entry and caveats rows. It does NOT meet the slot row, and that
is recorded rather than worked around.

# FROZEN -- final counts

## Gene clusters and named entry-exclusion genes

| class | catalogued | slot-resolvable | named eex gene | unnamed candidate |
|---|---|---|---|---|
| MPF_T | 7,436 | 6,548 | **3,862** (59.0%) | 2,655 |
| MPF_F | 11,120 | 7,803 | **2,541** (32.6%) | 5,262 |
| MPF_I | 3,732 | 3,607 | **2,112** (58.6%) | 1,495 |
| MPF_FA | 2,391 | none (no slot) | **11** (plasmid-level) | -- |
| **TOTAL** | **24,679 rows / 24,470 distinct** | 17,958 | **8,526** | **9,412** |

## By family

| family | model | n | class |
|---|---|---|---|
| Eex_IncN | NF033894 | 3,644 | MPF_T |
| TraS | PF10624 | 2,541 | MPF_F |
| ExcA | NF033891 | 2,112 | MPF_I |
| TrbK_RP4 | TIGR04359 | 218 | MPF_T |
| DUF4467 | PF14729 | 11 | MPF_FA |
| EexR | NF041429 | **0** | -- |

## What these numbers are NOT

* **NOT operons.** No promoter, terminator or co-transcription evidence anywhere
  in the pipeline. These are anchor sets in measured adjacency.
* **8,526 is a LOWER BOUND.** It counts members of six already-named families.
  The 9,412 unnamed slot occupants are a CANDIDATE CONTAINER, not absence.
* **Only MPF_T has a usable criterion** for sifting candidates (lipobox: 5
  families at >90%, zero at three negative slots, positive control recovered).
  MPF_F and MPF_I have none -- every compositional axis tested failed or is
  unavailable. MPF_FA has no validated slot at all.

## Cross-class check

All six families were scanned in all four classes. **No family appears at another
class's slot.** The v1 census found TraS on 28 MPF_T plasmids and Eex_IncN on 150
MPF_F plasmids -- those are elsewhere on the replicon, not at the slot.
`NF041429` (EexR) hits ZERO proteins across all 20,269 annotated plasmids.

## MPF_T figures superseded

Earlier logs recorded eex_occupied 46.2% / candidate 40.8% / slot_empty 13.0% on
v1's 5,272-plasmid slot_ready set. This release uses v1.1's 6,548 with the
corrected VirB5 rule: **59.0% / 40.5% / 0.5% unresolved**. `slot_empty` does not
appear here because this annotation pass classifies every resolvable +1 gene
rather than dropping anchors.

## Criterion search: CLOSED

| axis | MPF_T | MPF_F | MPF_I |
|---|---|---|---|
| lipobox | **works** (56.1% vs 2.8%) | INVERTED (4.1% vs 6.1%) | negligible (4.4%) |
| length | -- | no signal (P=0.458/0.538) | fails (204 vs 210) |
| KD | -- | untested | 2.3x on 40.5% background |
| TM count | unavailable (fails controls, 3 attempts) | same | same |
| family count | ANTI-CORRELATED, dead | -- | -- |
| in-slot/total | rejects the known family (0.095) | -- | -- |

Untried: pI / amino-acid composition (expected to behave like the other
compositional axes), and transcriptional independence (own promoter,
Rho-independent terminator) -- the only axis orthogonal to sequence and position.
