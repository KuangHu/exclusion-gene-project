# MPF_I (IncI) catalogue -- caveats

`catalogue_MPF_I_v1.tsv`, 3,732 x 46. Job 25972622.

## Entry criterion

`T4SS_I_traU` at GA. 3,732 plasmids, 5/5 seeds, 0.46% cross-class.
Chosen over `T4SS_I_traY` (3,922 admissions but 86.2% purity) and the higher-purity
but lower-yield `traT`/`traV`/`traE`. MOB-suite purity was NOT used as a gate:
MPF_I already showed it to be class-dependent.

## Two ExcA numbers, different denominators -- do not conflate

| statement | value | denominator |
|---|---|---|
| admitted plasmids carrying ExcA | 2,213 / 3,732 = **59.3%** | the MPF_I catalogue |
| PLSDB ExcA carriers recovered by the entry criterion | 2,213 / 2,293 = **96.5%** | all ExcA-positive plasmids |

The second is the recall figure previously reported as 96.6%. Neither implies the
other. They are recorded together because reading one as the other is exactly the
error that produced the withdrawn `VirB5 3.1%` figure.

## Slot

TraY+1. `slot_resolved` 3,625 (97.1%); `slot_ready__traY` 3,607 (96.7%).

`slot_ready__traY` uses the v1.1 three-part definition -- anchor present AND slot
occupant resolvable AND architecture callable -- identical in shape to
`slot_ready__virb5_virb6`, so the two mean the same thing when the catalogues meet
in a unified annotation.

## `slot_admission` is `undetermined` on every row -- by decision

`slot_maxkd` is a COLUMN, never a filter. Measured: median 2.26, Q1 1.97, Q3 2.54.

Hydrophobicity separates candidates from background by only 2.3x against a 40.5%
background rate. MPF_T's lipobox criterion (56.1% in-slot vs 2.8% background, 20x)
does NOT transfer -- MPF_F measures 1.8%. Carrying the number without a verdict is
the honest form. MPF_F and MPF_I are both `undetermined`; MPF_FA will follow.

## excB: four columns, never collapsed to one call

excB is an in-frame reinitiation ORF inside excA and shares excA's stop.

| column | meaning |
|---|---|
| `excA_hit` | NF033891 at GA |
| `excB_scanned` | was the genome actually retrieved and searched |
| `excB_inframe_orfs` | ALL in-frame candidates >=30 aa, not only 147 aa |
| `excB_n_candidates` | count; **blank, not 0, when unscanned** |
| `excB_has_147aa` | seed-matching length; **blank, not 0, when unscanned** |
| `excB_family_hit` | any exclusion family on the plasmid |
| `excB_call_mode` | `sixframe_inframe` 2,212 / `no_excA` 1,519 / `genome_not_retrieved` 1 |

`excB_has_147aa` = **1,067 of 2,212 SCANNED** -- not of 2,213, not of 3,732.

Reporting only the 147 aa seed length would be fitting the answer to the seeds, so
every candidate length is listed and `has_147aa` is a separate flag.

### `genome_not_retrieved` is not `no candidate`

`NC_011077.1` carries ExcA but is absent from the PLSDB nucleotide FASTA, so its
excB was never searched. An earlier build labelled it `excA_but_no_candidate` --
reporting "we did not look" as a negative result. It now has its own mode and its
count columns are blank rather than 0.

## Architecture

`layout` and `order_canonical` are carried from the frozen survey. Order
canonicalisation was NOT redone (out of scope). 58 plasmids are `uncallable`;
strata: order#1 1,544, split_2 1,195, order#3 439, order#4+ 305, order#2 132,
split_3 50, split_4 6, split_6 2, split_5 1.

`order#2` reproduces at 132, matching the prior analysis -- the stratum where
131/132 carry a 217 aa protein at TraY+1 across 22 Mash clusters, 9 unique
sequences, 99% relaxase-positive.

## Red lines observed

Exclusion families (`NF033891` ExcA included) are scanned AFTER the catalogue is
built and feed ONLY the post-hoc `excA_*` / `excB_*` block. They never enter the
entry criterion, the anchors, the architecture or the slot definition.
`assign_method = positional` never feeds `slot_status`.

## Not done, by scope

* seed representativeness -- open item
* architecture order canonicalisation -- layout only
* candidate-pool analysis -- the container is recorded, not analysed
* KD at unique-sequence level (175 unique vs 1,205 records) -- open
* no relaxase anchor (mandatory in every CONJScan class definition) -- open

## NOT an operon annotation

HMM hits and gene positions. No promoter, terminator or co-transcription evidence.
