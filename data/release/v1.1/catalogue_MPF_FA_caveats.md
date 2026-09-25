# MPF_FA (Gram-positive) catalogue -- caveats

`catalogue_MPF_FA_v1.tsv`, 2,391 x 29. Jobs 25974613, 25991904, 26020027,
26035629, 26038132, 26039112.

## Entry criterion

```
T4SS_virb4  AND NOT (MPF_T or MPF_F or MPF_I)  AND  >=2 FA/FATA profiles
```

**Cross-class contamination is EXCLUDED BY CONSTRUCTION, not measured.** The
`NOT` clause removes those plasmids by definition. No "0% cross-class" figure is
quoted, because that would read as "measured and found clean".

`T4SS_virb4` alone was REJECTED: 23,191 admissions at **79.5%** contamination
against a rule of <=30% fixed before the run (chance baseline 30.5%). It admits
99.99% of MPF_T and 99.6% of MPF_F but only **3.9% of MPF_I** -- the generic
profile fails on an entire class while saturating two others.

The `>=2` threshold was chosen AFTER seeing the background table (`>=1` gives only
2.1x over MPF_F; `>=2` gives 84x over MPF_F and 19.4x over MPF_T). It is post-hoc,
with one independent check: all five seeds carry >=2, and the seeds played no part
in choosing it.

## THERE IS NO SLOT COLUMN. THIS IS THE MAIN LIMITATION.

`slot_status = no_slot_defined` on all 2,391 rows.

The `conG+3` slot is confirmed on ICEBs1 and nowhere else.

**Seed check (job 26035629) -- 1 of 4 locatable seeds:**

| seed | slot occupant | strand vs anchor | family |
|---|---|---|---|
| ICEBs1 | `conJ` 126 aa | same | **PF14729 99.1** |
| pAD1 | 108 aa type II toxin-antitoxin (antitoxin MazE at +1) | same | none |
| pLS20 | 45 aa hypothetical | **opposite** | none |
| pAM373 | 39 aa hypothetical | **opposite** | none |
| pCF10 | anchor absent | -- | -- |

Strand orientation was verified against the GenBank records: there is no
direction bug, the offsets are taken correctly, and the result is genuinely
negative. Two of four occupants are ANTI-ORIENTED to the anchor and so cannot
plausibly be co-transcribed with it -- the offset is landing on whatever is
adjacent. pAD1's is an identified non-candidate, not an unknown.

**Class census (job 26038132) -- DUF4467 is effectively absent:**

| class | named family at the slot | recoverable across the class |
|---|---|---|
| MPF_F | TraS / `PF10624` | 2,725 |
| MPF_I | ExcA / `NF033891` | 2,213 |
| MPF_T | TrbK / `NF033894` | 224 |
| **MPF_FA** | DUF4467 / `PF14729` | **11 of 2,391 = 0.5%** |

Within those 11, the modal offset from `FA_orf15` is **-1**, not ICEBs1's **+3**.
The one case we know is not the modal case.

A slot column validated on a single element would not be a caveated column, it
would be an empty one. It is omitted rather than shipped with a warning.

### A wrong verdict was printed and is withdrawn

The census script reported "dominant offset -1 at 45.5% -> CONSISTENT: slot
supported". That applied a >=30% share rule calibrated on MPF_T/MPF_I (n in the
thousands) to **n = 11**, where 4-of-11 is also chance. The rule fired where it
had no power to discriminate -- the same defect as the vacuous VirB3 length band.
A minimum-n guard (n >= 100) was added to the script. **The slot is NOT
supported.**

## `tier` was replaced by `profile_support`

The specified tier was "entry hit + slot-flanking anchors present". With no slot
the second half is impossible.

`profile_support` counts FA/FATA profiles: `minimal` 757 (2), `medium` 593 (3-4),
`high` 1,041 (>=5). Median 4, max 16.

It is deliberately NOT called A/B/C/D. Those letters mean 7/7 structural core in
MPF_T, and every CONJScan Gram-positive class-specific profile is `accessory` --
there is no mandatory core to count, so importing the letters would import a
meaning this class has no basis for.

## `fa_subclass`

FATA 1,170 (48.9%), unassigned 726 (30.4%), FA 495 (20.7%). Same 2x-margin rule as
the seeds; `unassigned` where neither set reaches a 2x margin.

## CAVEAT: the seed check is on SEED FILES, not PLSDB

5/5 was measured on the GenBank seed records. **pCF10 is absent from PLSDB and
ICEBs1 is a chromosomal ICE.** A PLSDB seed check would read 0/5 for reasons
unrelated to the criterion -- the A12 failure mode already seen when MPF_F
reported 0/5 because PLSDB is RefSeq and the seeds were INSDC.

The criterion is validated on seeds and APPLIED to PLSDB. Those are not the same
population.

## Recorded, not explained

1,239 of the 4,747 conjunction residual (26.1%) carry NO FA/FATA profile and are
excluded by the `>=2` requirement.

## Relaxase

1,963 of 2,391 (82.1%) carry a `T4SS_MOB*` profile -- recorded in
`relaxase_profiles`. A relaxase is the third mandatory gene in every CONJScan
class definition and is still absent as an anchor from the other three
catalogues.

## Not done, by scope

architecture survey; seed representativeness; admission criterion; candidate-pool
analysis.

## NOT an operon annotation

HMM hits and gene positions. No promoter, terminator or co-transcription evidence.
