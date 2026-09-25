# The admission criterion: lipobox alone

## Decision

One criterion decides admission. Everything else is recorded and reported, and
decides nothing.

```
1. POSITION      immediately 3' of VirB5, after an architecture call (98.1%)
2. DECONTAMINATE phmmer vs confirmed members of ALL 14 anchor families
3. ADMIT         family-level lipobox rate >= 50%, reported under all three
                 definitions (strict / relaxed / structural)
4. ANNOTATE      named-Eex hits, Asp@+2, in-slot ratio -- recorded, NOT deciding
```

## Why one and not a conjunction

Every other axis failed a robustness check:

| axis | outcome |
|---|---|
| family count | **anti-correlated** -- VirB10 yields 62 families, target 35 |
| within-cluster presence/absence | dead -- only 8 mixed clusters |
| Asp@+2 | collapsed -- Asp is 5.1% of NF033894's own 525 members, and that family contains pKM101 eex, a confirmed Asp-carrying entry-exclusion protein |
| in-slot / total | null max 0.745 against real families at 0.771/0.812 -- a 0.026 margin |
| **lipobox** | **holds under all three definitions (9-20x), all four seeds score 100% under each, and the 84.2% occupancy correction is invariant to which definition is used** |

Adding a criterion whose null overlaps the signal costs true positives and buys no
precision. On this data that is measurable, not a stylistic preference.

The single-criterion claim is also falsifiable as stated, which a composite score
is not:

> In MPF_T conjugative plasmids, the position immediately 3' of VirB5 is occupied
> by a small lipoprotein at ~20x the rate of size-matched intergenic gaps on the
> same plasmids; of the lipoprotein families at that position, N are recognised by
> no named exclusion family.

## Decontamination is not a second criterion

It is data cleaning, and under a lipobox-only rule it is mandatory rather than
optional. **VirB7/TraN is itself a small lipoprotein, and its detection fails on
97% of plasmids** (PF20898 hits 206/7,436; PF06986 174/7,436). Every unrecognised
VirB7 therefore passes a lipobox test by construction -- an entire class the
criterion cannot see.

The same structural trap has been paid for three times already, because in a
colinear operon **anchor+1 IS the next anchor**:

    candidate pool   unrecognised VirB5    186 unique removed
    recovered set    unrecognised VirB5    190 unique removed
    VirB10 null      unrecognised VirB11   under test

So the check runs against confirmed members of all 14 anchor families, not against
the one expected, and by phmmer rather than by the family HMM -- the HMMs are what
missed these proteins (PF07996: 0/1221 candidates at GA; phmmer: 162/228).

## Two costs, stated

### 1. IncI2 is rejected by the criterion

Its dominant 85 aa occupant (472 of 570 records) has no hydrophobic stretch under
any definition; N-terminus `MKKNDDLITLPSHNSHALDYIIKKAFTPRK...` is charged. It is
**excluded by the criterion, not absent from the data**, and that distinction must
survive into any write-up. It goes to a separate case study, not the main line.

### 2. Scope narrows to MPF_T

The lipobox calibration rests on four MPF_T seeds. Elsewhere in the seed set:

* F `traS` -- UniProt-curated inner membrane, no recorded lipobox
* R64 `excA` -- reported in both membrane-bound and soluble forms
* ICEBs1 `yddJ` -- lipobox at Cys19; SXT `eex` unchecked

So "exclusion genes are lipoproteins" may be an MPF_T property rather than a
general one. The claim is therefore **"small lipoproteins in the VirB5-adjacent
slot of MPF_T plasmids"**, and extension to MPF_F/MPF_I requires recalibration.

## What this defers rather than answers

The skeleton-normalised divergence rate is the only relational axis and has the
thickest literature support (R64 vs R621a: all Tra/Trb >95% identical except ExcA
and TraY; RP4 vs R751: TrbK 37.8% against 75-92% for the rest of Tra2). It is not
used here. If the claim "these genes diverge anomalously fast" is ever made, that
measurement must be done separately -- it is a second paper's claim, not this
one's criterion.

`in-slot/total` is likewise demoted to a descriptive statistic. The pending null
decontamination is still worth reading, since it says whether 0.745 was an
artefact, but it no longer blocks Round 1.
