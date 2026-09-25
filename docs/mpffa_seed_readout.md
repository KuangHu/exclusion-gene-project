# MPF_FA §1.1 -- five Gram-positive seeds against every FA and FATA profile

Job 25974613. Five seeds, 323 proteins, 7 FA + 27 FATA + 14 mandatory + 6
exclusion profiles.

## Positive control passed BEFORE anything else was read

pCF10 is the FATA reference strain (`FATA_prgB/C/F/H/I/K/L` are derived from it),
so pCF10 must fire on FATA or the scan is broken. **8 of 27 FATA profiles fire on
pCF10.** The script raises before printing any other result if this fails.

## Class assignment -- MEASURED, not inferred from organism

| seed | FA | FATA | assignment |
|---|---|---|---|
| **ICEBs1** | **5** | 1 | **FA** |
| pCF10 | 0 | 8 | FATA |
| pAM373 | 3 | 2 | AMBIGUOUS |
| pLS20 | 1 | 1 | UNASSIGNED |
| pAD1 | 1 | 1 | UNASSIGNED |

Rule: `UNASSIGNED` if fewer than 2 hits on the stronger side; the class wins only
with a 2x margin; otherwise `AMBIGUOUS`. An earlier version called 1-vs-1 "both",
which overstates one hit per side as dual membership -- corrected.

**ICEBs1 is FA.** That is the question this step existed to answer, and it is a
measurement. Note pAM373 is *E. faecalis* like pCF10 and does NOT come out
cleanly FATA, so class does not follow organism -- which is why inferring ICEBs1
from *B. subtilis* would not have been safe.

`UNASSIGNED` is a real outcome, not a failure: every class-specific profile in
both sets is `accessory`, so a genuine MPF_FA element may hit almost none.

## The VirB4 convergence holds on a fifth class -- 5/5

| seed | `T4SS_virb4` |
|---|---|
| pCF10 | 279.8 |
| pLS20 | 257.2 |
| pAD1 | 147.6 |
| ICEBs1 | 111.7 |
| pAM373 | 91.0 |

`T4SS_virb4` is the ONLY profile firing on all five, and it is mandatory in every
CONJScan plasmid class definition. On ICEBs1 it corresponds to `conE`, annotated
"VirB4-like ATPase ConE", 831 aa. This is the entry-criterion candidate for §2.

Coupling protein: `T4SS_t4cp2` fires on 4 of 5 (not ICEBs1, which instead carries
`T4SS_tcpA` 544.6 and `T4SS_MOBT` 443.9). Relaxase coverage is split across
`MOBC/P1/P2/T` -- no single relaxase profile covers the class.

## Slot positive control ESTABLISHED

```
PF14729 (DUF4467)   ICEBs1 conJ   99.1 bits (GA 25)   126 aa   16495..16875
```

99.1 against GA 25 is decisive. The position is exactly `conG+3`:

| offset | gene | coords | aa |
|---|---|---|---|
| 0 | `conG` | 12478..14925 | 815 |
| +1 | `cwlT` | 14922..15911 | 329 |
| +2 | `ACGWZ1_02860` (= `yddI`) | 15926..16432 | 168 |
| **+3** | **`conJ`** | **16495..16875** | **126** |

### This control was NOT previously established

`PF14729.hmm` and `PF10624.hmm` **did not exist** in the HMM directory. The first
run of this script skipped missing profiles silently and reported nothing for the
control -- which read as "the control did not fire" when the truth was "the model
was never on disk". Both were fetched from Pfam-A (`hmmfetch`).

The script now ABORTS if any named exclusion model is missing. A missing model
read as a missing hit is the same error class as reading "we did not look" as a
negative result.

Before this, the only evidence that `conJ` carries DUF4467 was the RefSeq product
string "DUF4467 domain-containing protein" -- RefSeq's own annotation pipeline,
not our measurement. Structurally the same circularity as the R388 "EexN family
lipoprotein" case, though milder since DUF4467 is not a mechanism name. It is now
a measurement.

## THE SLOT CONTROL COVERS ONE ELEMENT OF FIVE

This is the structural weakness of MPF_FA and must not be smoothed over by the
fact that the control passed.

| class | named family at the slot | recoverable across the class |
|---|---|---|
| MPF_T | TrbK / `NF033894` | **224** plasmids |
| MPF_F | TraS / `PF10624` | **2,725** plasmids |
| MPF_I | ExcA / `NF033891` | **2,213** plasmids |
| **MPF_FA** | **DUF4467 / `PF14729`** | **1 element (ICEBs1 only)** |

MPF_T, MPF_F and MPF_I each have a named exclusion family recoverable at
class scale, so their slot definitions are validated against thousands of
independent instances. MPF_FA's control is a single seed.

No exclusion family fires on pLS20, pCF10, pAD1 or pAM373 -- not `PF14729`, not
`PF10624`, not any of the four NCBI families. `PF14729` on ICEBs1 `conJ` is the
whole of the evidence.

**Consequence:** the MPF_FA slot validation base is weaker than the other three by
roughly three orders of magnitude. `conG+3` is confirmed on ONE element. Whether
the offset generalises across the class is UNMEASURED, and a class-scale
`PF14729` census is required before any MPF_FA slot figure is quoted with the
confidence the other three carry.

## Negative results, stated

* `PF10624` (TraS) fires on NONE of the five seeds.
* No exclusion family fires on pLS20, pCF10, pAD1 or pAM373. The named-family
  positive control exists for **ICEBs1 only** -- one seed of five.
* 0 of 7 FA profiles fire on pCF10; 0 FA and 0 FATA class-specific profiles reach
  2 hits on pLS20 or pAD1.

## Counts corrected

FATA has **27** profiles (not 26 -- `FATA_cd419_1` sorts separately from the
`FATA_cd419a/b` pair). pAD1 contributes **81** CDS with translations, not 82; one
CDS carries no `/translation`.

## Open items (deliberately not done, per scope)

* seed representativeness -- not measured
* architecture / order canonicalisation -- layout only
* whether `pLS20`/`pAD1` belong to a Gram-positive class at all, or are simply
  not covered by either profile set
