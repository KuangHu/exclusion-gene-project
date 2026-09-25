# CONJScan reopens the anchor sets (post-release finding)

Written AFTER release v1 was sealed. The release files are unchanged; this records
what a v2 would have to revisit.

## What happened

MPF_I's seed readout scanned Pfam-A only and reported three components with NO
family at GA: TraU (DotO/IcmB, the VirB4 homolog), TrbC (DotL, the T4CP) and TraY
(DotA). The conclusion drawn was "MPF_I breaks the VirB4-ATPase convergence".

**That was a library artefact.** CONJScan has been on disk the whole time
(`funcannot_dbs/macsy_models/CONJScan`, 404 profiles, 17 of them `T4SS_I_*`) and
all three components are detected strongly and uniformly:

| CONJScan profile | R64 | R621a | pEK204 | pCVM29188 | ColIb-P9 |
|---|---|---|---|---|---|
| `T4SS_I_traU` (GA 413.1) | **1460** | 1459 | 1459 | 1455 | **1461** |
| `T4SS_I_traY` | 951 | 879 | 949 | 952 | 952 |
| `T4SS_t4cp2` | 105 | 104 | 104 | 105 | 105 |
| `T4SS_virb4` (generic) | — | — | — | — | — |

CONJScan's own `T4SS_typeI` definition says `T4SS_virb4` is **mandatory**, with
`T4SS_I_traU` as its exchangeable -- so the convergence holds across all four
classes: MPF_T VirB4, MPF_F TraC, MPF_FA ConE, MPF_I TraU.

Note the generic `T4SS_virb4` scores ZERO on TraU while the class-specific
`T4SS_I_traU` scores 1460. Same lesson as PF19044 (57% purity) vs PF11130 (93%)
being two domains of one protein: **the right protein with the wrong model reads
as absence.** Fifth instance of a structural-looking conclusion turning out to be
a model artefact, and the first where the fix was already on disk.

## Why this reaches back into the releases

CONJScan was used from MPF_T step one, but only for CLASS ASSIGNMENT. Its
component profiles were never used as anchors. **Both frozen catalogues are built
entirely on Pfam anchors.**

So the two-tier failure structure may be partly a library property:

| anchor | GA failure | Pfam model |
|---|---|---|
| VirB5 | 27.1% | PF07996 -- scores IncP TrbJ at 17.7 / 22.2 against GA 24.7 |
| VirB2 | 20.9% | PF04956 |
| VirB1 | 19.5% | PF01464 (SLT -- 22,605 admissions at 25.6% purity, promiscuous) |

VirB5's 27.1% is the root of a large amount of scaffolding: the
`E<=1e-5 AND aa>150` conjunct, the architecture fallback, every
"conditional on VirB5 detectable" label, and the 5,101 vs 7,014 denominator gap.
**If `T4SS_T_virB5` works, most of that is unnecessary.**

## The check to run first

Ten seeds (five MPF_T, five MPF_F) against `T4SS_T_*` and `T4SS_F_*` profiles,
scored beside the Pfam equivalents. Cheap, seed-level, no database scan.

Risk direction, and why it is lower than the phmmer recall that had to be
retracted: that failure was "more sensitive -> false positives inflate
completeness", and its precision control was structurally blind on short
peripheral anchors. CONJScan profiles are curated, class-specific models with
their own GA thresholds, so the failure mode is different -- but purity must still
be measured the way MPF_F's was, since that is what disqualified four cross-class
families there.

## Decision deferred

Run the ten-seed comparison first, then choose: continue MPF_I's entry criterion,
or revisit the two catalogues on CONJScan anchors. Release v1 stands either way --
it is internally consistent and its limits are documented.

## Second gap: no relaxase anchor at all

`T4SS_MOBB` and its nine exchangeables (MOBV, MOBT, MOBQ, MOBP1-3, MOBH, MOBF,
MOBC) are `presence="mandatory"` in EVERY CONJScan class definition. **Neither
frozen catalogue scores a relaxase.**

So tier A as defined here is not aligned with CONJScan's system definition: a
plasmid can be tier A by our completeness measure while lacking the relaxase
CONJScan requires. For the stated goal -- a trustworthy tra cluster record -- that
matters directly, because a conjugation system without a relaxase cannot mobilise
DNA and is not a complete system.

Measure it in the same pass (same scan, no extra cost):
* which MOB type each of the ten seeds carries (expect F -> MOBF/TraI,
  RP4 -> MOBP/TraI, R64 -> MOBP/NikB)
* how many tier A plasmids database-wide lack any relaxase

The second number decides whether the tier definition changes: ~95% carrying one
means add a column; ~20% missing means tier needs revising.

## The accessory labelling supports the two-tier result rather than contradicting it

All nine `T4SS_T_*` are `presence="accessory"`; only virb4 + t4cp1 + relaxase are
mandatory. That is a different axis from the measured two-tier structure:

| | what it measures |
|---|---|
| CONJScan mandatory/accessory | whether the system can function without it |
| our core 4.5-4.9% vs peripheral 8-27% | whether the component can be reliably detected |

VirB6/B8/B9/B10/B11 are **accessory in CONJScan yet fail on only ~4.8%** -- so they
are functionally replaceable but sequence-conserved. VirB1/B2/B3/B5 are accessory
AND poorly detected. The two axes are orthogonal and the combination is worth
keeping.

## Full specification for the next run

    ten seeds (five MPF_T + five MPF_F)
      x  9 T4SS_T_*  +  12 T4SS_F_*  +  10 T4SS_MOB*
      x  the corresponding Pfam families
    -> bitscore and GA side by side

then purity, measured as MPF_F's was, on any anchor where the two libraries
disagree materially -- VirB5 first.
