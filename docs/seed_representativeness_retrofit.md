# Seed representativeness, retrofitted to MPF_T and MPF_F (post-release)

Written after release v1 was sealed. Release files unchanged; this records which
numbers need qualifying and why.

## Why

MPF_I's five seeds all sat in ONE canonical order and ONE layout -- 1 of 118
orders, 45.7% of the class -- and that stratum is the one where ExcA is
essentially universal (99.7%), against 72.7% / 43.6% / 18.5% / 0% elsewhere. The
uniformity that looked like a property of MPF_I was a property of the seed set.
Neither frozen catalogue had this check.

## MPF_F: seeds broader than MPF_I's, but one large stratum unsampled

| seed | layout | order rank |
|---|---|---|
| F, R100 | contiguous | #1 (60.5%) |
| R27 | split_2 | #2 (12.9%) |
| IncC, SXT | split_2 | #3 (5.6%) |

**3 of 298 orders, 2 layouts, covering 79.0%** -- genuinely broader than MPF_I's
1 order / 45.7%.

Named-TraS coverage by stratum:

| stratum | n | slot_ready | named TraS / ready |
|---|---|---|---|
| order#1 | 6,248 | 6,070 | **36.5%** |
| order#4+ | 905 | 335 | 26.6% |
| split_2 | 2,441 | 948 | 22.8% |
| **split_3** | **1,302** | 415 | **2.4%** |

split_3 is 11.8% of the class, carries no seed, and TraS coverage collapses to
2.4%. The headline "32.5% of MPF_F slots hold a named Eex family" is therefore a
mixture dominated by order#1, not a uniform class property.

## MPF_T: the narrowest, and two controls are not in the catalogue

pKM101 (U09868.1, a tra-region fragment) and R46 (AY046276.1) are absent from
PLSDB entirely. The three present -- RP4, R751, R388 -- are **all `canonical`**.

| architecture | n | eex_occupied | candidate | slot_empty |
|---|---|---|---|---|
| canonical | 6,020 | **55%** | 31% | 14% |
| rearranged | 724 | **10%** | 87% | 4% |
| uncallable | 269 | **4%** | 92% | 4% |

Partial mitigation, and it is real: MPF_T's other two architectures were
DISCOVERED from the data rather than seeded -- IncI2 (+2 offset, 121 Mash
clusters) and pEC4115 (-1 offset, 9 independent clusters). So the slot rule there
has been tested outside the seeds' stratum, unlike MPF_F's split_3.

## The sharpest finding: two defining controls fail their own slot-ready test

**RP4 and R751 have `slot_ready__virb5_virb6 = 0`.** Both are `eex_occupied`, so
their slots are demonstrably filled -- but `slot_ready` requires PF07996 to detect
VirB5, and PF07996 scores IncP TrbJ at 17.7 (RP4) and 22.2 (R751) against GA 24.7.

The two plasmids whose TrbK defined the MPF_T slot are excluded by the criterion
built from them.

`T4SS_T_virB5` scores the same proteins at 106 and 122. So release v1's
**slot_ready = 5,272 (70.9%) is a PF07996 artefact in a directly demonstrable
way**, and the demonstration is that the controls are among the excluded.

## Numbers that need qualifying

| number | qualification |
|---|---|
| MPF_T slot_ready 5,272 (70.9%) | PF07996-limited; the defining controls RP4 and R751 are excluded by it |
| MPF_T "46.2% eex_occupied" | canonical-stratum figure; 10% in rearranged, 4% in uncallable |
| MPF_F "32.5% named TraS" | order#1-dominated; 2.4% in split_3, which is 11.8% of the class and unsampled by any seed |
| MPF_T occupancy 82.2% / 84.2% | measured on canonical plasmids, where all three present controls sit |

None of these is wrong. Each is a mixture reported as a single figure, measured on
a set the seeds do not evenly represent.

## Consequence for practice

Run the representativeness check BEFORE freezing, not after. It costs one pass
over an architecture survey that already exists, and it changes how every
stratified figure should be read.
