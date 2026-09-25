# The VirB5 `E<=1e-5` threshold was justified by E-values from a different search

## Status: the threshold is WITHDRAWN. Replaced by `PF07996 @ GA OR T4SS_T_virB5 @ GA`, both `aa>150`.

## The claim in `anchor_set.tsv`

The VirB5 row is the only non-GA row in the anchor set. Its stated rationale:

> GA IS INSUFFICIENT. At GA (24.7) only pKM101 (87.4) and R388 (122.9) pass;
> RP4 trbJ 17.7 and R751 22.2 both FAIL. By E-value all four pass
> (7.9e-07, 3.4e-08, 3.5e-28, 4.6e-39).

## What was measured (job 25913366, hmmsearch PF07996 -> the 1.1M cache)

| control | note bits | measured bits | note E | measured E | `E<=1e-5` |
|---|---|---|---|---|---|
| RP4 trbJ (BN000925.1) | 17.7 | **17.7** | 7.9e-07 | **0.22** | **FAIL** |
| R751 trbJ (NC_001735.4) | 22.2 | **22.2** | 3.4e-08 | **0.0093** | **FAIL** |

The bitscores reproduce exactly. The E-values differ by five to six orders of
magnitude at identical bitscores, which can only mean the two numbers came from
different searches. `E = P * N`: scanning one protein against Pfam-A gives
N ~ 30,134 models; searching one model against the proteome gives N = 1,099,911
sequences. The note's E-values have the shape of the former; the pipeline runs
the latter.

**Both controls the rule was introduced to rescue FAIL that rule under the search
the pipeline actually performs.**

## The threshold is also backwards

Measured boundary (job 25913252): in a 1,099,911-protein hmmsearch, `E = 1e-5`
falls at **~31.7 bits** -- well ABOVE the GA of 24.7.

| rule | plasmids (of 7,436) |
|---|---|
| PF07996 @ GA 24.7 | **5,545** |
| PF07996 @ `E<=1e-5` (~31.7 bits) | **5,368** |

Switching from GA to `E<=1e-5` **tightened** VirB5 by ~177 plasmids while being
documented as a loosening to recover divergent members. Every VirB5-conditional
figure in release v1 rests on the tighter rule.

## `E`-value thresholds are database-size dependent and must not be used

The same rule on the same sequences gives 5,368 in a 1.1M-protein search and
**5,341** in a 7.7M-protein search. Nothing about the proteins changed; only N did.
Bitscore thresholds (GA) are invariant -- confirmed here, since VirB1 (5,984),
VirB2 (5,884) and VirB3 (6,769) reproduced to the digit across both caches, and
the two caches were verified byte-identical on all 7,436 accessions
(0 disagreements in per-accession protein counts).

**Rule adopted: no anchor threshold may be an E-value.** GA, or an explicit
bitscore derived from measurement.

## Unreproduced residual

`anchor_set.tsv` records 5,424 for VirB5 (27.1% failure). The stated rule
yields **5,368** on the very cache v1 was built from. 56 plasmids are
unaccounted for by either the db-size effect or the rule as written. The
provenance of 5,424 is unknown; it is NOT reproducible from the recorded
threshold. Recorded, not explained.

## What replaces it

`T4SS_T_virB5` scores the two failing controls at **105.6** (RP4) and **122.3**
(R751) -- not marginal. The disjunction passes both purity controls that the
retracted phmmer recall failed:

* new admissions median 258 aa against the family's own 238 (the retracted recall
  gave 816 against 238)
* 91.4% inside the measured robust band
* 6.2% cross-class MPF_F

## Corrections to figures previously reported

| figure | as reported | corrected |
|---|---|---|
| VirB5 failure, disjunction | 27.1% -> 3.1% | **25.4% -> 7.5%** (GA both sides). With the withdrawn E-value rule it was 28.2% -> 8.9% |
| `slot_ready` | 70.9% -> 87.2% | **withdrawn**. Job 25913461 gives a CEILING of 91.5% (6,806/7,436 carry VirB5-by-disjunction AND VirB6). That is not `slot_ready`: it omits the callable-architecture requirement, so the true figure is lower. No replacement number is quoted until architecture is applied |
| "VirB5 enters the core tier" | asserted | **neither tier.** On GA both sides (E-value dropped) VirB5 is 25.4% -> **7.5%**, against a core band of 4.1-4.9% and a peripheral band of 13.2-20.9%. It is ~1.5x the worst core anchor and ~1.8x below the best peripheral one. A binary tier label is not supported at this value |

The `3.1%` arose from dividing a union counted over all 72,556 plasmids by the
7,436-plasmid MPF_T denominator; 435 of those union plasmids are not in the MPF_T
set. Numerator and denominator were over different universes.

## Why the A9 seal did not catch this

A9 verifies that a curated file has not CHANGED. It cannot verify that a number
in it was ever derivable. `anchor_set.tsv` was sealed with a threshold whose
justification could not be reproduced from the cache -- and the seal made that
look settled. A seal is an integrity guarantee, not a correctness one.

**Proposed guard (A13): every threshold in the anchor set must be re-derivable
from the cache by script, and any `E`-based threshold is rejected at load.**

## Outcome of the wider test (job 25913461)

Five MPF_T anchors were tested against their CONJScan counterparts under the same
two controls. **Only one disjunction besides VirB5 is admissible.**

| anchor | GA fail | union | admissible | why |
|---|---|---|---|---|
| VirB1 | 19.5% | 13.2% | **yes** | band 95.4%, cross-class 8.5% -- stays peripheral |
| VirB2 | 20.9% | 5.2% | NO | **cross-class 33.1% MPF_F** |
| VirB3 | 9.0% | 5.0% | NO | **cross-class 57.4% MPF_F** -- a majority |
| VirB5 | 25.4% | 7.5% | **yes** | band 91.3%, cross-class 6.6% |
| VirB6 | 4.9% | 4.3% | NO | **band 15.6%** -- slot-defining, would corrupt coordinates |

VirB2 and VirB3 appear to collapse, but their added plasmids are heavily F-type.
`T4SS_T_virB3` admitting a MAJORITY of MPF_F plasmids is not recovery of divergent
P-type VirB3. Both unions are NOT read; GA stands.

**The two-tier structure therefore survives this challenge**, and the reading that
pilus-associated components are harder to detect does not need withdrawing. But
the reason VirB2/VirB3 did not move is that their disjunctions were REJECTED on
purity, not that CONJScan agreed with Pfam. That distinction matters: the raw
union numbers look like a collapse and would have been reported as one without the
cross-class control.

**VirB1 is why the rejections are readable.** It passed both controls and still
stayed peripheral, which rules out "CONJScan profiles are simply more permissive"
as a blanket explanation for VirB2/VirB3. An expectation of no change is not a
control; running it anyway is what made the other three interpretable.

## A second safeguard written and not applied

The first version of the test script (job 25912704) computed its verdict from the
failure rate ALONE and never consulted the `in band` or `cross-class` columns it
had just calculated. It printed "now core-level" for VirB2, VirB3 and VirB6 --
two of which fail cross-class and one the band.

This is the same failure already recorded in `ga_threshold_failure.md`: the guard
existed, was computed, and was not consulted. It was reproduced inside the script
written to test that very file. The verdict logic now refuses to read a union from
a disjunction that failed either control, and writes `disjunction_admissible` and
`verdict` into the output TSV so a rejected union cannot be picked up downstream.

## One control was vacuous and must not be counted as a pass

Measured robust bands from each anchor's own Pfam hits:

| anchor | band (aa) | median | in-band | usable? |
|---|---|---|---|---|
| VirB1 | 76-304 | 186 | 95.4% | yes |
| VirB2 | 60-156 | 101 | 94.3% | yes |
| **VirB3** | **1-2139** | 106 | **100.0%** | **NO -- vacuous** |
| VirB5 | 198-290 | 238 | 91.3% | yes |
| VirB6 | 274-426 | 351 | 15.6% | yes (and fails) |

A band of 1-2,139 aa admits every protein in the cache, so VirB3's 100.0% is the
filter NOT FIRING, not the filter passing. Same defect as VirB10's 43-857 band.
PF05101 evidently hits as a domain inside proteins of wildly varying length.

VirB3 is still REJECTED -- cross-class 57.4% is decisive on its own -- but that
rejection rests on ONE control, not two. Recorded so the row is not later read as
having cleared both.

**Zero rejection by a filter is evidence the filter did not fire**, not evidence
of purity. A band control must be checked for degeneracy before its result is
read; a band wider than the cache's own length range carries no information.

## Reseal (2026-09-14)

`data/anchors/anchor_set.tsv` VirB5 row:

* `threshold` `E<=1e-5 AND aa>150` -> **`GA AND aa>150`**
* rationale replaced with the measurement above; the withdrawn rule and the reason
  it was withdrawn are kept IN the note, not deleted, so a reader cannot come away
  believing IncP TrbJ was ever rescued by an E-value rule
* header vocabulary updated; a standing prohibition on E-value thresholds added,
  to be enforced by A13
* measured at the new threshold: **5,545 of 7,436 = 25.4% failure**
* no binary tier label assigned to VirB5 (7.5% with the disjunction falls between
  the two bands); carry the failure rate itself

Resealed under A9. Pre-edit copy retained in the session scratchpad.

## A13 landed (2026-09-14), and where the withdrawn rule still lives

`assertions.assert_no_evalue_thresholds(proj)` rejects any E-value in the
`threshold` column of an anchor set. `decontaminate.assert_no_evalue_anchor()`
does the same for `ANCHOR_FAMILIES` in code.

Self-tested both directions, which is the discipline that was missing when the
rule reached v1: A13 PASSES the corrected `anchor_set.tsv` and REJECTS the
pre-reseal copy, naming the offending row.

`scripts/lib/decontaminate.py` corrected: `("VirB5", "anchor_T4SS", 1e-5, 150)`
-> `("VirB5", "anchor_T4SS", None, 150)`, plus a new `ANCHOR_DISJUNCTION` table
carrying the adopted CONJScan arms (VirB5, VirB1 only). **A consumer that scans
ANCHOR_FAMILIES alone now gets the GA arm and undercounts VirB5 by 1,336
plasmids** -- 5,545 rather than 6,881. That is stated at the definition.

### Still carrying the withdrawn rule -- historical, NOT corrected

Seven analysis scripts hardcode `("anchor_T4SS", 1e-5, 150)`:

`91_gap_occupant_lengths.py`, `92_slot_sixframe_recheck.py`,
`98_negative_slot_control.py`, `99_architecture_call.py`,
`100_inslot_calibration.py`, `101_round0_final.py`, `126_virb5_swap_purity.py`

These are left AS THEY RAN, deliberately: editing them would break reproduction
of results already reported. But every VirB5-conditional figure they produced was
computed under the withdrawn rule and is therefore an UNDERCOUNT. Specifically:

* the six-frame slot recheck (92) and the occupancy figures derived from it
* the architecture fallback validation at 98.1% (99)
* the in-slot calibration (100)
* the negative slot control (98)

None of these is re-derived here. **Any of them reused for v1.1-era conclusions
must be re-run against the disjunction first.** Recorded as outstanding, not
quietly inherited.
