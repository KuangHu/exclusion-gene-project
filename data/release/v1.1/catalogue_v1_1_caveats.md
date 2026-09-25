# Release v1.1 -- MPF_T catalogue, corrected VirB5 rule

**v1 is frozen and unedited.** v1.1 is a corrected rebuild. v1's numbers remain
reproducible under v1's own (wrong) rule -- that is what the freeze is for. This
file records the delta.

## The one change

| | v1 | v1.1 |
|---|---|---|
| VirB5 detection | `PF07996 @ E<=1e-5 AND aa>150` | `(PF07996 @ GA AND aa>150) OR (T4SS_T_virB5 @ GA AND aa>150)` |

Nothing else in the detection stack changed.

## Why the v1 rule was withdrawn

Full account: `docs/virb5_evalue_threshold_invalid.md`. In short, three
independent defects, any one of which is disqualifying:

1. **It failed its own controls.** The rule existed to rescue IncP TrbJ, which GA
   loses. Under the pipeline's own hmmsearch (job 25913366) RP4 trbJ scores
   17.7 bits **E = 0.22** and R751 trbJ 22.2 bits **E = 0.0093**. Both FAIL
   `E<=1e-5`. The rescuing E-values in the rationale (7.9e-07, 3.4e-08) came from
   a protein-vs-Pfam scan (N ~ 30,134 models), not a model-vs-proteome search
   (N = 1,099,911). The bitscores matched exactly; only the E-values did not.
2. **It was backwards.** `E = 1e-5` falls at ~31.7 bits in this search, ABOVE the
   GA of 24.7. It was documented as a loosening and was a tightening, costing
   ~177 plasmids (GA 5,545 vs E-rule 5,368).
3. **It was not reproducible.** E-values scale with database size: the same rule
   on the same sequences gives 5,368 in a 1.1M-protein search and 5,341 in a
   7.7M-protein one. Now prohibited outright (A13).

## Measured delta (job 25916674)

| | v1 | v1.1 |
|---|---|---|
| VirB5 present | 5,424 (27.1% failure) | **6,881 (7.5% failure)** |
| `slot_ready__virb5_virb6` | 5,272 (70.9%) | **6,548 (88.1%)** |
| plasmids whose VirB5 presence flipped | -- | **1,457** |
| `tier` distribution | -- | **identical (asserted)** |
| `architecture` | -- | **identical (asserted, build aborts on drift)** |

VirB5 call provenance, new column `virB5_source`:

| source | plasmids |
|---|---|
| `pfam_only` | 3,468 |
| `both` | 2,077 |
| `conjscan_only` | **1,336** |
| `none` | 555 |

`conjscan_only` is the disjunction's actual contribution: 1,336 plasmids PF07996
misses at GA entirely. RP4 and R751 are among them.

## What did NOT move, and why that was checked rather than assumed

`tier` and `core_completeness_ga` use
`CORE = VirD4, VirB4, VirB6, VirB8, VirB9, VirB10, VirB11`, which excludes VirB5.
`architecture` is called from VirB4/VirB6 strand order only. Both are therefore
VirB5-invariant *by construction* -- but "by construction" is exactly the kind of
claim that was wrong about the E-values, so the builder recomputes architecture
per plasmid and **aborts on any difference**, and compares the full tier
distribution and **aborts if it moves**. Both passed.

## Columns added in v1.1

* `virB5_source` -- `pfam_only` / `conjscan_only` / `both` / `none`
* `slot_ready__v1` -- v1's value, carried alongside, so the delta is auditable
  per plasmid and not only in aggregate

## Figures previously reported that this release corrects

| figure | previously stated | v1.1 |
|---|---|---|
| VirB5 failure under the disjunction | 27.1% -> 3.1% | **27.1% -> 7.5%** |
| `slot_ready` | 70.9% -> 87.2% | **70.9% -> 88.1%** |
| "VirB5 enters the core tier" | asserted | **no tier label** -- see below |

The `3.1%` and `87.2%` came from dividing a union counted over all 72,556 PLSDB
plasmids by the 7,436-plasmid MPF_T denominator. 435 union plasmids are not in the
MPF_T set. Numerator and denominator were over different universes.

## Anchor tiers are no longer binary

`data/anchors/anchor_tiers_v1_1.tsv` replaces the core/peripheral label with the
measured failure rate. VirB5 at 7.5% falls between the tightest group (4.1-4.9%)
and the next (9.0-20.9%); any cutoff between them is arbitrary and a binary column
would discard the number that is actually measured. The `group` field there is a
descriptive annotation and **must not be branched on** -- code needing a cut
declares its own and says so.

## Scope: this release corrects MPF_T only

`catalogue_MPF_F_v1.tsv` and `annotation_MPF_T_F_v1.tsv` are NOT reissued here.
The MPF_F entry criterion (`PF11130 TraC_F_IV`) is GA-based and unaffected by the
E-value defect. However `docs/mpf_f_entry.md` argues against a `TrbC_Ftype` entry
filter partly from the binary two-tier structure, and that argument's wording is
now stale even though its conclusion stands on independent grounds (the VirB4-ATPase
convergence across three classes, and 5/5 seeds). **The merged annotation still
carries v1's VirB5 column** and should be regenerated before use in any
VirB5-conditional analysis. Flagged, not silently reissued.

## Not carried over from the wider anchor test

Four of five anchor disjunctions tested in job 25913461 were REJECTED on purity
and are NOT applied in v1.1:

| anchor | rejected on |
|---|---|
| VirB2 | cross-class 33.1% MPF_F |
| VirB3 | cross-class 57.4% MPF_F (band control was vacuous, 1-2139 aa) |
| VirB6 | length band 15.6% |

Their apparent failure-rate collapses (VirB2 5.2%, VirB3 5.0%) are artefacts of
importing F-type systems and are not recorded as failure rates anywhere.
