# Iteration module

## Round 0 is already complete

It was not called that, but the machinery exists and is calibrated:

    position nomination -> phmmer all-vs-all -> Leiden -> 63 families with no
    named Eex member

with the threshold validated by seed control (E<=1e-05; the 1e-07 and 1e-09
"stricter" choices split the RP4/R751 TrbK pair, whose true edge is 6.0e-07),
percolation excluded by ground truth rather than by size statistics (the two seed
families never merge at any resolution tested; connected components at the same
threshold puts 80.3% of the input in one component), and family extension = 22
stable across a 20-fold resolution range.

So the question is not whether iteration can start. It is whether Round 1 should
run, and that depends on the precision controls.

## Precision controls: status

| # | control | status |
|---|---|---|
| a | parallel negative slot | **being built** (`98_negative_slot_control.py`) |
| b | within-cluster presence/absence | **DEAD** -- only 8 mixed clusters |
| c | seed freezing + drift measurement | specified below, not implemented |
| d | out-of-slot occurrence rate | **is Round 1's main output** |

(b) is unrecoverable. (a) and (d) are required.

## What Round 1 is actually for

Not expansion. **Validation.**

The candidates were nominated by POSITION. Building an HMM per family and
searching the whole database finds their homologs OUTSIDE the slot, which is the
test the nomination could not perform on itself:

| in-slot / total | reading |
|---|---|
| ~1 | slot-specific: a genuinely position-defined family |
| clearly < 1 | a general small-protein family; the positional nomination caught coincidences |

This criterion does not depend on the dead channel (b) and is not affected by
percolation. Finding additional members is a by-product.

## BLOCKER 1 (resolved): the nomination set was biased

Nominating "eex at +1 from VirB5" silently drops every plasmid where PF07996 fails
-- 27.1% of the VirB4+ set, and non-random: it enriches for divergent systems,
i.e. exactly where a novel family would be.

**Architecture can be called without VirB5.** In canonical order VirB4 precedes
VirB6; in the rearranged order VirB6 precedes VirB4 (`VirB6 VirB5 eex VirB3
VirB4 ...`). The sign of (VirB6 - VirB4) in transcription order therefore calls it
using VirB4 (0% failure -- it is the entry criterion, present by construction) and
VirB6 (4.9%). The eex offset then follows from VirB6 alone: -1 canonical, +2
rearranged.

Validated against the plasmids where VirB5 IS detectable, so the architecture is
already known independently:

| | called correctly |
|---|---|
| known canonical | 4808 / 4875 = **98.6%** |
| known inverted | 506 / 546 = **92.7%** |
| overall | **98.0%**, using no VirB5 |

Applied to the 1,699 VirB6-positive / VirB5-negative plasmids: 82.1% canonical,
7.7% rearranged, 10.2% uncallable.

**One correction was needed.** In canonical architecture eex is at -1 from VirB6
and VirB5 at -2, so when the slot is EMPTY the gene at -1 is VirB5 itself -- and on
these plasmids VirB5 is by definition undetected. The first pass nominated it:
length median 98, Q3 255, straight into the VirB5 band. phmmer against the 228
confirmed VirB5 proteins removed 190 unique / 387 records:

| | records | unique | median aa | lipobox (record, strict) | 60-100 aa |
|---|---|---|---|---|---|
| before filter | 1509 | 752 | 98 | 36.8% | 45.7% |
| after filter | **1122** | **562** | **77** | **43.3%** | **56.9%** |

Content check, record-level strict lipobox: main VirB5-nominated set 51.1%, this
recovered set 43.3%, VirB10 negative slot 0.7%. The recovered set sits with the
target, not the null.

**Round 1's nomination set gains 562 unique proteins** from the subpopulation most
likely to carry divergent families.

## Before Round 1 runs

1. **The negative slot** (running). Same code path at VirB8/VirB9/VirB10 -- anchors
   with comparable detection failure (4.9%/5.0%/4.7%), GA-level families, and no
   exclusion gene reported adjacent to them. If the negative slots also yield
   dozens of families, 63 is this pipeline's baseline yield between any conserved
   gene pair, not a signal. **This is the only thing that gives 63 a denominator.**

2. **Composition of the 63.** They are currently only known to be "not hit by a
   named Eex family" -- they could be transposases or resistance genes. Stratify by
   lipobox rate, the strongest discriminator available (67.8% in slot occupants
   against 0.5% in matched control gaps). Expect two modes; only the high one is a
   candidate set.

3. **Anchor rule fixed** -- DONE, see `slot_anchor_is_virb5.md`. The nomination
   step depends on it. Architecture-first: core anchors establish canonical vs
   rearranged, then eex at +1 from VirB5 (equivalently -1 from VirB6 in canonical).

## The negative slot, and what it forced

Same code path at anchors with no exclusion-gene hypothesis (`data/anchors/
negative_slot_control.tsv`):

| anchor | unique nominated | named Eex | families | lipobox>=50% families | overall lipobox |
|---|---|---|---|---|---|
| **VirB5 (target)** | 1098 | **372** | 35 | **13** | **51.1%** |
| VirB8 | 1950 | 0 | 11 | 0 | 5.6% |
| VirB9 | 2102 | 0 | 19 | 3 | 17.2% |
| VirB10 | 1890 | 0 | **62** | 3 | **0.7%** |

**Family count is not a signal, and is ANTI-correlated with one.** VirB10 yields 62
families -- nearly twice the target -- with the lowest lipobox rate measured. Any
report of "N novel families" describes this pipeline's baseline yield between two
conserved genes, not a discovery. The earlier figure of 63 (VirB6-based nomination
over 7,014 accessions) is not even comparable to the negatives; the comparable
target number is 35.

### Occupancy separates before any family count

Measured on adjacent anchor pairs, all orders (`negctrl_occ`):

| pair | role | pairable | occupancy |
|---|---|---|---|
| **VirB5\|VirB6** | **TARGET** | 5101 | **73.4%** |
| VirB8\|VirB9 | negative (adjacent) | 7002 | **0.7%** |
| VirB9\|VirB10 | negative (adjacent) | 6980 | 15.2% |
| VirB10\|VirB11 | negative (adjacent) | 5768 | 20.9% |

VirB8 and VirB9 are immediately adjacent in 6,954 of 7,002 plasmids. So "a gene
sits between two adjacent skeleton components" is NOT generic -- it is specific to
VirB5|VirB6 by ~100x. This is the earliest-readable evidence in the whole
analysis, and it does not depend on clustering, thresholds, or family counts.

**Two distant pairs (VirB3|VirB9, VirB4|VirB11) were also run and are USELESS.**
Both returned 100% occupancy, necessarily: their anchors are separated by known
skeleton components, so genes always sit between them. The control measured "are
there genes between two genes that are far apart". Recorded so it is not repeated.

The remaining gap: VirB8/9/10/11 form one contiguous stretch of skeleton, so their
agreement could mean "stable baseline" or "this stretch is uniform". A genuine
independent baseline needs an adjacent pair from a DIFFERENT functional module --
two Rep, Par or MOB genes that neighbour each other in their own operon. Not yet
built; it requires anchor families outside the T4SS set.

What separates the target is CONTENT:

* **372 named-Eex members at VirB5, zero at all three negatives.** The known
  families sit at this one skeleton position and nowhere else along it.
* **lipobox 51.1% vs 0.7-17.2%**, and 13 lipobox-rich families vs 0-3.

## Stop rules

The originally specified rule -- halt when negative-slot family count >= 50% of
target -- **fires immediately at Round 0** (11/35 = 31%, 19/35 = 54%, 62/35 =
177%) and is therefore wrong, not conservative. It was written on the assumption
that family count carries signal; the control shows it does not.

Restated on the quantities that do separate. Halt when ANY of:

* new members in a round = 0
* median identity of new members to the **r0 seeds** < 25% (inter-family background
  is 21-25%, so below this the round is recovering background)
* **lipobox-positive family count at the target falls to <= 2x the negative-slot
  maximum** (currently max 3, so the trigger is <= 6 against Round 0's 13)
* **a round's new members are majority non-lipobox** -- the criterion separating
  slot from matched control by 9-20x
* round > 3

## Evidence tiers

`E4_r0`, `E4_r1`, `E4_r2`, ... **Never merged, never upgraded.** A member found in
round 2 is permanently marked as such. This mirrors the existing rule that
`E4_positional` never upgrades: provenance is a property of the assertion, not a
stage it graduates from.

## Seed freezing

New members do NOT become seeds for the next round. Every round searches with
models built from the **r0 seeds only**, and reports the identity distribution of
its new members against those r0 seeds.

Without this, each round's model drifts toward whatever the previous round
admitted, and after three rounds the family is defined by the search rather than
by the seeds -- with no way to measure how far it moved. Freezing makes drift
observable: it is exactly the reported identity distribution.
