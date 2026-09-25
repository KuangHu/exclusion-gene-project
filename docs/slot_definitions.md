# Two slot definitions, reported side by side

They measure different things and must never be merged into one "occupancy"
number. Both are kept; each carries its conditioning in its label.

| label | definition | denominator |
|---|---|---|
| `slot_occupancy \| VirB5 detectable` | is there a gene between the VirB5 and VirB6 anchors | plasmids where **VirB5** is detectable |
| `neighbour_is_unassigned \| VirB6 detectable` | what is the gene immediately 5' of VirB6 | plasmids where **VirB6** is detectable |

## Why the second one exists

The GO on the core hypothesis came from one observation: all four MPF_T controls
carry a small lipoprotein **immediately 5' of VirB6**. "Between VirB5 and VirB6"
is a proxy for that, and it is only equivalent when the operon is in canonical
order. Removing the VirB5-before-VirB6 assumption recovered 554 plasmids where
VirB6 precedes VirB5; in those, the between-anchor interval is not the
5'-of-VirB6 position at all.

## The mistake this file exists to prevent

A first draft of `95_upstream_of_virb6.py` proposed dropping VirB5 entirely and
asking only "is there a gene 5' of VirB6". That is a tautology -- there always is
one, unless VirB6 sits at a record boundary. It would have produced ~100%
occupancy, and under v2 6.1 a ~100% occupancy retires the whole presence/absence
direction. A definitional identity would have killed a research direction.

This is the same shape as the `\r` bug: a non-biological cause producing a number
that lands exactly on a decision threshold. **Requiring VirB5 is what made "empty"
definable.** Dropping it did not remove the conditionality, it removed the
measurement.

So "empty" is decided by the neighbour's IDENTITY, not its existence:

| upstream neighbour hits | verdict |
|---|---|
| a canonical T4SS component (VirB1-11, VirD4, TraN) | `slot_empty` |
| a known Eex family | `eex_occupied` |
| nothing | `candidate` |

The requirement is now spread over thirteen component families rather than
resting on PF07996 alone, whose 21.4% failure rate is an order of magnitude worse
than VirB6/VirB8/VirB9 at 4.9%.

**Red line (v2 1.3) holds.** The Eex families label a gene that position has
already selected. They never bound a window, never pick an anchor pair, and never
decide which gene is examined.

## Validation 1 -- the inverted class is real, not a strand bug

If the pairing code were not strand-aware, "inverted" would come out near 50%,
since plasmid operons face either way relative to the deposited sequence. Measured:
**554 / 5,149 = 10.8%**. The code is strand-aware and the inverted plasmids are
genuine rearrangements. Independently, in six separate Mash clusters both anchors
score above GA (VirB5 213-215 aa at 34-43; VirB6 328-331 aa at 120-137) and a full
Pfam-A scan hits nothing else on either protein.

## Validation 2 -- the definitions agree exactly where the proxy is valid

1/16 deterministic-hash shard, 348 plasmids (A8: hash sampling, not contiguous):

| order | n | old verdict | new verdict |
|---|---|---|---|
| canonical | 260 | OCCUPIED | 212 `eex_occupied` + 45 `candidate` + 3 `slot_empty` |
| canonical | 51 | EMPTY | 48 `slot_empty` + 3 `candidate` |
| inverted | 36 | EMPTY | 35 `candidate` + 1 `slot_empty` |
| inverted | 1 | OCCUPIED | 1 `candidate` |

Concordance on canonical plasmids: **305/311 = 98.1%**. Complete divergence on
inverted plasmids: the old definition called 97.8% of them empty, the new one
calls 97.3% of them occupied by an unlabelled gene.

That is the signature of a proxy that is valid under canonical order and invalid
under inversion -- not of two competing measurements. The 8.6-point drop in record
occupancy (82.4% -> 73.8%) when the inverted class was admitted was a proxy
artefact, not a real occupancy difference.

## Standing caveat

`slot_occupancy | VirB5 detectable` should be quoted **only for canonical
plasmids**. Applying it to inverted ones measures the wrong window. It is retained
because it is the figure every earlier result was computed against.

## Full-scale result (7,473 rows over 7,014 accessions)

Labelled against 13 anchor families (VirB1-11, VirD4, both TraN models) and 6 Eex
families. `candidate` was never provisional on five anchors -- the full set was
used from the first run.

| verdict | records | Mash clusters | unique sequences |
|---|---|---|---|
| `slot_empty` | 974 (13.0%) | 346 (11.8%) | 320 (16.1%) |
| `eex_occupied` | 3453 (46.2%) | 1169 (39.9%) | 442 (22.3%) |
| `candidate` | 3046 (40.8%) | 1511 (51.5%) | 1221 (61.6%) |

**Report all three levels or none.** The named Eex families hold 46% of records but
only 22% of distinct proteins; `NF033894` alone is 3,236 records, 1,047 clusters,
396 unique sequences (unique/record 0.13, three times more redundant than the
candidate pool at 0.40). A record-level reading says the slot is mostly covered by
known families; a sequence-level reading says the opposite.

### By anchor order

| order | records | slot_empty | eex_occupied | candidate |
|---|---|---|---|---|
| canonical | 5080 | 841/292/261 | 3355/1123/411 | 884/550/485 |
| inverted | 589 | 19/11/11 | 0/0/0 | 570/115/30 |

*(records/clusters/unique sequences)*

The 570 inverted candidate records are **30 distinct proteins** -- one candidate
family, not 570 discoveries. Zero inverted plasmids carry a named Eex family
against 66% of canonical ones; that asymmetry is unexplained and is not to be
quoted as a result until it is.

## Known limits of the `candidate` class

`candidate` means "hit nothing", so it inherits the detection failure of every
family used to label it. From the completed ladder: VirB1 19.5%, VirB2 20.9%,
VirB5 27.1%, and VirB7/TraN 97.2%/97.7%. Some candidates are unrecognised
canonical components, and the VirB7-class ones are invisible for the reason the
retraction established -- the caller does not emit 48 aa ORFs. The error therefore
inflates the candidate pool rather than the occupied one.

## What the VirB5 conditioning costs, measured

Over the 7,436 VirB4+ plasmids, the between-anchor definition can only be
evaluated where a VirB5/VirB6 pair exists:

| class | n | % of examined |
|---|---|---|
| pairable (both anchors, same strand, <=6 genes apart) | 5493 | 73.9% |
| `no_VirB5_only` | 1286 | 17.3% |
| `no_VirB5_and_no_VirB6` | 307 | 4.1% |
| `gap_gt_6_genes` | 259 | 3.5% |
| `no_VirB6_only` | 60 | 0.8% |
| `opposite_strand_only` | 31 | 0.4% |

**SUPERSEDED IN PART -- see `docs/slot_anchor_is_virb5.md`.** The reasoning below
treats the VirB6-anchored definition as a clear gain. It is not: the exclusion
gene tracks VirB5 (84.8% pooled) rather than VirB6 (72.8%), and the VirB6 rule is
a proxy valid only under canonical order. The coverage figures stand; the
conclusion drawn from them does not.

**VirB5 absence causes 1,593 of the 1,943 exclusions -- 82%.** The
neighbour-identity definition needs only VirB6 and covers 7,014 accessions
against 5,493: **+1,521 plasmids, a 28% larger denominator**, bought by depending
on a 4.9%-failure anchor instead of a 27.1%-failure one.

This is the concrete reason both definitions are kept rather than one being
retired: they are not the same measurement on the same population.

## Correction: the candidate pool was contaminated by unrecognised VirB5

The cross-tab exposed it. Among rows where no VirB5 was detected anywhere on the
plasmid, **87.2% were called `candidate`**, against **17.4%** of canonical rows --
a 5x difference that is not plausibly biological. On a plasmid where VirB5 is
undetectable, the gene 5' of VirB6 is very often the unrecognised VirB5 itself.

Two tests, and the first one was inconclusive for a reason worth recording:

1. **PF07996 at a permissive threshold: useless here.** Zero of the 1,221 unique
   candidates reach GA. That cannot settle anything -- a VirB5 the model cannot
   see is invisible to the model by construction, which is exactly the 27.1%
   failure rate. (A first pass compared E-values across databases of different
   size, which is invalid; bit scores are the database-independent quantity.)

2. **phmmer against the 228 confirmed VirB5 proteins: decisive.**

| candidate length band | unique | phmmer E<=1e-3 vs real VirB5 |
|---|---|---|
| <60 | 83 | 5 (6.0%) |
| 60-100 | 486 | **0 (0.0%)** |
| 101-150 | 353 | 12 (3.4%) |
| 151-199 | 56 | 2 (3.6%) |
| **200-288 (VirB5 band)** | 228 | **162 (71.1%)** |
| >288 | 15 | 5 (33.3%) |

The 60-100 aa band is a clean negative control at 0/486. PF07996 misses these
completely -- not marginally -- and phmmer recovers them, the same pattern that
decided the tool choice earlier.

**186 unique proteins / 390 records reclassified** as
`slot_empty_virb5_recovered`.

| verdict | records | clusters | unique |
|---|---|---|---|
| `slot_empty` | 974 | 346 | 320 |
| `slot_empty_virb5_recovered` | 390 | 217 | 186 |
| `eex_occupied` | 3453 | 1169 | 442 |
| `candidate` | 2656 | 1327 | **1035** |

Candidate pool 1,221 -> **1,035 unique (-15.2%)**. Filtering by length alone would
have cut 228 (-18.7%) and discarded 66 proteins that sit in the VirB5 length band
with no VirB5 homology by either method.

**The conclusion survives the correction**: named Eex families still hold 46% of
records but only 442 of 1,477 non-empty unique proteins. The candidate pool is
still the larger half.

### The inverted asymmetry is NOT VirB5 contamination

Zero inverted plasmids carry a named Eex family against 66% of canonical ones. The
natural explanation -- that inverted "candidates" are really VirB5 -- is refuted:
the 30 unique inverted candidates have median length **85 aa** (Q1 84, Q3 102),
only 2 sit in the 200-288 aa band, and **0/30 hit any real VirB5 by phmmer**. They
are small proteins in the exclusion size band that no named family recognises.

**Followed up and DEMOTED.** The inverted class is 85% IncI2 / IncI2(Delta) and
97.6% MPF_T -- more MPF_T than the canonical set (87.4%), so they are not
mistyped. Gene order is identical across every IncI2 plasmid examined:

    ... 148aa | 85aa | 100aa | [85aa slot protein] | VirB6 (TrbL) | VirB5 (214aa) | 74aa ...

VirB5 sits **3' of VirB6**, not 5'. The relation "gene immediately 5' of VirB6" is
preserved across the rearrangement while VirB5's position is not -- direct
architectural support for the neighbour-identity definition over the
between-anchors proxy.

But the 85 aa occupant FAILS the architectural test for an exclusion protein:

| | RP4 trbK | R751 trbK | pKM101 eex | R388 cand | IncI2 85 aa |
|---|---|---|---|---|---|
| lipobox | LAG-C23 | VAG-C22 | LVA-C15 | LSA-C17 | **absent** (1/10 unique, 1/499 records) |
| max hydrophobicity (KD w=19) | +1.97 | +1.96 | +1.45 | +1.04 | **+0.83** |
| Pfam-A at GA | | | | | none |
| ExcA homology (R64 220 aa, R621a 204 aa) | | | | | no hit at E<=100 |

N-terminus `MKKNDDLITLPSHNSHALDYIIKKAFTPRK...` is charged and hydrophilic, with no
signal peptide. Against a slot-wide lipobox rate of 67.8% vs 0.5% in matched
control gaps, this sits on the wrong side of the strongest discriminator
available. 471 of 499 records are the identical sequence, which fits a conserved
IncI2 gene of some other function better than a divergent exclusion protein.

**Positional candidate, architectural mismatch. Not a discovery lead.** The
literature framing that made it attractive -- IncI1/IncI-gamma carry excA/excB
targeting TraY on an MPF_I system, so an MPF_T IncI2 would be structurally
independent -- makes the position interesting but is not evidence about this
protein, and the direct tests are negative. The 0%-vs-66% Eex asymmetry therefore
remains unexplained.
