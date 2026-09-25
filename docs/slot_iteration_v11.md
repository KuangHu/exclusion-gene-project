# MPF_T trbJ+1 slot iteration on the v1.1 nomination set

Jobs 26048996 (nominate+decontaminate), 26051683 (Round 0), 26054729-33 (Round 1),
26067590 (three-layer counts). Negative slots traversed the identical path.

## Two calibration values from earlier work are NOT usable in this run

### 1. `in-slot/total` is not on the documented scale

The Eex_IncN cluster -- which CONTAINS pKM101 eex and R388, i.e. a known true
exclusion family -- scores **0.095** here. Its documented reference (NF033894) is
**0.771**.

The denominator differs: this run counts hits across all 72,438 cached plasmids;
the reference was measured on a smaller universe. **The calibration values 0.745
(null max), 0.771 and 0.812 (upper) therefore CANNOT be compared to any number in
this run.** Only target-vs-negative, measured identically, is valid.

Anyone reaching for 0.745 as a threshold on these numbers will be wrong.

### 2. The seed-control pass band is narrower than documented

The instruction stated a pass band of "1e-03 to 1e-05 x resolution 0.005-0.1,
spanning two orders of magnitude, so threshold choice is not load-bearing".
Measured on this 681-sequence pool:

| E | res 0.005 | 0.01 | 0.05 | 0.1 |
|---|---|---|---|---|
| 1e-03 | PASS | PASS | PASS | PASS |
| 1e-05 | PASS | **PASS** | FAIL eex_split | FAIL eex_split |
| 1e-07 | FAIL trbK_split | FAIL | FAIL | FAIL |

Resolution 0.05 -- the specified value -- splits the Eex_IncN pair. **Adopted
1e-05 x 0.01.** E was held at 1e-05 because Round 1 must search at the clustering
threshold.

Cause: seed separation is TIGHTER on this background than documented (within
TrbK_RP4 2.3e-07 vs 6.0e-07; within Eex_IncN 1.8e-09 vs 4.8e-09). The pool is 681
sequences, not the 4,410 the original separation was measured on, and E-values
scale with database size -- the same property that invalidated the VirB5 anchor
threshold (A13). **Every E-value cutoff here is a property of the pool it was
measured on.**

## The sentinel caught a criterion that would have destroyed the result

`in-slot/total` does not merely have a thin margin. **It REJECTS the known
Eex_IncN family at 0.095**, whose lipobox is 96.2%.

Diagnosis: a successful exclusion family is widespread across the database, so its
in-slot fraction is low BY CONSTRUCTION. The criterion selects for RARITY, not for
function.

Without the positive control, the conjunction `in-slot/total > 0.745 AND
lipobox > 50%` yields cl39 and cl71 at the target and **zero at all three negative
slots** -- a clean-looking specificity result that would have been reported. One
of its two conjuncts excludes the family we already know is real.

### Using a rejected criterion for extreme-value exclusion is a different act

cl34 is set aside as a generic lipoprotein family: 4,576 hits DB-wide, 66 in-slot,
ratio 0.014.

**That reasoning uses `in-slot/total`, which was just rejected.** The distinction:
0.014 is 7x below the known family's 0.095, with only 66 of 4,576 hits at the
slot. Using the axis to exclude an EXTREME OUTLIER is defensible; using it as a
THRESHOLD is not. Stated explicitly so the rejected criterion is not quietly
reinstated.

## What does not work, confirmed

| method | evidence |
|---|---|
| family count | VirB10+1 gives **50** novel families vs the target's **40** (previously 62 vs 35). Second independent confirmation of anti-correlation. Any "N novel families" claim is baseline yield. |
| `in-slot/total` alone | null reaches 1.000; rejects the known family at 0.095 |
| ratio 1.000 as evidence | the 12 target families at ratio 1.0 have **median hits_total = 5**. Trivially 1.0 on rare families. |

## What works: lipobox alone, which passes the positive control

Families with hits_total >= 20:

| cut | TARGET | VirB8+1 | VirB9+1 | VirB10+1 |
|---|---|---|---|---|
| >50% | 7 | 0 | 1 | 1 |
| >75% | 6 | 0 | 1 | 0 |
| **>90%** | **5** | **0** | **0** | **0** |

At >90%: five target families, **zero in all three negative slots**, and the known
Eex_IncN family is among the five -- recovered, not missed.

| cluster | r0 uniq | hits | in-slot | ratio | lipo% | med aa | identity |
|---|---|---|---|---|---|---|---|
| 39 | 10 | 28 | 27 | 0.964 | 100.0 | 66 | novel |
| 46 | 13 | 260 | 144 | 0.554 | 98.8 | 108 | novel |
| **25** | 88 | 1755 | 166 | 0.095 | 96.2 | 77 | **Eex_IncN, positive control** |
| 34 | 28 | 4576 | 66 | 0.014 | 94.4 | 84 | novel, DB-wide -- set aside |
| 13 | 57 | 159 | 111 | 0.698 | 91.8 | 114 | novel |

**Result: four novel families, one of which (cl34) is set aside, plus recovery of
the known family.** Three defensible novel families: cl39, cl46, cl13.

## Caveats on the three

* **cl39 has only 28 hits**, barely over the >=20 line, and this run established
  that small families score trivially high on any in-slot ratio. RESOLVED by the
  Mash layer below: 16 independent clusters, not one lineage.
* **cl46 (108 aa) and cl13 (114 aa)** sit above the 70-85 aa band of the known
  MPF_T occupants (Eex_IncN median 77, cl34 84). Not grounds for exclusion, but
  they may be a different kind of thing.
## Three-layer counts (§8), job 26067590

Mash single linkage at d <= 0.007.

| family | records | Mash clusters | unique | collapse | identity |
|---|---|---|---|---|---|
| cl39 | 27 | **16** | 10 | 1.7x | novel |
| cl13 | 106 | **69** | 58 | 1.5x | novel |
| cl46 | 108 | **46** | 13 | 2.3x | novel |
| cl25 | 173 | **135** | 92 | 1.3x | Eex_IncN, positive control |
| cl34 | 65 | 57 | 31 | 1.1x | set aside |
| **union** | **469** | **315** | **198** | 1.5x | |

**Nothing is clonal.** Collapse is 1.1-2.3x, against IncI2's 181 -> 25 (7x). The
cl39 concern does not hold: its 27 records span **16 independent Mash clusters**.

**cl46 is the strongest signal**: 108 records over 46 Mash clusters but only 13
unique sequences -- a near-invariant protein on diverse plasmid backbones, the
signature of functional constraint rather than clonal expansion.

The known family behaves as a real family should (173 / 135 / 92), which is an
independent check that the counting is sane.

NOTE ON UNITS: an earlier figure of 514 in-slot records was PROTEIN-LEVEL hits
from Round 1. The accession-level union is **469**. Both correct, different units;
469 is the figure to pair with the Mash counts.

d = 0.007 achieves DEDUPLICATION ONLY, not independence -- the graph percolates at
d = 0.016-0.018. Cluster-level hypergeometric tests still violate independence;
cross-family enrichment is valid only at PTU/Inc-group level.

## Decontamination is where the signal lives

All four slots nominate 1,618-1,856 unique sequences. After removing anchors:

| slot | nominated unique | kept unique | kept records |
|---|---|---|---|
| VirB5+1 | 1,618 | 681 | 1,375 |
| VirB8+1 | 1,730 | 20 | 27 |
| VirB9+1 | 1,856 | 354 | 772 |
| VirB10+1 | 1,772 | 446 | 859 |

**57.9% of the target's own nominations were anchors**, dominated by VirB6 (446) --
which IS VirB5+1 whenever no gene is inserted. All 14 families contributed
contaminants, including VirD4 (22), VirB11 (25) and VirB3 (11); checking only the
expected VirB6 would have left ~490 anchor proteins in the pool.

## Slot definition used

**VirB5 +1 in the VirB5 strand frame, all architectures** (`slot_anchor_is_virb5.md`).
The instruction specified "+1 for canonical/IncI2, -1 for pEC4115"; that -1 was
measured relative to VIRB6. The -1 occupant was recorded separately for 6,529
plasmids and is available if the pEC4115 -1 is real in the VirB5 frame.

## Two code defects found, both of which would have produced plausible wrong numbers

1. **Round 1 would have built family HMMs from the medoid alone** while its
   docstring claimed a medoid -> hmmalign -> rebuild bootstrap, because Round 0
   never saved cluster membership. `in-slot/total` is a RATIO; a weak profile
   shrinks numerator and denominator unequally, biasing every family by an amount
   depending on its divergence. Fixed: Round 0 emits `member_seqs`.
2. The seed-control resolution (above).

## The negative slots are conservative but not independent

VirB8+1, VirB9+1 and VirB10+1 share the operon with the target. They control for
"adjacent to a T4SS gene" but NOT for "small lipoprotein anywhere on a plasmid".
Since lipobox is a CONTENT criterion rather than a positional one, a same-operon
negative is a weaker test for it than for a positional criterion.

The genuinely independent baseline -- an adjacent gene pair from a non-T4SS module
(Rep/Par/MOB) -- requires new HMMs and STILL DOES NOT EXIST. The 0/0/0 at
lipobox>90% is real, but it is not the strongest available test of this criterion.

## Standing output

**Three defensible novel families: cl39, cl46, cl13.** Plus cl34 set aside as an
extreme outlier, and cl25 recovered as the known Eex_IncN.

Reported as families, not as record counts: the criterion that selected them
(lipobox) is validated against four seeds under three definitions at 100%, passes
the positive control, and gives zero at all three negative slots -- but the
negative slots are same-operon, and no family here has been shown to DO anything.
These are candidates with a defensible selection procedure, not confirmed
exclusion genes.

---

# Negative slots for MPF_F and MPF_I (jobs 26068616-23)

First negative controls these two classes have ever had. Before this, it was not
possible to test whether ANY criterion separates in MPF_F or MPF_I.

Anchor pairs chosen from MEASURED backbone adjacency:
MPF_F has only two usable pairs (the sole other pair >=30% is `TraH -> TraG_N`,
the target's own side); MPF_I has six, three used.

| class | slot | role | uniq nom | survival | len med | lipobox | named family |
|---|---|---|---|---|---|---|---|
| MPF_F | TraG_N+1 | **TARGET** | 1209 | **67.8%** | 177 | 4.1% | **TraS 63** |
| MPF_F | TraU+1 | negative | 1026 | 19.1% | 132 | 1.5% | NONE |
| MPF_F | TrbC_Ftype+1 | negative | 1401 | 11.8% | 126 | **6.1%** | NONE |
| MPF_I | traY+1 | **TARGET** | 492 | **92.7%** | 204 | 4.4% | **ExcA 204** |
| MPF_I | traM+1 | negative | 459 | 4.6% | 107 | 0.0% | NONE |
| MPF_I | traP+1 | negative | 227 | 2.2% | 169 | 0.0% | NONE |
| MPF_I | trbA+1 | negative | 530 | 11.7% | 210 | 1.6% | NONE |

**Both slot definitions validated**: named exclusion families appear at the target
slots and at ZERO of the five negatives -- the MPF_T pattern (372 vs 0/0/0)
reproduced in both classes.

MPF_I's slot is the cleanest of the three classes (92.7% survival; MPF_T's own
target loses 57.9% to decontamination).

## `pair_occupancy` -- a column name that must not be compared across classes

MPF_T's "pair occupancy" means **is there a gene BETWEEN the two anchors** (an
insertion rate): target 73.4% vs negatives 0.7 / 15.2 / 20.9%.

The MPF_F/MPF_I run computed something different under the same name: **is there
any CDS at anchor+1**, which is trivially ~100% everywhere (94.2-100.0% on all
seven slots) and carries no information.

**The quantity that carries MPF_T's meaning here is the decontamination survival
rate** -- the fraction of nominations that are NOT the next anchor. Do not compare
`pair_occupancy` between the MPF_T run and this one; they are different
quantities.

## lipobox is INVERTED in MPF_F

MPF_F target **4.1%** vs negative `TrbC_Ftype+1` **6.1%**. Not weak -- pointing the
wrong way. MPF_I is 4.4% vs 0/0/1.6%: right direction, negligible magnitude
against MPF_T's 56.1%.

Fourth measured instance of "criteria do not transfer across classes", and the
strongest form: the criterion giving 20x separation in MPF_T is inverted in MPF_F.

# Compositional criteria tested for MPF_F / MPF_I

| axis | result |
|---|---|
| lipobox | INVERTED in MPF_F (4.1% vs 6.1%); negligible in MPF_I |
| length | **DEAD in both.** MPF_F's apparent signal was an artefact -- see below. MPF_I fails outright (204 vs trbA+1's 210) |
| KD hydrophobicity | previously 2.3x against a 40.5% background in MPF_I -- insufficient alone |
| **TM segment count** | **UNAVAILABLE -- fails its controls** (below) |
| pI, amino-acid composition | not yet tested |

## TM segment count: unavailable, after three gate attempts

No TMHMM / Phobius / TMbed installed (registered downloads). No torch/transformers,
so DeepTMHMM cannot run locally; the BioLib route executes REMOTELY and would send
sequences off-cluster -- not done without explicit approval.

**EMBOSS `tmap` is installed and REJECTED**: returns `HitCount: 0` on the R64 ExcA
N-terminus (which contains an unmistakable TM stretch) AND on a hydrophilic
control. It implements Persson & Argos, which needs an MSA. A predictor returning
zero on both a known positive and a known negative is inert, not conservative.

A local Kyte-Doolittle window counter was then gated three times:

1. **Synthetic negatives** (poly-E/K, charged-random, proline-rich): passed 50 of
   60 settings. That is a bad GATE, not a good counter -- any hydrophobicity
   method separates those.
2. **Real proteins, wrong positive**: VirB6 (polytopic, median 351 aa) vs soluble
   ATPases. Failed everywhere. But the positive was wrong for the task: the
   targets are SMALL membrane-associated proteins (ExcA ~220, TraS ~163, TrbK 69),
   carrying one or two TM segments, not many.
3. **Named exclusion families vs SIZE-MATCHED proteins from the same plasmids**
   (essential: TM count scales with length, so unmatched negatives would separate
   on length alone and the counter would look good for the wrong reason), sweeping
   the count threshold k as well as the window parameters:

   | win | thresh | minlen | k | pos >=k | neg >=k |
   |---|---|---|---|---|---|
   | 15 | 2.0 | 15 | 1 | 83.1% | 22.5% |
   | 15 | 1.6 | 18 | 1 | 83.3% | 24.5% |
   | 19 | 1.2 | 15 | 2 | 73.2% | 18.9% |

   No setting meets the pre-registered bar (positives >=80%, negatives <=20%).

There IS signal -- 83.1% vs 22.5% is 3.7x enrichment -- but it misses the bar, and
a 22.5% false-positive rate on size-matched same-plasmid controls cannot call
individual proteins. Against MPF_T's lipobox (56.1% vs 2.8%, 20x) it is an order
of magnitude weaker. **Recorded as UNAVAILABLE, not as a weak criterion.** The bar
was not relaxed after seeing that it nearly passed.

## MPF_F length: the signal was an artefact of the negative slots (job 26079034)

The negative slots retain only 12-19% of nominations after decontamination, so
their survivors are small leftover ORFs at tightly packed anchor junctions -- a
SELECTED REMAINDER, not a representative background. Re-tested against the control
design used for MPF_I hydrophobicity: proteins from the SAME plasmids.

| set | n | Q1 | median | Q3 |
|---|---|---|---|---|
| ALL proteins on the same 7,803 plasmids | 1,258,260 | 106 | **183** | 293 |
| slot occupant, NAMED (TraS) | 2,534 | 159 | 163 | 165 |
| slot occupant, UNNAMED | 5,269 | 172 | 176 | 189 |
| slot occupant, ALL | 7,803 | 163 | 173 | 182 |

P(slot occupant longer than a random protein from its OWN plasmid):

| set | P | reading |
|---|---|---|
| NAMED (TraS) | **0.458** | slightly SHORTER than background |
| UNNAMED | **0.538** | indistinguishable |

Both ~0.50. The background median (183) is HIGHER than the slot occupants' (173).

**The 177-vs-126/132 comparison did not measure what it appeared to.** Length has
no discriminative power in MPF_F. The disjoint named/unnamed IQRs (159-165 vs
172-189) remain real as a difference BETWEEN the two occupant sets, but neither
set is remarkable against its own plasmids.

# Standing position for MPF_F and MPF_I

**No compositional criterion has been shown to work in either class.**

| axis | MPF_F | MPF_I |
|---|---|---|
| lipobox | INVERTED (4.1% vs negative's 6.1%) | negligible (4.4% vs 0/0/1.6%) |
| length | no signal (P = 0.458 / 0.538) | fails (204 vs 210) |
| KD hydrophobicity | untested | 2.3x on a 40.5% background |
| TM segment count | unavailable (fails controls) | unavailable |
| pI, amino-acid composition | untested | untested |

What the negative slots DID establish is that both slot DEFINITIONS are sound --
named families at the targets, zero at all five negatives. That is the
precondition for a criterion search, not a result of one.

The remaining untried axis with a mechanistic rationale is §13's transcriptional
independence (own promoter / Rho-independent terminator), which is independent of
both sequence and position and unaffected by percolation. It requires new motif
scanning and has not been attempted.
