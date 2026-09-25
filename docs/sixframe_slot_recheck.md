# Six-frame recheck of the "empty" slots

## Question

pKM101 traN is a real 48 aa gene that Pyrodigal -p anon does not call. Every
"empty slot" is therefore suspect in one direction: an absence the caller
manufactured. Under v2 6.1 the occupancy figure decides whether the whole
presence/absence direction survives, so this had to be measured, not argued.

## Answer

**Real, confirmed, and bounded.** Canonical occupancy moves 82.2% -> 84.2%.

## How the number was reached, including three wrong turns

### Geometry first: 90% of empty slots cannot hide anything

A 40 aa ORF needs 123 bp. Over 1,449 empty slots:

| gap between anchors | n | % |
|---|---|---|
| anchors overlap (<0) | 148 | 10.2% |
| 1-20 bp | 1112 | 76.7% |
| 21-50 bp | 38 | 2.6% |
| 51-122 bp | 5 | 0.3% |
| **123-300 bp** | **141** | 9.7% |
| 301-1000 bp | 5 | 0.3% |

Only **146 (10.1%) are searchable at all**, and 145 of those are canonical --
inverted "empty" slots are abutting anchors, so this correction cannot touch them.
The 51-122 / 123-300 discontinuity (5 rows then 141) is not a random distribution.

### Wrong turn 1 -- searching the whole window

Returned a "novel" 58 aa ORF in 11 of 13 occupied slots. Those windows are
250-280 bp and already hold an 85 aa gene: frame shadows, every time.

### Wrong turn 2 -- restricting to non-coding space

Killed the shadows, but **pKM101 traN overlaps the next called CDS by 11 bp, and
that overlap is why Pyrodigal misses it.** The fix would have excluded the one
gene the exercise exists to find. Replaced by an overlap FRACTION filter (<=50%
of the ORF inside a called CDS): traN is 92% outside and passes, a shadow is ~100%
inside and does not. Genes wholly nested in a called gene stay invisible -- the
same limit v2 8 records for StORF, and ColE1 mbeD inside mbeA is a real instance,
not a theoretical one.

### Wrong turn 3 -- an asymmetric control

The first matched control scanned ~2 kb on the slot side against ~200 bp per
control gap. Caught before the result was accepted; notably the flaw favoured the
conclusion the run was contradicting, not the one it was supporting.

### Controls that gate the run

* POSITIVE: pKM101 traN 7899..8045 recovered exactly (48 aa, ATG) in all 16 shards.
  A control that can fail for an interesting reason.
* BASELINE: six-frame re-finds every Pyrodigal occupant. **4,415 recovered, 0
  missed.** An earlier version missed 249/2,972 because the window ran
  anchor-end+1 to anchor-start-1 and occupants routinely overlap an anchor.

## The discriminator: lipobox, not hit rate

Hit rate alone is weak. Symmetric, size-matched: slot 89.7% vs control 31.7%
(2.83x), but only **1.89x** in the dominant 211-240 bp bin. ORF length barely
separates (median 51 vs 48).

The lipobox separates completely, and length does not explain it:

| ORF length | slot lipobox | control lipobox | ratio |
|---|---|---|---|
| 40-44 | 11.9% | 3.3% | 3.6x |
| 45-49 | 6.2% | 2.9% | 2.2x |
| 50-54 | 78.7% | 2.7% | 29x |
| 55-59 | 65.4% | 2.1% | 31x |
| 60-69 | 73.9% | 1.7% | **44x** |
| 70-84 | 89.5% | 3.1% | 29x |

Control rate reweighted to the slot length distribution: **2.4% against 56.1%**.
The signal lives at >=50 aa; the 40-49 aa bins look like chance ORFs, which is
what they should look like.

Defining a credible occupant as **>=50 aa with a lipobox**:

* slot gaps: 99/146 = 67.8%
* matched control gaps: 37/8,142 = **0.5%** false-positive rate
* chance-corrected: 98 of 146

## Direct confirmation

Of six recurrent slot proteins, one is an already-named exclusion gene the caller
failed to emit:

| orf | aa | records / clusters | family | phmmer vs controls |
|---|---|---|---|---|
| orf4 | 59 | 2 / 1 | **NF033894 Eex_IncN, bit 59.8 (GA 33)** | R388 cand, E 1.6e-05 |
| orf1 | 57 | 11 / 2 | TIGR04359 bit 24.1 (GA 23) -- marginal | **no hit** to either trbK |

Within-family reference on the same aligner: RP4/R751 trbK E 1.3e-09, pKM101 eex /
R388 E 1.1e-11. orf1's HMM hit clears GA by 1.1 bits and is not corroborated by
sequence comparison. **One confirmed, one pending.** The conclusion rests on orf4.

All six carry a lipobox at the canonical position (LVAC@22, ITAC@17, VTGC@18 x3,
LAAC@4). An earlier "no lipobox" reading was a regex anchored on `^M`; a six-frame
ORF opened at GTG/TTG translates to V or L.

Clonality matters: the top recurrent protein appears in 47 records but **1 Mash
cluster** -- one clone, not 47 observations. The genuinely spread family is the
72 aa VTGC protein, three near-identical variants across **9 clusters**.

## Effect

Report three numbers with their labels, never one. The correction applies to a
single, named definition and to one subpopulation within it.

| figure | value | what it is |
|---|---|---|
| `slot_occupancy \| VirB5 detectable` | **82.2%** (4031/4905) | as measured, canonical only |
| `+ six-frame recovery` | **84.2%** (4129/4905) | lipobox-based, chance-corrected |
| coverage of the correction | **145 of 146** searchable slots are canonical | see below |

**The six-frame correction does not cover the inverted class, and the numbers
make that easy to misread.** Admitting inverted plasmids raised the empty count
from 874 to 1,449, but the searchable count stayed at 146 -- because the 575
additional inverted empty slots are almost all abutting anchors, with no room for
a 40 aa ORF. Only **1 of 575** inverted empty slots is searchable at all. So:

* the correction is applied to `slot_occupancy | VirB5 detectable`, canonical
* it is NOT applied to `neighbour_is_unassigned | VirB6 detectable`
* it says nothing about the inverted class, which remains unmeasured in this
  respect rather than measured-and-clean

Adding a numerator obtained by one method to a denominator defined by another is
exactly the error to avoid here; the label carries the conditioning.

82.2% does not collapse. The manufactured absence is real and worth ~2 points on
the canonical, VirB5-conditional figure.

## Scale is not established

One case (orf4) proves the phenomenon exists. It does not size it. The ~2 points
comes from 99 credible occupants, of which exactly **one** is confirmed by a
family model; the other 98 rest on the lipobox criterion, which is a strong
discriminator (0.5% false-positive rate against 8,142 matched gaps) but not a
family assignment. Treat 84.2% as an upper-ish estimate of a real effect whose
magnitude is supported by one confirmed instance.
