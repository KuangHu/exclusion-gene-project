# MMseqs2 fails the seed control. Percent identity is the wrong criterion.

Validation run before clustering anything, using the four slot occupants as a
built-in positive/negative control:

* within-family: RP4↔R751 `trbK` (TrbK_RP4) and pKM101↔R388 (Eex_IncN)
* between-family: 21–29%, i.e. background

## Result: MMseqs2 recovers at most 1 of 2 known families

| setting | families recovered |
|---|---|
| `--min-seq-id 0.30 -c 0.5` (proposed) | **0 of 2** — four singletons |
| `--min-seq-id 0.25` / `0.20` | **0 of 2** |
| `-c 0.8` | **0 of 2** |
| `-s 7.5` (max sensitivity) | **1 of 2** — TrbK_RP4 joins; Eex_IncN does not |
| `search --exhaustive-search 1 -s 7.5 -e 1000` | **1 of 2** — same |

At default sensitivity the k-mer prefilter returns **zero** cross-hits among four
proteins known to form two families. These are 69–76 aa; the prefilter has too
little to work with. Removing the prefilter entirely does not rescue
pKM101↔R388.

**The proposed 30% / 50% parameters would have split both known families and
made every occupant a singleton** — which would have looked like "no families
found" rather than a tool failure.

## The deeper problem: percent identity does not separate these families

| pair | relation | % identity | global score | local score |
|---|---|---|---|---|
| pKM101 ↔ R388 | **same family** | 41.2% | **93** | **104** |
| RP4 ↔ R751 | **same family** | 37.8% | **104** | **104** |
| RP4 ↔ R388 | cross | **28.9%** | **11** | 49 |
| R751 ↔ R388 | cross | 16.7% | — | — |

Within-family identity is 37.8–41.2% and cross-family reaches **28.9%** — the
ranges nearly touch, so no `--min-seq-id` cut separates them cleanly.

**Alignment score does separate them, by an order of magnitude:** 93–104 within
family versus 11 cross-family (global, BLOSUM62 −11/−1). The cross-family local
hit is a 45-residue partial; the within-family hits are full-length (80–82).

## Consequence for the design

Do not cluster on `--min-seq-id`. Cluster on **alignment score / E-value**, with
the seed pairs as calibration:

* within-family: global score ≥ 93, near-full-length alignment
* cross-family: global score ≤ 11, partial alignment only

At the expected n (hundreds to a few thousand occupants) an **all-vs-all
Smith-Waterman is entirely feasible** — 439 occupants is 96k pairs, 3,000 is
4.5M — so the k-mer prefilter that causes the failure is not needed at all.

## Why this was catchable

The four seeds are a ready-made control with known ground truth on both sides.
Running them first cost minutes and prevented adopting parameters that would have
reported "no families" from a set that provably contains two.
