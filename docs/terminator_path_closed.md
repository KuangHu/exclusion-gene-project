# Terminator detection around exclusion genes -- CLOSED

Abandoned by decision. Recorded because the negative is measured, not assumed.

## The question

Does an exclusion gene sit in its own transcriptional unit? A unit needs
boundaries, so both flanks matter and they mean different things: an UPSTREAM
terminator blocks read-through from the upstream operon (decoupling); a
DOWNSTREAM one bounds the gene's own transcription (protects neighbours).

## The answer: there is usually no room for one

Measured on all **6,517** MPF_T plasmids with a resolved slot -- coordinates only,
no folding, no thresholds:

| region | n | Q1 | median | Q3 | >=30 bp | overlap <=0 |
|---|---|---|---|---|---|---|
| upstream of slot | 6,499 | 2 | **9** | 13 | **13.9%** | 22.7% |
| downstream of slot | 6,497 | 2 | **9** | 14 | **11.1%** | 20.9% |
| anchor-anchor (reference) | 6,502 | 2 | 5 | 11 | 1.0% | 17.4% |

A Rho-independent terminator needs roughly 30 nt (stem + loop + poly-U).

Counting "no space" as "no terminator":

```
upstream    86.1% have no terminator -- no room    13.9% testable (~900 plasmids)
downstream  88.9% have no terminator -- no room    11.1% testable (~720 plasmids)
```

**This is a statement about ARCHITECTURE, not about detection sensitivity.** It
rests on gene coordinates, not on folding, scoring or any threshold. Over 20% of
slot occupants overlap their neighbour outright.

## The one positive signal, left unexploited

The slot is **14x more spacious** than the operon interior: 13.9% of slot flanks
reach 30 bp against 1.0% of adjacent-anchor junctions (median 9 bp vs 5 bp). So
spacing does loosen at the slot -- it is a boundary, just not one that usually
carries a hairpin. The ~900-plasmid subset with room is a properly powered test
set if anyone returns to this.

## Unexplained

Named and unnamed occupants differ in SHAPE, not amount:

```
eex_occupied  up median 11 bp, 10% >=30 bp
candidate     up median  4 bp, 20% >=30 bp
```

Lower median but twice the fraction with room -- bimodal, not shifted. No
explanation; "candidate" may simply be heterogeneous.

## What the pilot did NOT establish

The first pilot (job 26127660) reported "SEPARATES" and that verdict was wrong:

* **n = 7 upstream, 9 downstream** across 120 plasmids. The +20.7 median gap was
  computed on seven sites.
* **Zero sites reached the 95th percentile in ANY class**, including the positive
  control -- so by the stated "top 5% of this plasmid's spacers" criterion,
  nothing passed.
* Downstream scored 22.5 points BELOW the negative control. The two sides
  straddled the control in opposite directions, which is what noise at n~8 looks
  like.

The verdict rule checked a median difference and never checked n. That is the
same small-n threshold failure as `in-slot/total`, `Asp@+2`, `TM>=3`, and reading
only the pre-registered column -- the sixth instance.

The low n was NOT an indexing bug: the slot gene resolved on 6,517 of 6,517. The
pilot was correctly reflecting that most slots have no scorable spacer.

## Also incomplete when closed

* only one scorer (ViennaRNA hairpin dG + poly-U); **TransTermHP was downloaded
  and working but never wired in**
* the 80 nt folding window was chosen, never validated by a sweep
* upstream and downstream were scored with `side="both"`, so the directional
  distinction the question depends on was never actually computed
* Rho-DEPENDENT termination has no sequence signature this method can detect, so
  even a clean negative would not have excluded termination

## Tools installed and left in place

`ViennaRNA 2.7.2` (pip), `TransTermHP v2.09` at
`/global/scratch/users/kh36969/bin/transterm_hp_v2.09` with `expterm.dat`.

## Data

`data/anchors/slot_flank_gaps.tsv` (6,517 rows) -- the coordinate measurement,
which stands on its own regardless of this path being closed.
`data/anchors/terminator_pilot.tsv` (146 sites) -- the underpowered pilot.
