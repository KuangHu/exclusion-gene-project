# What the ANI clustering does and does not give: a label correction

**d = 0.007 (99.3% ANI) achieves DEDUPLICATION. It does NOT give statistically
independent observations.** Every downstream statistic must be labelled
accordingly.

## Why percolation settles this

Clustering was expected to do two jobs at once:

1. remove repeated sequencing → needs a very tight threshold (~99.5% ANI)
2. supply independent observations → needs a threshold reflecting evolutionary
   independence, which would be far looser

The percolation result proves **job 2 is unobtainable by ANI clustering**. The
graph becomes connected between d = 0.016 and d = 0.018 (ANI 98.2–98.4%), and by
95% a single cluster holds 53% of the database. Any threshold loose enough to be
biologically meaningful fuses the database into one blob.

So there is no ANI threshold that delivers both. d = 0.007 delivers the first.

## The consequence that must not be lost

**A hypergeometric test at the cluster level still violates independence.** The
27,679 singletons and 6,757 non-singleton clusters share deep phylogenetic
structure — every IncP plasmid has a common ancestor, so they are not an
independent sample. Clustering removed resequencing artefacts, nothing more.

## Three reporting layers

| layer | unit | n | valid for |
|---|---|---|---|
| record | plasmid record | 72,557 | raw counts; always reported |
| cluster (d=0.007) | ANI cluster | 34,436 | **deduplicated occupancy frequency** |
| PTU / Inc group | plasmid taxonomic unit | ~10–30 in MPF_T | **the only layer where cross-family enrichment testing is valid** |

The n at the third layer is small, but it is the honest n.

Note the four-plasmid GO result is *already* a family-level claim, not a
record-level one: two families (`TrbK_RP4`, `Eex_IncN`) occupy the same slot with
between-family identity of 21–25%, i.e. background. That is an n=2 statement
about families and was never a statement about records.

## COPLA is now needed, and this is when

Previously deferred as "classification background for final candidates". The
PTU layer is what makes it necessary. Running COPLA over all 72k is infeasible,
but running it on the **VirB4⁺ subset only** is tractable — and that subset is
exactly the scope needing PTU labels.

## A note for later, not now

Percolation under single linkage is partly a chaining artefact and does not by
itself prove that 98% is biologically wrong. If a looser, biologically meaningful
stratification is ever needed, the fix is to change the **linkage** (centroid /
greedy, or Leiden on the ANI graph), not the threshold. Not required for the
current design.
