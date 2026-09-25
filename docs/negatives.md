# Negative controls — OPEN GAP

The entry-only scope cut (`docs/scope.md`) removed both negative controls, and
they have **not been replaced**. Until they are, the extractor has **no
false-positive arm**, and precision is unmeasured.

## What was lost and why

| former negative | why it was a negative | why it no longer works |
|---|---|---|
| **pAM373** | *sea1*-less pheromone plasmid — the comparator for pAD1 *sea1* | *sea1* is SURFACE exclusion; pAD1 is out of scope, so pAM373 has nothing to be negative against |
| **RSF1010** | IncQ mobilizable, no self-exclusion | still true, but it carries no T4SS at all, so it never exercised the VirB6-adjacency logic that entry exclusion turns on |

Both caught real bugs while they were live — pAM373 caught the
`aggregation substance` false positive that broke `prgA`. Losing that arm is a
genuine regression in test coverage, not a simplification.

## What a replacement needs to be

A conjugative element with a **complete, functional T4SS including a VirB6
homolog** but **no entry exclusion gene**. That is the only shape that exercises
the actual failure mode: the neighbourhood scan reaching a VirB6 gene and
inventing an exclusion factor beside it.

Candidates from the reference protocol's own negative list, all systematically
exclusion-free and all VirB6-bearing:

- **Ti VirB/D4** (*Agrobacterium*) — but NOT clean: the octopine *trb* region is
  IncP-like and carries an untested *trbK* homolog. Already filed as a candidate,
  not a negative.
- ***Helicobacter* Cag** — VirB6-bearing, no exclusion reported.
- ***Legionella* Dot/Icm** — MPF_I, so it also covers the non-VirB6 axis that
  R64/R621a occupy.

Cag and Dot/Icm are the two to pin. Neither is compromised the way Ti is.

## Until then

Any recall figure quoted from `30_extract_pairs.py` is **recall only**. Do not
quote precision, and do not treat a clean run as evidence the patterns are tight
— three false positives (`surface exclusion`, `aggregation substance`,
`entry exclusion`) were each found by inspection or by a negative, never by
recall.
