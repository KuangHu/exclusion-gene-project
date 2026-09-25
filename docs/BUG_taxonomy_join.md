# BUG: the PLSDB taxonomy join returns nothing

## Symptom

`scripts/163_asp2_acceptance.py` stratified 563 lipobox-positive named sequences
by host taxon and put **all 563 into "other / unknown"**. The Enterobacteriaceae
bucket was empty, which cannot be right -- PLSDB is dominated by clinical
Enterobacteriaceae.

## Cause (not yet fixed)

The join reads `taxonomy.csv` and looks for a genus under one of:

```python
g = r.get("genus") or r.get("GENUS") or ""
a = r.get("NUCCORE_ACC") or r.get("accession") or ""
```

One or both column names are wrong for this file. The lookup silently returns ""
for every accession, and `"" in ENTERO` is False, so every sequence falls to the
"other" bucket. **No error is raised** -- the stratification reports a clean,
complete-looking table that is entirely one bucket.

That is the same silent-failure shape as the anchors leaking into the candidate
pool: a step that does nothing still produces a plausible number.

## Why it did not matter this time

Asp@+2 was rejected on the overall rate: Asp is the 6th most common residue at
mature +2 at **8.0%**, and the two known true positives pass at 45.7% (TrbK) and
**6.5%** (Eex_IncN). No host stratification could rescue a rule with no
enrichment, so the conclusion stands regardless.

## Why it must be fixed before the next use

The join is required by three pending lines of work:

* host distribution of each exclusion family
* the Gram-positive ICE line (host phylum is the selector)
* any statement about whether a rule is species-restricted

A silent empty join in any of those would produce a distribution figure that is
wrong in a way nobody would notice from the output.

## Fix (cause confirmed)

Two errors, both in the same three lines.

**1. Wrong column names.** The real columns are `TAXONOMY_`-prefixed:

```
TAXONOMY_UID, TAXONOMY_taxon_rank, TAXONOMY_taxon_name, TAXONOMY_superkingdom,
TAXONOMY_phylum, TAXONOMY_class, TAXONOMY_order, TAXONOMY_family,
TAXONOMY_genus, TAXONOMY_species, ...
```

`r.get("genus")` and `r.get("GENUS")` both miss `TAXONOMY_genus`.

**2. There is no accession column in `taxonomy.csv` at all.** It is keyed by
`TAXONOMY_UID`. The accession lives in `nuccore.csv`, which carries both
`NUCCORE_ACC` and `TAXONOMY_UID`. So the join is two-step and was written as
one-step:

```
nuccore.csv:  NUCCORE_ACC -> TAXONOMY_UID
taxonomy.csv: TAXONOMY_UID -> TAXONOMY_genus / TAXONOMY_family / TAXONOMY_phylum
```

`r.get("NUCCORE_ACC")` against `taxonomy.csv` returns None for every row, so the
dictionary was empty and every lookup fell through to "".

## Guard to add

Binding the right names is not enough -- the failure was silent, and that is the
part worth fixing. Any taxon join must assert that a large fraction of
accessions resolve, and raise otherwise. A join that resolves 0% currently
produces a clean-looking table with one bucket holding everything.
