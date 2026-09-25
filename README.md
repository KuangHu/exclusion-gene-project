# Entry-exclusion genes in bacterial conjugative elements

Curated catalogues and a discovery pipeline for **entry-exclusion (EEx)** genes
across the eight MPF classes of conjugative plasmids and ICEs.

Entry exclusion is the mechanism by which a plasmid-carrying cell resists a
second transfer of the same plasmid. Only a handful of EEx genes are
characterised — TrbK (IncP), TraS (F/R100), ExcA (IncI), YddJ (ICEBs1) — and they
are not homologous to each other. Four of the eight MPF classes have no known EEx
gene in any organism.

Run on PLSDB 2024_05_31_v2 (72,556 plasmids) and an ICEberg-derived set
(9,602 elements), on the Lawrencium cluster (SLURM).

## Status

See **[`docs/STATUS_2026-09-24.md`](docs/STATUS_2026-09-24.md)** for the current
state: which classes are frozen, what each method produced, and what failed.

| class | plasmids | ICEs | known EEx | detected in |
|---|---|---|---|---|
| MPF_T | 7,436 | 1,173 | TrbK | 59.0% |
| MPF_F | 11,120 | 238 | TraS | 32.6% |
| MPF_I | 3,732 | 15 | ExcA/ExcB | 58.6% |
| MPF_FA | 495 | 601 | DUF4467/YddJ | 2.0% |
| MPF_FATA | 1,984 | 514 | none (surface exclusion only) | — |
| MPF_G | 109 | 1,147 | **none known** | — |
| MPF_B | 148 | 285 | **none known** | — |
| MPF_C | 117 | 5 | **none known** | — |

## Layout

```
scripts/        numbered analysis steps; each docstring states what it tests
  lib/          frozen ORF caller, assertions (A1-A14), anchor decontamination
slurm/          one job script per analysis step
config/         frozen criteria, pinned accessions -- CURATED, never machine-written
docs/           per-module write-ups, including the methods that failed
data/release/   frozen catalogues
data/anchors/   measurement outputs
data/seed/      pinned reference records (RP4, F, R100, R64, ICEBs1, ...)
```

Not in git, fetched or regenerated instead: `data/CONJScan/` (clone of
[macsy-models/CONJScan](https://github.com/macsy-models/CONJScan)), `data/hmm/`
(Pfam / NCBIfam / InterPro), `logs/`, and the CDS cache on scratch.

## Conventions this project enforces

These exist because each was learned from a failure that produced a plausible but
wrong number.

- **Sentinels verify recall; they never set thresholds.** Thresholds come only
  from target-vs-negative distribution comparison, with the full distribution
  printed rather than the tested cell.
- **Every benchmark carries a baseline expected to score perfectly.** If it does
  not, stop before reading any other result.
- **No E-value thresholds on anchors** (assertion A13). `E = P * N` is
  database-size dependent; use bitscore/GA cutoffs.
- **Exclusion-gene family models never enter anchors or length priors** —
  post-hoc classification only.
- **Report record / cluster / unique-sequence counts**, never records alone.
- **Gene names are not join keys.** Use `(element_id, gene_name)`; note that
  PLSDB is keyed by RefSeq accession, not GenBank.

## Known identifier traps

- `PF10624` ("TraS") detects **R100-type only** — 1 of 14 `traS`-carrying seed
  plasmids. F, R64 and R100 TraS share no detectable similarity (`phmmer` E>100).
- `PF20084` and `TIGR04359` both mean "TrbK" and have **zero overlap**: one finds
  the lipid-anchored form, the other the single-pass form.
- `NF033891` (ExcA) is described as "surface exclusion" in NCBI, inherited from a
  1994 title; the 2013 TraY-swap experiment establishes it as entry exclusion.

## Reference

Guglielmini et al. (2014) *Key components of the eight classes of type IV
secretion systems involved in bacterial conjugation or protein secretion.*
Nucleic Acids Research 42(9):5715-5727. Source of the eight-class scheme and the
CONJScan profiles used throughout.
