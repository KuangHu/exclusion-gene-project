# pKPC_UVA01 `trbK-like` — two competing E4 candidates, unresolved

`kamruzzaman2022` (PMID 35007138, PMC8923210) reports the gene but **no aa length
and no coordinates**, so it must be computed from the accession. Doing that on
`CP009465.1` yields two candidates that satisfy *different* halves of the RP4
model, and nothing in the record decides between them.

`trbL-like` (VirB6) is confirmed: `AJE44582.1`, 324 aa, `CP009465.1:21534..22508`,
annotated "conjugal transfer protein TrbL".

## Candidate A — positional

`AJE44581.1` · 152 aa · `21029..21487` · **47 bp upstream of trbL-like** · "hypothetical protein"

Fits the curated operon model `trbJ-like–trbK-like–trbL-like` and mirrors RP4,
where `trbK` sits 6 bp from `trbL`. **But**: no lipobox. Its N-terminus reads
`...SLSANA|AD`, with no cysteine near the cleavage site; its cysteines are at
49, 57, 78, 128. RP4 TrbK is a lipoprotein with `SLAG|C23`.

## Candidate B — architectural

`AJE44587.1` · 63 aa · **3,476 bp from trbL-like** · "hypothetical protein"

```
MKKIILILTTVALISG|C STTMKGGSGQLFELSERHVIYRDLTVPEGLGNEPQRQNVTAPGIAP
└─ signal 1-16 ──┘ └────────── mature 17-63, 47 aa ──────────────┘
```

The **only** ORF in the ±6 kb window carrying a canonical lipobox. Shares three
architectural features with RP4 TrbK: `MKK` start, lipobox cysteine, and a mature
form of **exactly 47 aa**. **But**: it is not adjacent to `trbL-like`, so it does
not fit the curated operon order.

## Why sequence cannot break the tie

Aligned against RP4 TrbK, candidate B is **17.3% identical over the precursor and
13.2% over the mature peptide** — background for peptides this short. §4 predicts
exactly this ("small, fast-evolving, most have no Pfam family"), so low identity
neither supports nor refutes. The shared 47-aa mature length is suggestive but,
at 13% identity, could be coincidence.

## Resolution required

Neither candidate may enter the seed table. Both are E4. The deciding evidence is
in `kamruzzaman2022` itself — its figures or supplementary data should give the
ORF or its coordinates. Until then `pKPC_UVA01.exclusion.sequence` stays
`NOT_RECOVERED`.

If candidate B is confirmed, the curated operon model `trbJ-like–trbK-like–trbL-like`
is **wrong** for this element and `context.adjacency` must be re-derived — its
`supports_discovery_heuristic: true` is currently asserted by homology to RP4,
not measured.
