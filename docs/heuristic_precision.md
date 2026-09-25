# The §4 discovery heuristic has ~11% precision. Measured, not estimated.

Run `scripts/64_decoy_analysis.py`. Output `data/exports/decoy_analysis.tsv`.

## The question

Inside the ±6 kb window around the VirB6-family target, how many ORFs pass §4
step 3's length filter (40–260 aa) but are **not** the exclusion gene?

## The answer

| element | target | ORFs in window | pass 40–260 aa | **decoys** | true d (bp) | true rank by distance |
|---|---|---|---|---|---|---|
| F | traG | 15 | 11 | **10** | 33 | **1/11** |
| ICEBs1 | conG | 16 | 8 | **7** | 1570 | 3/8 |
| R621a | traY | 11 | 5 | **4** | 76 | 2/5 |
| R64 | traY | 12 | 6 | **5** | 71 | 2/6 |
| RP4 | trbL | 20 | 16 | **15** | 7 | **1/16** |
| SXT | traG | 12 | 7 | **6** | 33 | **1/7** |

**Length filter alone: 6/53 = 11% precision.** Nine decoys per true positive.
**Nearest-ORF by distance: 3/6 rank the true gene first.**

## F, the worked case

The ±6 ORF window around `traG` contains nine ORFs of 40–260 aa. One is `traS`:

```
artA  104aa  d=2905   trbJ  113aa  d=1737   traT  244aa  d= 576
traQ   94aa  d=2540   trbF  126aa  d=1370   trbH  239aa  d=3716
trbB  181aa  d=2008   traS  173aa  d=  33   traX  248aa  d=9722
                            ^ the only true positive
```

`trbH` and `trbF` are in the same operon, on the same strand, inside the window,
and in the length band. **Nothing in the implemented filter excludes them.**

## What this means for the schema

`supports_discovery_heuristic` currently answers only *"is the exclusion gene
inside the window?"* — which 53 other ORFs also satisfy. It is a **recall**
statement and must never be read as precision. The exports carry that warning in
their header.

## What is missing

§4 step 3 specifies four conjuncts. Two are implemented, two are not:

| conjunct | status |
|---|---|
| length 40–260 aa | implemented |
| short intergenic distance to VirB6 | implemented (as the ±6 kb window) |
| **lipobox (LipoP/SignalP) OR 1–2 TM helices (DeepTMHMM)** | **NOT implemented** |
| **no confident InterProScan assignment** | **NOT implemented** |

The two missing conjuncts are the discriminating ones. `trbB`, `trbJ`, `trbF`,
`trbH` all have confident functional assignments — trbB is VirB11-family, trbJ is
VirB5-family — so an InterProScan negative-assignment filter should remove most
of the decoy set. That is the next thing to build, and this table is its
benchmark: **any candidate filter must be scored against these 53 ORFs, not just
against the 6 true positives.**

Distance is the strongest implemented signal but is not sufficient on its own:
ICEBs1 puts the true gene third at 1570 bp, behind two closer decoys.
