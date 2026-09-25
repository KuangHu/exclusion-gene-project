# Hand-verification pass — stage 3 candidates

**Run:** 2026-08-31 · `scripts/30_extract_pairs.py` · 28/41 target genes recovered · negative arm PASS · 0 length-window flags

Every row below is a CANDIDATE awaiting human sign-off. On sign-off the `coords` + `aa_sha256` pair is frozen into the seed table and re-extraction becomes a regression test.

## Acceptance gates

| Gate | Result |
|---|---|
| R100 `traS` recovered | **PASS** — `/gene` in GenBank, `/product` in RefSeq |
| RP4 `trbK` found | **PASS** — 69 aa precursor, `CAJ85697.1` |
| RP4 partner left NULL (no `trbL` in seed) | **PASS** |
| `sfx` restricted to IncC/IncA | **PASS** — pED208 collision not matched |
| RSF1010 zero hits | **PASS** |
| pAM373 zero hits | **PASS** (after the prgA fix) |

## Verified candidates

| element | gene | side | route | accession | aa | protein_id | aa_sha256[:12] |
|---|---|---|---|---|---|---|---|
| F | `traS` | exclusion | gene | `AP001918.1` | 173 | `BAA97970.1` | `0610d8b0fcad` |
| F | `traT` | exclusion | gene | `AP001918.1` | 244 | `BAA97971.1` | `2b93ec7a44ff` |
| F | `traG` | partner | gene | `AP001918.1` | 938 | `BAA97969.1` | `94d11378fb19` |
| ICEBs1 | `yddJ` | exclusion | product | `CP171645.1` | 126 | `YGG23124.1` | `268dccf2ab1e` |
| ICEBs1 | `conG` | partner | gene | `CP171645.1` | 815 | `YGG21655.1` | `ba5d713f59e3` |
| R100 | `traS` | exclusion | gene | `AP000342.1` | 159 | `BAA78881.1` | `d54e540dd46b` |
| R100 | `traT` | exclusion | gene | `AP000342.1` | 243 | `BAA78882.1` | `6a39d3790984` |
| R100 | `traG` | partner | gene | `AP000342.1` | 940 | `BAA78880.1` | `9cf4ff1716fd` |
| R16a | `eexC` | exclusion | product | `NZ_KX156773.1` | 78 | `WP_000718002.1` | `9f7cef892ad4` |
| R16a | `traG_C` | partner | product | `KX156773.1` | 1204 | `AOB42072.1` | `3e4ca91f937f` |
| R16a | `traN` | partner | product | `KX156773.1` | 932 | `AOB41975.1` | `0bdc2dc8ca87` |
| R391 | `traG` | partner | product | `AY090559.1` | 1189 | `AAM07996.1` | `86e988505cff` |
| R621a | `excA` | exclusion | gene | `AP011954.1` | 204 | `BAK64483.1` | `cb6b79a76aa8` |
| R621a | `traY` | partner | gene | `AP011954.1` | 745 | `BAK64484.1` | `9177558bfefb` |
| R64 | `excA` | exclusion | gene | `AP005147.1` | 220 | `BAB91650.1` | `43ab07ed1f6b` |
| R64 | `excB` | exclusion | gene | `AP005147.1` | 147 | `BAB91651.1` | `ee52017e2755` |
| R64 | `traY` | partner | gene | `AP005147.1` | 745 | `BAB91652.1` | `4108d588a23d` |
| RP4 | `trbK` | exclusion | gene | `BN000925.1` | 69 | `CAJ85697.1` | `23dde58dab57` |
| SXT | `eex` | exclusion | product | `KJ817376.1` | 143 | `AKA21199.1` | `4419c74f454d` |
| SXT | `traG` | partner | gene | `KJ817376.1` | 1189 | `AKA21200.1` | `5098b5fd46bb` |
| pAD1 | `sea1` | exclusion | gene | `CP046109.1` | 891 | `QKR98071.1` | `241d2685baae` |
| pCF10 | `prgA` | exclusion | gene | `AY855841.2` | 891 | `AAA65847.1` | `cfff30495211` |
| pED208 | `traS` | exclusion | gene | `CP146035.1` | 186 | `WWQ71496.1` | `978de7c2ab0f` |
| pED208 | `traT` | exclusion | gene | `CP146035.1` | 245 | `WWQ71497.1` | `72ee3706b9a6` |
| pED208 | `traG` | partner | gene | `CP146035.1` | 965 | `WWQ71495.1` | `d2cfb5374deb` |
| pKM101 | `eex` | exclusion | gene | `U09868.1` | 75 | `AAA86454.1` | `fd988b8cd6b2` |
| pVCR94 | `traG_C` | partner | product | `CP033514.1` | 1204 | `AYV08555.1` | `3e4ca91f937f` |
| pVCR94 | `traN` | partner | gene | `CP033514.1` | 936 | `AYV08496.1` | `f4e6946387f5` |

## Cross-record assertions

| assertion | result |
|---|---|
| RP4 TrbK: `BN000925.1` vs `CP152305.1` | **EQUAL** — 69 aa, `23dde58dab57d3c6`. No substitutions accumulated in the 2024 resequencing |
| RP4 TrbL: `BN000925.1` vs `CP152305.1` | **EQUAL** — 528 aa, `b05af4436bbcd68c` |
| RP4 trbK–trbL adjacency | **CONFIRMED** — `CAJ85697.1`/`CAJ85698.1`, consecutive |
| pED208 TraS across 3 records | **EQUAL** — 186 aa, `978de7c2ab0f28db` |
| pED208 TraT across 2, plus `traTp` | **EQUAL** — 245 aa, `72ee3706b9a67647`; proves `traTp` is a synonym rather than assuming it |
| SXT TraG 606–608 | **`P-G-E` = S group.** ICDC-1307 pin correct despite the `EexR1` gene name |
| pKPC_UVA01 trbK-like | **UNRESOLVABLE by SHA** — absent from both records; tie needs the primary paper |

## Misses — manual coordinate curation required

Confirmed absent from the record text, not pattern failures. §4 predicts exactly this: these genes are small, fast-evolving and mostly have no Pfam family.

| element | gene(s) | note |
|---|---|---|
| ColE1 | mbeD | absent from record text |
| R16a | sfx | absent from record text |
| R27 | eexA, eexB | absent from record text |
| R391 | eex | absent from record text |
| R6-5 | traS, traT, traG | `traT` IS in `X52553.1` but only as a `mat_peptide`, not a CDS product — needs coords |
| R621a | excB | absent from record text |
| Ti | trbK-like | absent from record text |
| pKPC_UVA01 | trbK | absent from both duplicate records |
| pLS20 | ses | absent from record text |
| pVCR94 | eexC, sfx | absent from record text |
