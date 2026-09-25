# MPF_FA precheck: what CONJScan actually provides for Gram-positive

Run BEFORE building anything, because the MPF_I lesson was that a Pfam-only
readout produced a conclusion that had to be retracted ("MPF_I breaks the VirB4
convergence" -- false; CONJScan detects TraU, TraY and TrbC).

## There is no `T4SS_FA_*` component profile

CONJScan splits Gram-positive into TWO classes, with class-specific profiles under
prefixes that are NOT `T4SS_FA_`:

| class | n profiles | names |
|---|---|---|
| FA | 7 | `FA_orf13 orf14 orf15 orf17a orf17b orf19 orf23` |
| FATA | 26 | `FATA_prgB C F Ha Hb Ia Ib Ic K L` (pCF10, *Enterococcus*); `FATA_trsC D F G J L` (pGO1, *Staphylococcus*); `FATA_cd411 cd419_1 cd419a cd419b cd424` (*C. difficile*); `FATA_gbs1346 1347 1350 1354 1365 1369` (*S. agalactiae*) |

Searching for `T4SS_FA_*` by analogy with `T4SS_T_*` and `T4SS_I_*` returns
nothing. The naming convention does not hold across classes.

## Every class-specific profile is `accessory`

From `definitions/Plasmids/T4SS_type{FA,FATA,T}.xml`:

| class | mandatory genes |
|---|---|
| FA | `T4SS_virb4`, `T4SS_tcpA`, `T4SS_MOBB` |
| FATA | `T4SS_virb4`, `T4SS_t4cp1`, `T4SS_MOBB` |
| T | `T4SS_virb4`, `T4SS_t4cp1`, `T4SS_MOBB` |

FA declares 3 mandatory / 7 accessory; FATA 3 mandatory / 27 accessory.

### Consequence 1: MPF_T's tiering does NOT transfer

MPF_T tiers on `CORE = VirD4, VirB4, VirB6, VirB8, VirB9, VirB10, VirB11` -- seven
structural components treated as expected-present. **CONJScan treats every
Gram-positive structural component as optional.** A `core_completeness_ga` column
built the MPF_T way would count profiles that the reference definition itself does
not require, and its tiers would not mean what MPF_T's mean.

This is the fourth measured instance of "criteria do not transfer across classes",
and the first one caught BEFORE building rather than after:

| instance | measurement |
|---|---|
| MPF_T lipobox 56.1% vs 2.8% background | 20x |
| MPF_F lipobox | 1.8% |
| MPF_I hydrophobicity 2.3x, background 40.5% | weak |
| **MPF_FA tiering** | **no core exists to tier on** |

### Consequence 2: no CONJScan slot anchor for ICEBs1

MPF_I's slot was anchored on `T4SS_I_traY`, a dedicated profile. ICEBs1's `conE`
(the "VirB4-like ATPase"), `conG` and `conJ` have no dedicated profile in either
FA or FATA. `conE` should be reachable via the mandatory `T4SS_virb4`, but the
slot anchor (`conG+3`) has no direct equivalent and will need to come from the
seed work, not from CONJScan.

## The VirB4 convergence is stronger than recorded

`T4SS_virb4` is mandatory in EVERY CONJScan plasmid class definition, not only the
three noted in `mpf_f_entry.md` (MPF_T VirB4, MPF_F TraC, MPF_FA ConE). The entry
criterion choice is consistent with the reference tool's own mandatory set across
all classes.

## Standing gap this makes sharper

The third mandatory gene in every class is a RELAXASE (`T4SS_MOBB`, with the other
`T4SS_MOB*` as exchangeables). **Neither catalogue has a relaxase anchor.** It is
mandatory in every CONJScan class definition and absent from ours. Recorded before
as an open item; MPF_FA is where it stops being optional, because with only three
mandatory genes a missing relaxase removes a third of the definition.

## Which class is ICEBs1?

NOT determined here. ICEBs1 is *Bacillus subtilis*; the FATA profiles are
Enterococcus/Staphylococcus/Clostridioides/Streptococcus and the FA profiles are
`orf*` of unstated origin. Assigning ICEBs1 to FA or FATA is a measurement
(scan the seed against both profile sets), not a guess, and has not been made.
