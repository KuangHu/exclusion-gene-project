# Catalogue v1 -- what it is, and four things it is not

`data/catalogue_MPF_T_v1.tsv`, 7,436 rows. GA is the primary measure throughout.

## Frozen numbers

| | |
|---|---|
| tier A (7/7 core anchors) | **6,717 (90.3%)** |
| tier B (6/7) | 344 (4.6%) |
| tier C (5/7) | 29 (0.4%) |
| tier D (<5/7) | 346 (4.7%) |
| `entry_criterion_only` (<=2/7, admitted on VirB4 alone) | **282, all in tier D** |
| `slot_ready__virb5_virb6` | 5,272 (70.9%); 75.2% of tier A |
| architecture | canonical 6,020 / rearranged 777 / uncallable 639 |

`tier` is skeleton completeness ONLY. Architecture-callability lives in
`slot_ready__*`, one column per slot, because MPF_F has two slots and folding
architecture into tier would bake an MPF_T assumption into a class-agnostic table.

VirB1 is excluded from tiering: of 1,452 plasmids lacking it, 1,118 (77.0%) carry
all seven core anchors, so its absence is real. Gating on it would drop ~14% of
plasmids for lacking an optional peptidoglycan hydrolase.

## 1. VirB7 is permanently unreliable, and six-frame does not fix it

GA detects VirB7 on 206/7,436 (2.8%). phmmer recall recovered exactly ZERO,
because the proteins are not in the protein set at all.

The 206 that ARE callable are median 56 aa (Q1 43, Q3 56) -- the same short range
as pKM101's uncalled 48 aa traN -- and 94% overlap no neighbouring CDS while
pKM101's overlaps by 11 bp. So length is not the discriminator and overlap is,
which predicted that six-frame with an overlap<=50% rule should recover them.

**It did not.** Six-frame finds novel ORFs in 551 of 1,882 VirB6|VirB8 gaps
(29.3%) at median 57 aa, but only **2.3% carry a lipobox** -- indistinguishable
from the ~2.8% background. pKM101 traN DOES carry one (`ISAC@15`). Right length,
right position, wrong architecture: these are chance ORFs.

Consequence: the VirB7 position cannot be scored, and the MPF_F TraN slot inherits
the problem, since `PF06986` is itself only a domain model (median 39.8% target
coverage even on short hits).

## 2. `n_unresolved_gaps` is not yet usable

v1 counts genes between the FIRST and LAST anchor. Recomputing on ADJACENT anchor
pairs did NOT fix it -- median 11 either way, Q3 32 -> 36, max 2,815 -> 4,400.

The cause: adjacency in the anchor ORDERING is not adjacency on the genome. Gap
sizes reach **1,838,810 bp**. A bp ceiling (~10 kb) is required before this column
means anything, and no downstream analysis should consume it until then.

## 3. Three architectures, and the invariant is adjacency not side

| architecture | order | eex offset from VirB5 |
|---|---|---|
| canonical | VirB5 -> [eex] -> VirB6 | +1 |
| IncI2 inverted | VirB6 -> VirB5 -> [eex] | +1 (69%), **+2 (30%)** |
| pEC4115 | [eex] -> VirB5 -> ... | **-1** |

The architecture difference is **VirB6's position**; only pEC4115 flips the eex
side. pEC4115 is 11 records across **9 independent Mash clusters**, so not clonal.
The IncI2 +2 subgroup (151 records) is unresolved and blocks finalising
"eex is adjacent to VirB5".

## 4. Fusions are real but mostly harmless here

`VirB3|VirB4` fusion on **2,359 plasmids (31.7%)**, median 917 aa, which fully
explains VirB3's bimodal length (Q1 99, Q3 916). It does NOT corrupt
`core_completeness` -- VirB3 is not one of the core seven -- and it sits upstream
of VirB5|VirB6, so slot arithmetic is unaffected. `VirB4|VirB4b` at 6,847 is not a
fusion: PF03135 and PF27097 are two models of one protein.

VirD4 is unimodal (600-699 dominant, 500-599 shoulder, thin tail below 500), so
its 20.8%-outside-robust-band is a one-sided tail, not a second subtype.

## Second-pass recall is an ANNOTATION, not a measure

Its precision control is structurally blind where it claims most: the conflict
filter rejects proteins already assigned to another family at GA, which catches
VirB5's 816 aa VirB4 contaminants (72% rejected) but cannot see VirB2's plausible
~120 aa ones (0.3% rejected). Zero rejection is evidence the filter did not fire.
See `docs/ga_threshold_failure.md` for the full retraction.
