# Gene architecture of all eight MPF classes

Read off the published prototype of each class in coordinate order, job 26175259
(`scripts/173_architecture_all8.py`). Profiles are CONJScan's; exclusion families
are scanned LAST and post-hoc and never enter a class definition.

**n = 1 per class.** This is an architecture readout, not a population statistic.
Nothing here may be quoted as "MPF_x operons look like this" without the
population behind it. What it is good for is locating the exclusion gene relative
to the machinery, which no count can show.

`|` marks a break of more than 3 unrecognised genes. `**X**` is an exclusion family.

---

## The eight

**MPF_T — RP4** (74 CDS, 13 recognised)
```
virB11 virB2 virB3 virB4 virB8 virB9 virB10 virB5 **TrbK** virB6 | t4cp2 MOBP1
```
TrbK (69 aa) sits between virB5 and virB6 — **the VirB5+1 slot, confirmed from
architecture rather than from the slot logic that defined it**.

**MPF_F — F** (105 CDS, 16 recognised)
```
traL traE traK traB | traV virB4 traW traU trbC traN traF | traH traG t4cp2 MOBF
```
The readout shows no exclusion family, which is a MODEL GAP, not biology — see the
TraS correction below. The genes between `traG` and `t4cp2` (= traD) are:
```
97 traH 458aa | 98 traG 938aa | 99 traS 173aa  <- TraG+1 | 100 traT 244aa | 101 traD 717aa
```
**traS is at exactly TraG+1**, confirming the MPF_F slot definition, with traT
(OM surface exclusion, category 2) immediately behind it.

**MPF_I — R64** (134 CDS, 21 recognised)
```
MOBP1 t4cp2 trbB trbA **ExcA** **ExcA** traY traW traV traU traT traR traQ traP traO traN traM traL | traK traI | traE
```
TWO adjacent ExcA hits, 220 aa and 147 aa — the overlapping `excA`/`excB` pair with
translational reinitiation, which is why ExcA has no single topology and why its
2 TM was never a prediction failure. Directly adjacent to traY.

**MPF_FA — ICEBs1** (25 CDS, 9 recognised)
```
orf23 tcpA MOBT | orf13 orf17b virB4 orf15 orf14 **DUF4467**
```
DUF4467 = YddJ, 126 aa, immediately after orf14 (cwlT).

**MPF_FATA — pCF10** (57 CDS, 10 recognised)
```
prgB prgC prgF prgHb prgIc virB4 prgK prgL | t4cp2 | MOBP1
```
No entry-exclusion family, as expected: FATA's characterised exclusion is SURFACE
exclusion (PrgA/Sea1), which sits on the category-2 discard list.

**MPF_G — ICEHin1056** (62 CDS, 19 recognised)
```
tfc2 tfc3 tfc5 t4cp2 tfc7 tfc8 tfc10 tfc11 tfc12 tfc13 tfc14 tfc15 virB4 tfc18 tfc19 tfc22 tfc23 tfc24 | MOBH
```
Eighteen contiguous genes — the tightest block of all eight. No exclusion family.

**MPF_B — CTnDOT** (19 CDS, 14 recognised)
```
MOBP1 | traE traF virB4 traH traI traJ traK traL traM traN traO traP traQ
```
Thirteen contiguous. No exclusion family. (This record is the transfer REGION, so
the flanks are absent by construction.)

**MPF_C — pCC7120a** (386 CDS, 13 recognised)
```
MOBF | MOBV | alr7204 alr7205 virB4 alr7207 alr7208 alr7209 alr7210 alr7211 alr7212 t4cp2 MOBP1
```
Eleven contiguous inside a 386-CDS plasmid. No exclusion family.

---

## What is common to all eight

**virB4 is embedded mid-block in every class without exception.** It is the
structural centre, which is also why it is useless as a class discriminator —
`T4SS_virb4` is mandatory in all eight CONJScan definitions.

Where an exclusion gene is known, it sits **at the edge of the machinery block, not
inside it**: TrbK at the virB5/virB6 junction, traS immediately after traG, ExcA
immediately before traY, DUF4467 at the end of the ICEBs1 block. That is the
positional claim the slot work rests on, and it holds in all four cases.

---

## Exclusion gene status, by class

| class | exclusion gene | evidence | located here |
|---|---|---|---|
| MPF_T | **TrbK** (69 aa) | E1, lipoprotein, C-term truncation abolishes | yes, virB5+1 |
| MPF_F | **TraS** (173 aa on F) | E1, blocks transfer post-MPF | yes, traG+1 — but NOT by PF10624 |
| MPF_I | **ExcA/ExcB** (220/147 aa) | E1, donor TraY is the target | yes, traY-adjacent |
| MPF_FA | **DUF4467 / YddJ** (126 aa) | E1, ConG E288K resistance | yes, end of block |
| MPF_FATA | none (surface exclusion only) | PrgA/Sea1 is category 2 | n/a |
| MPF_G | **none known** | — | — |
| MPF_B | **none known** | — | — |
| MPF_C | **none known** | — | — |

---

## CORRECTION: PF10624 is an R100-specific TraS model, not a TraS family model

`config/anchor_candidate_criterion.yaml` lists
`{model: PF10624, family: TraS, evidence: "F TraS entry exclusion"}`.
The model's own DESC is "Plasmid conjugative transfer entry exclusion protein TraS"
(157 aa, GA 27.4). Measured against every seed plasmid carrying an annotated `traS`:

```
R100        traS 159 aa   PF10624 @ GA = 312.7 bits   HIT
F           traS 173 aa   PF10624 @ GA = none, and none at E<=10
R64         traS  62 aa   none
ColIb-P9, R621a, pED208 (x3), pEK204, pCVM29188   none
```

**1 of 14.** And the three are not divergent homologues — `phmmer` of R100 TraS
against F TraS and R64 TraS finds **no detectable similarity even at E<=100**.
Three unrelated proteins share one gene name.

This is the third instance of the same identifier problem already recorded here:
`PF20084` vs `TIGR04359` (two models, one TrbK name, two architectures), and
`NF033891`'s "surface exclusion" DESC inherited from a 1994 title.

**What it changes:**
* PF10624's category-1 entry stands — R100 TraS is a genuine, characterised entry
  exclusion protein. Its RECALL is what was overstated.
* Any figure of the form "TraS in N plasmids" derived from PF10624 means
  "**R100-type** TraS", and should be relabelled.
* MPF_F's slot occupants are NOT covered by a TraS model. The F, R64 and
  ColIb-P9 proteins are unmodelled, which makes them candidates in the pool rather
  than archived knowns — the opposite of how they have been treated.
* The pending "3 TraS members at 417 aa scoring higher on TraT" are spurious
  against a 157-aa model and should be dropped, not reassigned.

**What it does NOT change:** Rule B. Its multi_pass exemplars are convergent and
phylogenetically unrelated by design — Humbert/Burrus established exactly that for
TraS/ExcA/EexS. This result is that principle showing up one level down, inside the
gene name itself.
