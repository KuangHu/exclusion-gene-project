# Rules A, B, C -- what was done, what was tested, what was not

One criterion. A/B/C are values of the `anchor_mode` FIELD, not three gates.
Frozen in `config/anchor_candidate_criterion.yaml` v1.0 (sha `a90fea17`).

```
ANCHOR_CANDIDATE = in slot
                 AND length <= 250 aa
                 AND (lipobox_strict OR tm_count >= 1)
                 AND not hit by category 2 / category 4
```

Everything else -- TM count, orientation, post-anchor length, the other two
lipobox definitions -- is RECORDED, never used for admission.

---

# The three modes, side by side

| | **A** lipid | **B** multi_pass | **C** single_pass |
|---|---|---|---|
| rule | lipobox ∧ TM=0 | ¬lipobox ∧ TM≥2 | ¬lipobox ∧ TM=1 |
| exemplars | TrbK, YddJ, Eex_IncN | TraS, EexS, EexR, ExcA | **none** |
| evidence tier | **E1** | **E1** | **E4** |
| MPF_T slot share | 45.8% | 1.0% | 23.8% |
| independent positives | 3 | **2** (EexS/EexR are one 143 aa pair) | 0 |
| negative-slot test | passed (0 hits) | **failed at MPF_T** | passed (6-81x) |

No anchoring evidence at all: 28.5% of the MPF_T slot. Ambiguous (lipobox ∧
TM≥1): 0.9%.

---

# RULE A -- lipid-anchored

## Tested and passed

**Literature panel, 3/3.** TrbK (69 aa), YddJ (126), Eex_pKM101 (75) are
lipobox-positive under all three definitions. TMbed independently gives all three
a signal peptide and ZERO TM -- two methods agreeing from different evidence.

**Negative slots, zero hits.** Named families appear at the target slot and at
none of the three MPF_T negative slots.

**Per-family rates on the named set** (`class_A_named.tsv`): Eex_IncN 89.9%
strict, TrbK 87.8%, DUF4467 55.6%.

## Tested and FAILED -- and the failure is the important part

**Asp@+2 (Lol avoidance) does not scale.** On the panel it looked mechanistic:
TrbK and Eex both Asp, TraT Gly, 3/3 in the right direction. At n=563:

```
Y 21.1%  S 17.1%  K 12.3%  G 10.8%  E 9.6%  D 8.0%  F 7.6%
```

Asp is the **6th** most common residue at mature +2. The two known true positives
pass at **45.7%** (TrbK) and **6.5%** (Eex_IncN). A hard gate would have
discarded 93.5% of the largest named family. **Rejected; not in the criterion.**

**`structural` lipobox is unreliable.** It fires on TraS at 90.1% where strict
gives 6.8%, and TraS is curated inner-membrane. It has a principled basis (TraT's
h-region is 3/8 hydrophobic vs TrbK's 6/8) but n=1 on the decisive case.
Recorded, not promoted.

## The known hole

**TraT passes Rule A.** It is lipobox-positive, and it is an OUTER-membrane
lipoprotein -- surface exclusion, not entry. It is blocked by `PF05818` on the
category-2 list, i.e. **by name, not by the rule**. Residual risk: an UNNAMED
outer-membrane lipoprotein passes A and nothing downstream in this module
catches it.

Also note TMbed gives TraT **1 TM** where cryo-EM shows **3** (alpha-helical OM
insertion is a rare topology the tool is not trained on). So the 1-vs-0 TM cell
cannot be used to separate A from OM lipoproteins either.

## NOT tested

**The 842-sequence recall acceptance was never run.** Criterion: every E1/E2
protein must be admitted; one failure stops scale-up. **This is the outstanding
blocker before Rule A scales.**

---

# RULE B -- multi-pass

## Tested and passed

**All four exemplars behave.** TMbed: TraS 4 TM, EexS 4, EexR 4, ExcA 2. EexS's
4 TM matches the published TMHMM/SPLIT prediction exactly -- the stop condition
("if TMbed cannot call >=3 TM on small IM proteins, halt") was not triggered.

**Stable on the close pair.** EexS and EexR give near-identical segment
boundaries, so the tool is not unstable across near-neighbours.

**Convergent, not homologous.** TraS, ExcA and EexS are phylogenetically
unrelated (Humbert/Burrus) yet all inner-membrane. `anchor_mode` classifies a
CONVERGENT architecture -- which is exactly why sequence homology carries no
cross-class signal and the structural layer is mandatory.

## Tested and FAILED -- but the test was mis-specified

**Negative-slot acceptance at MPF_T: 0.5% target vs 1.3-2.0% negatives,
enrichment 0.3-0.4x.** Pre-registered rule says Class B does not stand.

**That test was in the wrong class.** All B exemplars live in MPF_F (TraS) and
SXT ICEs (EexS/EexR). It was measured against the MPF_T slot, whose known
occupants are Class A lipoproteins. The valid conclusion is narrow: **Class B
does not enrich at the MPF_T slot** -- which is what the lipobox result already
said.

## Corrections made

**The TM>=3 gate was dropped.** It came from n=2 independent positives that both
measured exactly 4. ExcA's 2 TM would have been killed by it. TM is now an
integer field, not a gate.

**ExcA's 2 TM is not a failure.** It has DUAL LOCALISATION -- inner-membrane and
cytoplasmic soluble forms, `excAB` overlapping with translational reinitiation. A
protein with two states has no single topology. Same class as the CrcB
calibration failure.

**ExcA is entry, not surface, exclusion.** NCBI's `NF033891` DESC says "surface",
inherited from the Furuya & Komano 1994 title; the 2013 TraY-swap experiment
settles it as entry. **The MPF_I classification and its 2,112 slot occupants
stand.** Seventh identifier-confusion case.

## NOT tested

**The MPF_F slot (TraG_N+1), where B's only independent positive lives.**
Requires negative slots drawn from WITHIN MPF_F -- using MPF_T negatives would
measure differences between MPF classes, not between slots.

**Prediction to archive before running it:** TraG_N+1 should enrich for
multi_pass. If it does, that is bidirectional evidence that anchoring mode
partitions by MPF class.

---

# RULE C -- single-pass

## What it is

**162 sequences, 285 plasmids, 213 Mash clusters.** Median 67 aa, TM starting at
residue 5, 42 aa post-anchor domain. **100% conform** to both architecture
conditions (TM within N-terminal 40 aa, >=30 aa after).

Composition: **73 DUF2749** (PF10907), **30 PF20084**, 59 with no Pfam hit.

## Tested and passed

**Artifact test -- not signal-peptide misreading.** Only 5 of 145 TM=1 target
occupants are lipobox-positive. The obvious confound (a lipoprotein signal
peptide called as a TM helix) is excluded.

**Enrichment.** 24.2% at the MPF_T slot vs 4.0% / 0.3% / 0.0% at the three
negative slots -- **6x to 81x**. The only significant signal in that table.

**Clonality (4b).** 213 Mash clusters over 285 plasmids -- not a clonal
expansion.

**Full set (4c).** Rerun on all 681, not the 600 sample.

## The near-miss

**Blocking rule 4a said "kill it".** 45.1% of C hits DUF2749, above the 30%
"shares one family -> not novel" threshold, so the rule would have added DUF2749
to the exclusion list and discarded C.

**That rule was wrong.** It conflated MODELLED with CHARACTERISED. DUF2749 is
"Protein of unknown function"; CDD notes it appears to derive from the Trb
operon, function unknown. Discarding it would have thrown away the strongest
novel-family candidate -- and by the same logic would have discarded `PF14729`,
itself a DUF and the MPF_FA positive control.

**Corrected:** rerun only if the hit family is EXPERIMENTALLY CHARACTERISED. A
third category was added -- `named_uncharacterized` -- which stays in the pool.

## PF20084 resolved

Zero overlap with `TIGR04359` (0/49, 0/35), and now explained: **TIGR04359 finds
lipid-anchored TrbK (87.8% lipobox), PF20084 finds single-pass TrbK (11.4%)**.
Two models, one name, two architectures. Never merge.

But `PF20084`'s DESC is "conjugative transfer **REGION** protein" -- POSITIONAL,
no functional claim. So it is **not** a sentinel, and C stays at **E4**.

## NOT tested -- C has no orthogonal evidence at all

Every A and B conclusion rests on characterised proteins. C rests on 162
predictions. Promotion path, by cost:

1. **structural** -- 3Di/Foldseek the 162 against TrbK/Eex folds (pipeline exists)
2. **partner pairing** -- VirB6/TrbL allelic diversity on C-carrying vs
   non-C-carrying plasmids
3. **complex prediction** -- C representatives x same-plasmid TrbL, with a known
   pair as positive control

**Independent corroboration already in hand:** an unrelated annotation team
recorded DUF2749 as "appears to derive from the Trb operon" with no knowledge of
this framework -- an outside observer noting the same positional association.

---

# Cross-cutting: what the tooling can and cannot do

**Orientation is unusable.** The i/o calibration FAILED: CrcB/FluC is a
dual-topology protein assembling as an antiparallel homodimer, so it has no fixed
orientation and cannot anchor an `o` label. TMbed also produced a topologically
impossible string on it (H,H with no intervening h). RstR passed (0 TM), so the
zero-TM side is calibrated. **Orientation is recorded, never read** -- which is
why demoting it from the Class B criterion mattered.

**Non-alternating orientation strings are now a free QC flag**, from that failure.

---

# Method failures this module, all the same shape

| # | what | shape |
|---|---|---|
| 1 | `in-slot/total` | threshold rejected a known family (Eex_IncN at 0.095) |
| 2 | `Asp@+2` | 3/3 on a panel, collapsed at n=563 (8.0% background) |
| 3 | `TM >= 3` | gate set from n=2 independent positives |
| 4 | "has a Pfam name -> not novel" | treated MODELLED as CHARACTERISED |
| 5 | read only the pre-registered column | TM=1 enriched 6-81x, in the data, unprinted |
| 6 | terminator pilot "SEPARATES" | verdict on n=7, zero sites above the 95th percentile |

**Every one was a plausible-looking tightening whose cost was true positives.
None crashed. Each produced a number that would have survived review.**

Principle: **sentinels verify RECALL and never set thresholds; thresholds come
only from target-vs-negative distribution comparison, and the FULL distribution
must be printed, not the tested cell.**

---

# Outstanding

| item | blocks |
|---|---|
| 842-sequence recall acceptance (all E1/E2 must pass) | scaling A and B |
| MPF_F slot test with MPF_F-internal negatives | Rule B validity |
| C's orthogonal evidence (structural first) | C leaving E4 |
| category-4 models (Namtar/Attar/Ishtar/AbjA) | unnamed immunity proteins fall through |
| taxonomy join fix (needs `nuccore.csv` bridge) | any host-distribution claim |
