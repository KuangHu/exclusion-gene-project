# State of the four conjugation classes

Where each class stands: what is solid, what is provisional, what is absent.

---

# Summary table

| | MPF_T | MPF_F | MPF_I | MPF_FA |
|---|---|---|---|---|
| plasmids catalogued | 7,436 | 11,120 | 3,732 | 2,391 |
| entry criterion | `PF03135` VirB4 | `PF11130` TraC_F_IV | `T4SS_I_traU` | virb4 conjunction |
| seeds | 4/4 | 5/5 | 5/5 | 5/5 (seed files) |
| slot | **VirB5+1** | **TraG_N+1** | **TraY+1** | **none** |
| slot_ready | 6,548 (88.1%) | 7,803 (70.2%) | 3,607 (96.7%) | -- |
| negative slots | 3 | 2 | 3 | -- |
| named eex at slot | 3,862 | 2,541 | 2,112 | 11 (plasmid-level) |
| content criterion | **lipobox** (56.1% vs 2.8%) | **none** (inverted) | **none** (2.3x) | -- |
| ICEs of this class | 919 | 1,323 | **10** | 838 |

---

# MPF_T (IncP / trb) -- the most complete

**Solid.** Entry `PF03135`, all four controls admitted at GA. Slot VirB5+1,
validated on three negative slots (VirB8+1, VirB9+1, VirB10+1) with 372 named-family
hits at the target and **0/0/0** at the negatives. Three architectures resolved
(canonical / IncI2 / pEC4115); the invariant is eex adjacent to **VirB5**, not
VirB6 -- VirB6 is what relocates.

**The only class with a working content criterion.** Lipobox: 56.1% in-slot vs
2.8% background, 20x, holding under three definitions and 100% on four seeds.
Five families at >90% lipobox, zero at all three negative slots, positive control
(Eex_IncN) recovered.

**Corrected during this work.** The VirB5 `E<=1e-5` threshold was withdrawn -- its
own controls failed it (RP4 TrbJ E=0.22, R751 E=0.0093 under the pipeline's search;
the rescuing E-values came from a protein-vs-Pfam scan). Replaced by
`PF07996 @ GA OR T4SS_T_virB5 @ GA`, both `aa>150`. VirB5 failure 25.4% -> 7.5%;
slot_ready 70.9% -> 88.1%.

**Open.** Three candidate families (cl39, cl46, cl13) selected by lipobox, none
shown to do anything. cl46 is the strongest: 13 unique sequences over 46 Mash
clusters -- near-invariant protein on diverse backbones.

---

# MPF_F (F / tra) -- complete catalogue, no content criterion

**Solid.** Entry `PF11130` TraC_F_IV, 5/5 seeds, 93% of admissions independently
MPF_F. Slot TraG_N+1. Two negative slots (TraU+1, TrbC_Ftype+1), both with **zero**
named-family hits against 63 unique TraS at the target.

Only two negative slots exist: of the three anchor pairs reaching 30% adjacency,
one is the target's own side (`TraH -> TraG_N`, 96.1%).

**No content criterion, and lipobox is INVERTED.** Target 4.1% vs a negative slot's
6.1%. This is the strongest of the four measured instances of "criteria do not
transfer across classes": the criterion giving 20x separation in MPF_T points the
wrong way here.

Length also has no power: slot occupants are median 173 aa against a same-plasmid
background of 183, and P(slot longer than a random protein from its own plasmid) =
0.458 for TraS, 0.538 for unnamed. The apparent 177-vs-126 signal was an artefact
of negative-slot survivors being a selected remainder (12-19% survival).

---

# MPF_I (IncI) -- the cleanest slot

**Solid.** Entry `T4SS_I_traU`, 5/5 seeds, 0.46% cross-class. Slot TraY+1 with
**92.7% decontamination survival** -- the cleanest of the three classes (MPF_T's
own target loses 57.9%). Three negative slots (traM+1, traP+1, trbA+1), all with
**zero** named-family hits against 204 unique ExcA at the target.

`slot_ready__traY` 3,607 / 3,732 = **96.7%**, the highest of any class.

**Two ExcA figures, different denominators, both correct:**
59.3% (2,213/3,732 of admitted plasmids carry ExcA) and 96.5% (2,213/2,293 of all
PLSDB ExcA carriers are recovered by the entry criterion). Neither implies the other.

**No content criterion.** Hydrophobicity separates at only 2.3x against a 40.5%
background; `slot_admission` is `undetermined` on every row and `slot_maxkd` is
carried as a column, never a filter. Length fails outright (target 204 aa vs a
negative slot's 210).

**MPF_I is plasmid-specific.** Only 10 MPF_I ICEs exist, against 3,732 plasmids.

---

# MPF_FA (Gram-positive) -- entry criterion only, NO SLOT

**Solid.** Entry `T4SS_virb4 AND NOT(T/F/I) AND >=2 FA/FATA profiles` -> 2,391.
The `>=2` requirement discriminates at 84x over MPF_F background and 19x over
MPF_T; all five seeds pass it independently of how it was chosen.

`T4SS_virb4` alone was REJECTED at 79.5% cross-class contamination against a
pre-registered 30% bar. It admits 99.99% of MPF_T and 99.6% of MPF_F but only
**3.9% of MPF_I** -- a generic profile failing on an entire class.

**NO VALIDATED SLOT. This is the binding limitation.**
`conG+3` holds on ICEBs1 alone: 1 of 4 locatable seeds, with two occupants
ANTI-ORIENTED to the anchor and one an identified toxin-antitoxin system. The
class census found DUF4467 on **11 of 2,391** plasmids (0.5%), and within those 11
the modal offset from `FA_orf15` is **-1**, not ICEBs1's **+3** -- the one known
case is not the modal case.

`slot_status = no_slot_defined` on all 2,391 rows. A slot column validated on one
element would be empty, not caveated.

**Where to fix it: ICEs, not plasmids.** DUF4467 is 11 in plasmids and **98 in
ICEs** (all MPF_FA, all Firmicutes). Gram-positive conjugative elements are
predominantly ICEs. The plasmid database does not contain enough material.

---

# Cross-class results

**The VirB4-family ATPase is the entry criterion in every class**, and is mandatory
in every CONJScan plasmid class definition. MPF_T `PF03135`, MPF_F `PF11130`,
MPF_I via traU, MPF_FA conE.

**Criteria do not transfer across classes -- four measured instances:**

| instance | measurement |
|---|---|
| MPF_T lipobox | 56.1% in-slot vs 2.8% background -- 20x |
| MPF_F lipobox | **INVERTED** -- 4.1% target vs 6.1% negative |
| MPF_I hydrophobicity | 2.3x against a 40.5% background |
| MPF_FA tiering | no mandatory core exists to tier on |

**Plasmids and ICEs do not share an exclusion repertoire:**

| family | plasmids | ICEs |
|---|---|---|
| TraS | 2,541 | 2 |
| ExcA | 2,112 | 2 |
| Eex_IncN | 3,644 | 474 |
| TrbK | 218 | 10 |
| DUF4467 | 11 | **98** |
| EexR | **0** | **38** |

Only Eex_IncN has substantial presence on both, and even there ~8:1.

**Specificity, bounded by five negative sets** spanning 86 to 184,668 elements:
ICE 74.3% -> ICE-NA 21.0% -> IME 0.5% -> prophage 0.165% -> AICE 0.0% ->
**curated phage genomes 0.000%** (0 of 6,467, including P1, lambda, T4, Mu, P22).

---

# What is NOT established, for any class

* **These are gene clusters, not operons.** No promoter, terminator or
  co-transcription evidence anywhere in the pipeline.
* **8,526 named exclusion genes is a LOWER BOUND** -- members of already-named
  families only.
* **Surface exclusion is outside the frame.** The VirB6-anchored slot is an
  entry-exclusion instrument; `sfx` abuts traN, PrgA abuts PrgB, Sea1 abuts Asa1.
  Adding TraT and PrgA/Sea1 models yielded 7 sequences. A TraN-anchored frame
  would be a separate instrument.
* **No candidate has been shown to do anything.** The three MPF_T families are
  selected by a controlled procedure, not demonstrated.
* **The candidate pool is 1,734 unique** (was 2,514 before anchor decontamination),
  banded 947 large / 470 medium / 164 small / 153 singleton. The large band is
  unidentified -- resolving it needs AFDB/PDB, not our own 801 named proteins.

---

# Data

| file | rows |
|---|---|
| `data/release/v1.1/catalogue_MPF_T_v1.1.tsv` | 7,436 |
| `data/release/v1/catalogue_MPF_F_v1.tsv` | 11,120 |
| `data/release/v1.1/catalogue_MPF_I_v1.tsv` | 3,732 |
| `data/release/v1.1/catalogue_MPF_FA_v1.tsv` | 2,391 |
| `data/release/v1.1/annotation_MPF_T_F_v1.1.tsv` | 18,381 |
| `data/release/v1.1/slot_occupant_families.tsv` | 20,349 |
| `data/anchors/mpf_in_ice.tsv` | 9,602 |
| `data/anchors/mpf_in_refseq_phage.tsv` | 6,467 |

Plasmid release mirrored at `hukuang/Exclusion_project` (HuggingFace, private).
