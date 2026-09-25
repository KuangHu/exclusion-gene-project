# Entry-exclusion genes in conjugative elements: plasmids and ICEs

Final results. Two element classes carry the conjugation machinery and the
exclusion genes that go with it: **plasmids** (PLSDB) and **ICEs**. Phage work is
closed and appears here only as the control that bounds the specificity.

One pipeline, one set of frozen criteria, applied to both:

```
entry criterion -> architecture -> slot -> catalogue -> POST-HOC exclusion families
```

| class | entry criterion |
|---|---|
| MPF_T | `PF03135` CagE_TrbE_VirB (VirB4) |
| MPF_F | `PF11130` TraC_F_IV |
| MPF_I | `T4SS_I_traU` |
| MPF_FA | `T4SS_virb4` AND NOT(T/F/I) AND >=2 FA/FATA |

Gene calling Pyrodigal anon mode (`ANON=True, MIN_GENE=90, CLOSED=True,
TRANSL_TABLE=11`, A1). Search at gathering thresholds only; **E-value thresholds
are prohibited** (A13) because they scale with database size.

---

# 1. PLASMIDS -- PLSDB 2024_05_31_v2

**24,470 distinct plasmids catalogued. 8,526 named entry-exclusion genes at the slot.
9,412 unnamed slot occupants as a candidate set.**

| class | catalogued | slot-resolvable | named eex gene | unnamed candidate |
|---|---|---|---|---|
| MPF_T | 7,436 | 6,548 | 3,862 (59.0%) | 2,655 |
| MPF_F | 11,120 | 7,803 | 2,541 (32.6%) | 5,262 |
| MPF_I | 3,732 | 3,607 | 2,112 (58.6%) | 1,495 |
| MPF_FA | 2,391 | none (no slot) | 11 (plasmid-level) | -- |
| **total** | **24,679 rows** | 17,958 | **8,526** | **9,412** |

| family | model | n | class |
|---|---|---|---|
| Eex_IncN | NF033894 | 3,644 | MPF_T |
| TraS | PF10624 | 2,541 | MPF_F |
| ExcA | NF033891 | 2,112 | MPF_I |
| TrbK_RP4 | TIGR04359 | 218 | MPF_T |
| DUF4467 | PF14729 | 11 | MPF_FA |
| EexR | NF041429 | **0** | -- |

## Slots, validated by negative controls

Named exclusion families appear at the target slot and at ZERO of five negative
slots, in all three classes that have one:

| class | slot | target | negative slots |
|---|---|---|---|
| MPF_T | VirB5+1 | 372 | 0 / 0 / 0 (VirB8+1, VirB9+1, VirB10+1) |
| MPF_F | TraG_N+1 | TraS 63 unique | 0 / 0 (TraU+1, TrbC_Ftype+1) |
| MPF_I | TraY+1 | ExcA 204 unique | 0 / 0 / 0 (traM+1, traP+1, trbA+1) |

MPF_FA has **no validated slot**: `conG+3` holds on ICEBs1 alone, and DUF4467
occurs on 11 of 2,391 admitted plasmids at an offset that does not reproduce.

## Candidate families, MPF_T only

Selection by lipobox >90% with >=20 hits -- zero families at all three negative
slots, positive control (Eex_IncN) recovered:

| family | records | Mash clusters | unique | lipobox | median aa |
|---|---|---|---|---|---|
| cl46 | 108 | 46 | **13** | 98.8% | 108 |
| cl13 | 106 | 69 | 58 | 91.8% | 114 |
| cl39 | 27 | 16 | 10 | 100.0% | 66 |
| *cl25* | *173* | *135* | *92* | *96.2%* | *77* | *(= Eex_IncN, control)* |

Nothing is clonal (collapse 1.1-2.3x). **cl46 is the strongest**: 13 unique
sequences across 46 Mash clusters -- a near-invariant protein on diverse
backbones, the signature of functional constraint.

**These are candidates, not discoveries. None has been shown to do anything.**

---

# 2. ICEs -- 9,602 elements

## The controls came free with the dataset, and all passed

| element_type | n | any MPF | relaxase | expectation |
|---|---|---|---|---|
| **ICE** | 2,634 | **74.3%** | 96.5% | self-transmissible -> needs MPF |
| **IME** | 1,490 | **0.5%** | 92.8% | mobilizable -> relaxase, no MPF |
| **AICE** | 86 | **0.0%** | 5.8% | FtsK/SpoIIIE, no T4SS |
| NA | 5,391 | 21.0% | 73.3% | untyped in source |

**ICE vs IME is 150-fold.** Relaxase, present in ~95% of both, detects none of it.
This is the defining biological difference between the two element classes, and
MPF recovers it.

**ICEBs1** (`NC_000964:529362..549932`) types **MPF_FA** with 6 FA/FATA profiles
and a DUF4467 hit on `conJ` -- the exact protein the MPF_FA slot control rests on.

**This is the evidence the frozen criteria are not overfit to PLSDB.**

| class | n |
|---|---|
| MPF_F | 1,323 |
| MPF_T | 919 |
| MPF_FA | 838 |
| MPF_I | 10 |

**MPF_I is essentially absent from ICEs** (10, vs 3,732 in plasmids) -- IncI
conjugation is plasmid-specific.

## Exclusion genes: 623 elements (6.5%), 93.1% of them MPF-positive

| family | ICEs | plasmids |
|---|---|---|
| Eex_IncN | 474 | 3,644 |
| DUF4467 | 98 | 11 |
| EexR | **38** | **0** |
| TrbK_RP4 | 10 | 218 |
| TraS | 2 | 2,541 |
| ExcA | 2 | 2,112 |

## Two recorded hypotheses became measurements

**EexR.** `mpf_f_entry.md` recorded: *"SXT/R391-family elements being chromosomal
ICEs and under-represented in a plasmid database -- is a HYPOTHESIS, not a
measurement."* EexR fires **0 times across 20,269 plasmids** and **38 times in
ICEs** -- all `element_type=ICE`, all MPF_F, concentrated in *Proteus* (27),
*Shewanella* (4), *Vibrio* (4). Confirmed by changing dataset, not by relaxing a
threshold.

**DUF4467.** 13 in plasmids -> **98 in ICEs**, all MPF_FA, all Firmicutes. MPF_FA,
the weakest class in the plasmid release and the one with no validated slot, is
far better represented in ICEs. **The MPF_FA slot question should be attacked in
ICEs, not plasmids.**

---

# 3. Specificity, bounded by the phage controls (line closed)

| set | n | MPF-positive |
|---|---|---|
| ICE (element_type=ICE) | 2,634 | **74.3%** |
| ICE (NA) | 5,391 | 21.0% |
| ICE (IME) | 1,490 | 0.5% |
| prophage island calls (16/32 shards) | 184,668 | 0.165% |
| ICE (AICE) | 86 | 0.0% |
| **curated phage genomes (RefSeq viral)** | **6,467** | **0.000%** |

Five independent negative sets spanning 86 to 184,668 elements. **Zero false
positives in 6,467 complete phage genomes**, including P1 (the phage-plasmid,
94,800 bp -- plasmid replication but not self-transmissible, correctly negative),
lambda, T4, Mu and P22.

**Relaxase is a length artefact; MPF is not.** In phages, relaxase hit rate rises
with genome size -- 1.0% (<50 kb), 6.9% (50-100 kb), **26.1%** (100-200 kb) --
while MPF stays at exactly **0** in every band. MOB profiles pick up phage
nucleases at GA; the MPF criteria never fire on a non-conjugative element.

The 304 MPF-positive prophages (0.165%) are median 123 kb with 95% relaxase -- an
ICE signature. TIGER/Islander call islands from integrase + att sites, which ICEs
satisfy too. Read as **~0.17% ICE contamination in the prophage set**, not as
conjugative phages.

---

# 4. Standing limits

* **NOT operons.** No promoter, terminator or co-transcription evidence anywhere.
  Anchor sets in measured adjacency.
* **8,526 is a lower bound** -- members of six already-named families. The 9,412
  unnamed slot occupants are a candidate container, not absence.
* **Only MPF_T has a usable content criterion** (lipobox). In MPF_F it is
  INVERTED (target 4.1% vs a negative slot's 6.1%); in MPF_I it is negligible.
  Length, KD and TM count all fail or are unavailable. The criterion search is
  closed.
* **MPF_FA has no validated slot.**
* Untried: transcriptional independence (own promoter, Rho-independent
  terminator) -- the only axis orthogonal to both sequence and position.

# 5. Data

| file | rows |
|---|---|
| `data/release/v1.1/catalogue_MPF_{T,I,FA}*.tsv`, `data/release/v1/catalogue_MPF_F_v1.tsv` | 24,679 |
| `data/release/v1.1/annotation_MPF_T_F_v1.1.tsv` | 18,381 |
| `data/release/v1.1/slot_occupant_families.tsv` | 20,349 |
| `data/anchors/mpf_in_ice.tsv` | 9,602 |
| `data/anchors/mpf_in_refseq_phage.tsv` | 6,467 (control) |

Plasmid release mirrored at `hukuang/Exclusion_project` (HuggingFace, **private**).
