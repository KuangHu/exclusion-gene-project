# The four MPF classes searched in ICEs and prophages

Jobs 26079721 (ICE, complete) and 26081820 (prophage, **stopped at 16/32 shards**).

Both datasets were typed with the FROZEN entry criteria from
`data/release/v1.1`, the frozen Pyrodigal settings (`ANON=True, MIN_GENE=90,
CLOSED=True, TRANSL_TABLE=11`, A1) and gathering thresholds only (A13), so the
counts are directly comparable to the plasmid release rather than being a
parallel pipeline.

| | |
|---|---|
| MPF_T | `PF03135` anchor_CagE_TrbE_VirB (VirB4) |
| MPF_F | `PF11130` mpff_TraC_F_IV |
| MPF_I | `T4SS_I_traU` |
| MPF_FA | `T4SS_virb4` AND NOT(T/F/I) AND >=2 FA/FATA profiles |

---

# 1. ICE / ciMGE database -- COMPLETE

9,602 elements, 629,173 proteins called.

## The controls came free with the dataset, and all passed

| element_type | n | any MPF | relaxase | expectation |
|---|---|---|---|---|
| **ICE** | 2,634 | **74.3%** | 96.5% | self-transmissible -> should carry MPF |
| **IME** | 1,490 | **0.5%** | 92.8% | mobilizable -> relaxase but NO MPF |
| **AICE** | 86 | **0.0%** | 5.8% | FtsK/SpoIIIE, no T4SS |
| NA | 5,391 | 21.0% | 73.3% | untyped in the source table |

**ICE vs IME is 150-fold.** This is the defining biological difference between the
two element classes, and relaxase -- present in ~95% of BOTH -- cannot detect it.
MPF can. AICE at 0.0% is the second independent negative.

**Positive control: ICEBs1** (`NC_000964:529362..549932`, 20,571 bp) types as
**MPF_FA** with 6 FA/FATA profiles, and its DUF4467 hit is `conJ` -- the exact
protein the MPF_FA slot control rests on. **PASS.**

This is the first evidence the frozen criteria are not overfit to PLSDB.

## Classes

| class | n | share |
|---|---|---|
| MPF_F | 1,323 | 13.8% |
| MPF_T | 919 | 9.6% |
| MPF_FA | 838 | 8.7% |
| MPF_I | 10 | 0.1% |
| MPF_T;MPF_F | 8 | -- |
| MPF_T;MPF_I | 1 | -- |

**MPF_I is essentially absent from ICEs** (10, vs 3,732 in plasmids). IncI
conjugation is a plasmid-specific system.

## Two recorded hypotheses became measurements

**EexR (`NF041429`).** `mpf_f_entry.md` recorded: *"SXT/R391-family elements being
chromosomal ICEs and under-represented in a plasmid database -- is a HYPOTHESIS,
not a measurement."*

| | PLSDB plasmids | ICEs |
|---|---|---|
| EexR | 4 db-wide, **0** in any slot_ready set | **38** |

All 38 are `element_type=ICE`, all are MPF_F, concentrated in *Proteus* (27),
*Shewanella* (4), *Vibrio* (4). The family that fired zero times across 20,269
plasmids is present and class-consistent in the right dataset. **Hypothesis
confirmed by changing dataset, not by relaxing a threshold.**

**DUF4467 (`PF14729`).** 13 in plasmids -> **98 in ICEs**, all MPF_FA, all
Firmicutes. MPF_FA -- the weakest class in the plasmid release, the one with no
validated slot -- is far better represented here. **The MPF_FA slot question is
better attacked in ICEs than in plasmids.**

## Exclusion genes in ICEs

623 of 9,602 elements (6.5%) carry a named exclusion family; **93.1% of those are
also MPF-positive**, so the exclusion genes track the conjugation systems rather
than scattering.

| family | n |
|---|---|
| Eex_IncN | 474 |
| DUF4467 | 98 |
| EexR | 38 |
| TrbK_RP4 | 10 |
| TraS | 2 |
| ExcA | 2 |

---

# 2. Prophages -- STOPPED at 16/32 shards

**Stopped by request after 184,668 of ~368,615 sequences (50.1%).**

The partial result is a reliable estimate: the 16 completed shards are
independent and their MPF rates range only **0.123%-0.232%** (mean 0.165%).
Stopping at half the dataset cost essentially no precision.

| | value |
|---|---|
| sequences scanned | 184,668 (50.1%) |
| MPF-positive | **304 (0.165%)** |
| relaxase-positive | 5,916 (3.204%) |
| classes | MPF_F 179, MPF_FA 98, MPF_T 26, MPF_I 1 |
| exclusion families | Eex_IncN 67, DUF4467 13, TrbK_RP4 13, ExcA 3, EexR 1 |

## Prophages are the negative control, and they behave like one

**0.165% against ICEs' 74.3% -- 450-fold.** With IME (0.5%) and AICE (0.0%), the
typing now has three independent negative sets at three different scales.

## The 304 MPF-positive prophages are most likely MISANNOTATED ICEs

| | MPF-positive | other informative rows |
|---|---|---|
| median length | **123,353 bp** | 26,483 bp |
| relaxase-positive | **95%** (288/304) | -- |

**4.7x longer than typical, and 95% carry a relaxase.** Long + MPF + relaxase is
an ICE signature, not a phage one. TIGER/Islander call mobile islands from
integrase + att sites, which ICEs satisfy just as well as prophages do.

**The reading is NOT "prophages are conjugative."** It is that the prophage set
carries ~0.17% contamination from conjugative elements, and this search
identifies them. The run doubles as QC on that dataset.

Relaxase at 3.2% is 20x the MPF rate -- the same asymmetry seen in ICEs, where
relaxase could not separate ICE from IME but MPF could. Relaxase-only elements are
mobilizable in trans or carry MOB relics; the full MPF machinery is what is
genuinely rare.

---

# 3. Wording

These are **gene clusters, not verified operons**. No promoter, terminator or
co-transcription evidence anywhere in the pipeline -- anchor sets in measured
adjacency, exactly as in the plasmid release.

# 4. Outputs

| file | rows |
|---|---|
| `data/anchors/mpf_in_ice.tsv` | 9,602 (complete) |
| `/global/scratch/.../mge_search/phage/phage_mpf.shard0{00..15}.tsv` | 16 shards, informative rows only |

Shards 016-031 were never run. Resuming is a matter of resubmitting
`slurm/150_phage.sh` with `--array=16-31`.
