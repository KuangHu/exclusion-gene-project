# MPF_FA §2 -- entry criterion

## Adopted

```
T4SS_virb4  AND NOT (MPF_T or MPF_F or MPF_I)  AND  >=2 FA/FATA profiles
    -> 2,392 plasmids (3.3% of PLSDB)
```

Jobs 25991904 (contamination), 26020027 (residual identity).

**Cross-class contamination is EXCLUDED BY CONSTRUCTION, not measured.** The
`NOT (T or F or I)` clause removes those plasmids by definition. Reporting this as
"0% cross-class" would read as "measured and found clean", which it is not, and no
such figure is quoted anywhere.

## `T4SS_virb4` alone was REJECTED

Rule fixed before the run: contamination <= 30% usable, > 30% not.

| criterion | admitted | %db | contamination |
|---|---|---|---|
| `T4SS_virb4` | 23,191 | 32.0% | **79.5%** |
| `T4SS_t4cp1` | 27,039 | 37.3% | 75.9% |
| `T4SS_t4cp2` | 29,571 | 40.8% | 73.0% |
| `T4SS_tcpA` | 3,917 | 5.4% | 45.5% |

The three existing catalogues cover 22,079 of 72,438 plasmids (30.5%), so 30.5% is
the chance baseline. At 79.5% `T4SS_virb4` is ENRICHED for already-classified
plasmids, not neutral.

### The generic/class-specific gap is systematic, not a single-seed artefact

`T4SS_virb4` admits:

| class | admitted | share of that class |
|---|---|---|
| MPF_T | 7,435 / 7,436 | **99.99%** |
| MPF_F | 11,071 / 11,120 | **99.6%** |
| MPF_I | 145 / 3,732 | **3.9%** |

The known observation was that generic `T4SS_virb4` scores ZERO on R64 while
class-specific `T4SS_I_traU` scores 1,460. At class scale the generic profile
fails on essentially the ENTIRE IncI class while saturating two others. A generic
mandatory profile is not a safe entry criterion by default.

## The `>=2` threshold, and how it was chosen

Requiring at least TWO FA/FATA profiles. Measured against MPF_T and MPF_F
background samples (n=2,000 each), which is what makes the residual rate readable:

| threshold | residual | MPF_T bg | MPF_F bg | vs MPF_T | vs MPF_F |
|---|---|---|---|---|---|
| >=1 | 73.9% | 7.8% | 35.2% | 9.5x | **2.1x** |
| **>=2** | **50.4%** | **2.6%** | **0.6%** | **19.4x** | **84.0x** |
| >=3 | 34.4% | 2.5% | 0.0% | 13.8x | inf |

**`>=1` is nearly worthless**: MPF_F plasmids carry one FA/FATA profile 35.2% of
the time. Without the MPF_F background the residual's 73.9% would have read as
strong evidence. The discrimination is at `>=2`.

**The threshold was selected AFTER seeing this table.** That is post-hoc. It has
one independent check, which the seeds provide because they played no part in
choosing it:

| seed | FA | FATA | total | `>=2` |
|---|---|---|---|---|
| ICEBs1 | 5 | 1 | 6 | PASS |
| pCF10 | 0 | 8 | 8 | PASS |
| pAM373 | 3 | 2 | 5 | PASS |
| pLS20 | 1 | 1 | 2 | PASS |
| pAD1 | 1 | 1 | 2 | PASS |

5/5. `pLS20` and `pAD1` sit exactly at 2 -- the same two the class-assignment rule
left UNASSIGNED.

## CAVEAT 1: the seed check is on SEED FILES, not PLSDB

5/5 was measured on the GenBank seed records (job 25974613), NOT on PLSDB.
**pCF10 is absent from PLSDB and ICEBs1 is a chromosomal ICE.** A PLSDB-based seed
check would read 0/5 for reasons that have nothing to do with the criterion --
the A12 failure mode, already seen when MPF_F's entry reported 0/5 because PLSDB
is RefSeq and the seeds were INSDC.

The criterion is therefore validated on seeds and APPLIED to PLSDB. Those are not
the same population.

## CAVEAT 2: the slot control is a single element

`PF14729` fires on ICEBs1 `conJ` and on nothing else.

| class | named family at the slot | recoverable across the class |
|---|---|---|
| MPF_T | TrbK / `NF033894` | 224 |
| MPF_F | TraS / `PF10624` | 2,725 |
| MPF_I | ExcA / `NF033891` | 2,213 |
| **MPF_FA** | DUF4467 / `PF14729` | **1** |

No exclusion family fires on pLS20, pCF10, pAD1 or pAM373. This is MPF_FA's
weakest link and is not offset by the entry criterion being sound.

## Recorded, not explained

1,239 of the 4,747 residual plasmids (26.1%) carry NO FA/FATA profile at all and
are excluded by the `>=2` requirement.

## Single-profile alternates, not adopted

`FA_orf14` is the cleanest single profile -- 1,260 admitted, **0.3%**
contamination, 3/5 seeds (ICEBs1, pLS20, pAD1). Rejected as sole entry because it
reaches only 3/5.

Three-way disjunctions reaching 5/5 exist (e.g. `FATA_prgHb + FA_orf13 +
FA_orf14`, worst arm 3.9%) and were **deliberately not adopted**: with 5 seeds and
34 profiles a 3-way cover is close to guaranteed by chance, so its 5/5 would be
fitted to the seed set rather than measured. It would also mix FA and FATA arms,
merging two classes CONJScan keeps separate.
