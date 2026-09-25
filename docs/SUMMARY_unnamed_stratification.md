# Stratifying the unnamed slot occupants

Of the slot occupants carrying no named exclusion family, how many are genuinely
new? The starting number was 9,412 records and had no interpretable meaning until
split.

---

# 1. The answer

**9,412 records are 1,734 distinct proteins.** Mutually exclusive by
construction -- no category double-counts.

| class | unique | records |
|---|---|---|
| entry exclusion (named) | 832 | 8,561 |
| surface exclusion (named) | 7 | 12 |
| both entry + surface | 3 | 3 |
| **ANCHOR family -- was miscounted as candidate** | **739** | 2,378 |
| **unassigned candidate** | **1,734** | 6,973 |
| total | 3,315 | |

Previously reported as 2,514; the difference is 739 anchor proteins plus 41
newly-detected known-family members.

---

# 2. Two defects, both mine

## 2.1 A quarter of the "candidates" were T4SS machine components

The MPF_T slot iteration decontaminated nominations against all 14 anchor
families (removing 57.9%). The annotation pass that produced
`slot_occupant_families.tsv` classified occupants using ONLY the six exclusion
families and skipped decontamination. **One concept, two implementations, two
standards.**

Leaked into the pool: VirB6 460, VirB5 318, VirB8 125, VirB9 121, VirB4b 82,
VirD4 65, VirB4 56, VirB1 49.

Two independent methods agree on the magnitude -- a full Pfam-A scan gives
562/2,297 = **24.5%**, phmmer decontamination gives **22.3%**. The figure does not
depend on either tool's threshold.

**Fixed at the definition layer**: decontamination is now part of what "slot
occupant" means, not a step some call sites perform. Fixing the definition rather
than the call site is what prevents recurrence.

## 2.2 Model choice, not divergence

TrbK was scanned only as `TIGR04359`. Pfam's own TrbK model `PF20084` finds 35
members where TIGR04359 finds 49 -- same family, two models, different coverage.
Both are now used. 41 sequences moved from "unnamed" to named on this basis.

---

# 3. Surface exclusion: a FRAME BOUNDARY, not a coverage gap

Adding `PF05818` (TraT) and a `PrgA_Sea1` profile (pCF10 prgA + pAD1 sea1, both
891 aa) yields **7 unique sequences, 12 records.**

That is not a weak model. TraT is detectable and demonstrably present in F and
R100 -- it simply does not occupy VirB5+1. The position data agrees: **sfx abuts
traN, PrgA abuts PrgB, Sea1 abuts Asa1.**

**Conclusion: the VirB6-anchored slot frame is an ENTRY-exclusion instrument.**
Surface exclusion requires an independent TraN-anchored frame. This is the
applicable boundary of the method, not a defect in it.

Since the project scope is entry exclusion, the practical consequence is small:
the 7 surface hits are flagged and removed, and a TraN frame is an optional
extension rather than a hole to fill.

**Recorded as not covered:** IncC `sfx` (`sequence: status: NOT_RECOVERED` in
`data/elements/IncC.yaml`) and pLS20 `ses` (absent from the GenBank record;
`gene_targets.tsv` notes "essentially no sequence family").

---

# 4. What each layer removed

| layer | removed | note |
|---|---|---|
| known family, sequence (HMM E<=1e-3) | 217 | saturates at E<=1e-2 |
| known family, structure (3Di) | 612 | TrbK 323, TraS 262, ExcA 15, Eex_IncN 12 |
| anchor contamination | 739 | defect 2.1 |
| isolated on a single element | 160 | only 9 overlap the structural layer |

## 4.1 Sequence search is exhausted; structure is not

Relaxing the HMM cutoff from 1e-2 to 1 adds **7 sequences out of 2,514**. The
hold-out shows why: recovery is *identical* at 1e-5, 1e-3 and 1e-2 for every
family. The failure is profile-level, not threshold-level, so threshold relaxation
is closed as a strategy.

Structure then found **3x more** known-family members than sequence (612 vs 217).

## 4.2 The result this line rests on

Two unrelated experiments give the same ordering across four families:

| family | hold-out recovery @10% | structural recruits |
|---|---|---|
| **TrbK_RP4** | **75.0%** (worst) | **323** (most) |
| TraS | 93.8% | 262 |
| Eex_IncN | 94.7% | 12 |
| ExcA | 100.0% (best) | 15 |

TrbK -- RP4/R751 at 41% identity while the rest of Tra2 is 75-92% -- is hardest
for sequence and richest for structure. ExcA is the reverse.

**Consequence.** A genuinely new family is by definition sparsely sampled, i.e.
TrbK-shaped. **So the sensitivity applicable to the candidate pool is 50-75%, not
the 94.9% weighted average.** The weighted figure hides precisely the case that
matters.

## 4.3 3Di calibration

TM-score was impossible: ProstT5 yields 3Di SEQUENCES, not coordinates, so there
is no `_ca` file and TM-score is undefined without folding all 3,098 proteins. The
threshold was instead calibrated on our own labelled data.

| | |
|---|---|
| same-family pairs | 217,080, median E **3.55e-04** |
| different-family pairs | 104,558, median E **4.21** |
| cut E<=2.11e-02 | keeps 90% same-family, admits 2.9% different-family |

Four orders of magnitude between medians; the cut is sharp, so the negative result
(73.4% with no structural hit) is readable.

---

# 5. Claims that did not survive

| claim | outcome |
|---|---|
| "9,412 candidates" | 2,514 unique, then **1,734** |
| isolated ORFs are most of the remainder | **7%** (161/2,297) |
| cluster 21 is a large novel family (193 unique, 1,196 Mash clusters) | almost certainly **VirB6** -- TrbL was both the largest Pfam hit (324) and the largest contaminant (460) |
| the two filters are redundant | nearly independent: 612 structural, 160 isolated, **9 both** |

---

# 6. Method errors, recorded

1. hold-out used GA on freshly built HMMs -- they carry no gathering cutoff
2. **hold-out held out the LARGEST cluster** -- profile built from 5 sequences,
   tested on 380, produced a plausible **8.9%** that measured nothing. Fixed to
   hold out the divergent minority with `assert len(keep) > len(hold)`
3. `easy-search` given FASTA instead of the prebuilt 3Di databases
4. CPU foldseek build asked for a GPU (CUDA build is a separate 608 MB package)
5. TM-score requested from coordinate-free databases
6. **Pfam-A never scanned until this round** -- "unnamed" had only ever meant "no
   hit to six families", never "no hit to anything"
7. jobs sized by guess: 16 CPU/96 GB/12 h for a 2,200-sequence job
   (unschedulable), then a 5.6 h scan against a 4 h wall (killed)

**Item 2 is the dangerous class: a successful run producing a meaningless
number.** A crash is safe by comparison.

---

# 7. Next, in priority order

1. **Rerun the size bands on the corrected 1,734** -- REQUIRED. The existing bands
   (large 894 / medium 460 / small 172) were computed on the contaminated pool.
   The **medium band, 5-30 unique**, is the plausible magnitude for a real family:
   TrbK is 218 in plasmids, DUF4467 98, EexR 38. Large band is probably category
   3; small band is noise.
2. **Structural search against AFDB/PDB** -- higher priority than a TraN frame.
   The structural layer so far compared only against our own 801 named proteins,
   so "no structural hit" means "not like a known exclusion family", NOT "no
   functional assignment". Only an AFDB search separates genuinely uncharacterised
   proteins from ones with a different known function. The ProstT5 3Di
   representation already exists, so the cost is low.
3. **Slot divergence rate** -- independent of 1 and 2. Compute each slot's
   divergence relative to its backbone across homologous element groups, with
   TrbK/TraS slots as positive controls and VirB6/VirB11 as negatives. This gives
   the 1,734 a RANKING rather than a binary label.
4. **TraN-anchored frame** -- optional under an entry-only scope.
