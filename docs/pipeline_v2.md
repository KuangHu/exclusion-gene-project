# Pipeline v2 — MPF_T slot-occupancy scan

**Status:** design frozen, not yet run. Supersedes v1, which rested on the
assumption that exclusion genes hide in intergenic gaps — refuted by the
four-plasmid GO test (`docs/slot_hypothesis_GO.md`).

**Core hypothesis (passed 4/4 positive controls):** in MPF_T conjugative
plasmids the entry-exclusion gene occupies the position immediately upstream of
the VirB6 homolog.

| plasmid | Inc | coords | aa | annotation | lipobox | mature |
|---|---|---|---|---|---|---|
| RP4 | IncP-α | 27450..27659 | 69 | `trbK` | `LAG↓C23` | 47 |
| R751 | IncP-β | 26340..26567 | 75 | `trbK` | `VAG↓C22` | 54 |
| pKM101 | IncN | 6533..6760 | 75 | `eex` | `LVA↓C15` | 61 |
| R388 | IncW | 5883..6113 | 76 | unnamed, "EexN family lipoprotein" | `LSA↓C17` | 60 |

**Why this is a real result:** two families that homology cannot connect occupy
the same slot. Within `TrbK_RP4` (TIGR04359) 37.8%; within `Eex_IncN` (NF033894)
40.7%; **between families 21–25% = background.** The slot is conserved across
families that homology detection cannot link.

---

## 0. Database

**PLSDB only (72,360 dereplicated complete plasmids). No IMG/PR at this stage.**

The reason is not annotation quality — we re-run Pyrodigal, so the deposited
annotation is irrelevant. The reason is **completeness**:

* a truncated or mis-assembled record that loses the VirB6 anchor simply does
  not enter the set (acceptable)
* a record that keeps VirB6 but loses the upstream gene **manufactures an empty
  slot** (disastrous, and indistinguishable from a genuinely empty one)

The second failure lands exactly along the axis the signal is most sensitive to.
IMG/PR's 699,973 records include many metagenome-predicted plasmids with no
completeness guarantee. Add it only after an occupancy baseline exists from
PLSDB, with a completeness filter, and **report the two sources separately**.

### 0.1 Statistical denominator: cluster first

**72,360 records are not 72,360 independent observations.** Clinically important
plasmids are sequenced repeatedly; 500 near-identical IncP plasmids all carrying
`trbK` is approximately one independent observation. Raw frequencies would be
dominated by oversampling and the hypergeometric denominator would be wrong.

**Decision: Mash/sourmash ANI clustering (threshold starting at 95%) is the
statistical unit.**

* report occupancy **both per record and per cluster**
* COPLA/PTU is too expensive (72k pairwise ANI + HSBM); use it only as
  classification background for final candidates
* PlasmidFinder Inc typing is too coarse and leaves many plasmids untyped

**Changing the threshold re-runs every statistic. Fix it before running.**

---

## 1. Determining the anchor set (prerequisite — do this first)

**Do not copy any enumerated Pfam accession. Determine it empirically.**

Confirmed:

* `PF04610` = TrbL (VirB6) — GA threshold already validated on 20 seed proteins
* `PF03135` = CagE_TrbE_VirB (VirB4)
* `PF03524` = CagX (VirB9)
* `PF07996` = "Type IV secretion system proteins"
* `PF07916` = TraG_N (VirB6 of MPF_F)

**RESOLVED 2026-09-02** by `scripts/80_determine_anchors.py` — see
`data/anchors/anchor_set.tsv`. Full set below; two findings changed the design.

### 1.1 Empirical procedure (minutes)

```
1. take the Pyrodigal proteins of RP4 / R751 / pKM101 / R388
2. scan against the complete Pfam-A.hmm at GA thresholds
3. record which families hit each known component:
     IncP:  TrbB TrbC TrbD TrbE TrbF TrbG TrbI TrbJ TrbL TrbN
     IncN:  TraB TraC TraD TraE TraF TraN TraO TraG TraL TraM
4. the families that hit = the anchor set
5. record which components have NO usable Pfam family -> blind spots for
   slot definition later
```

The anchor set obtained this way carries its own positive-control validation.

### 1.1a Result

| VirB | Pfam | acc | threshold | role |
|---|---|---|---|---|
| VirB1 | SLT | PF01464 | GA | annotation |
| VirB2 | TrbC | PF04956 | GA | annotation |
| VirB3 | VirB3 | PF05101 | GA | annotation |
| VirB4 | CagE_TrbE_VirB | **PF03135** | GA | **entry criterion** |
| VirB4 | VirB4_TrbE_N | PF27097 | GA | secondary / sensitivity check |
| VirB5 | T4SS | **PF07996** | **E ≤ 1e-5** | **empty-slot test** |
| VirB6 | TrbL | **PF04610** | GA | **slot downstream boundary** |
| VirB7 | — | — | — | **BLIND SPOT** |
| VirB8 | VirB8 | PF04335 | GA | annotation |
| VirB9 | CagX | PF03524 | GA | annotation |
| VirB10 | TrbI | PF03743 | GA | annotation |
| VirB11 | T2SSE | PF00437 | GA | annotation |

**Finding 1 — the entry criterion survives, verified.** `PF03135` admits all four
controls at GA (RP4 97, R751 90, pKM101 99, R388 121), so §3.1's choice holds.
`PF27097` (VirB4_TrbE_N) also hits all four and scores *higher* on IncP (RP4 114
vs 97); it is not needed for entry but is kept as a sensitivity check on the
VirB4-negative remainder.

**Finding 2 — VirB5 cannot use GA, and this was a live blocker.** At GA (24.7)
`PF07996` detects the VirB5 position only in pKM101 (87.4) and R388 (122.9);
**RP4 `trbJ` scores 17.7 and R751 22.2, so both IncP controls fail.** By E-value
all four are strongly significant: 7.9e-07, 3.4e-08, 3.5e-28, 4.6e-39. Threshold
therefore set to **E ≤ 1e-5** for this anchor only.

**The four scores are BIMODAL, not a continuum** — and this changes what the
threshold is doing:

```
RP4    18  ┐ below GA 24.7, E-value only
R751   22  ┘
                     <- 65-point gap, nothing in between
R46    87  ┐ comfortably above GA
R388  123  ┘
```

PF07996's ability to detect the VirB5 position **fails qualitatively on IncP**
rather than degrading gradually, and GA sits squarely in the empty gap. So this
is not the routine "threshold slightly too strict" case: on IncP the family is
essentially non-functional and only residual signal remains.

Consequence: **on IncP, `E ≤ 1e-5` carries the entire VirB5 determination with
one order of magnitude of margin** (RP4 at 7.9e-07). The `aa > 150` conjunct is
therefore not a safety net on IncP — it is the primary discriminator.

**Reporting requirement:** every VirB5 call must record whether it was reached by
GA or by E-value alone, and occupancy must be stratified on it. Derivable
post-hoc from the recorded bitscore against GA = 24.70, so no re-run is needed.
An IncP subgroup with different occupancy could be biology or could be VirB5
detection failure, and the two must be separable.

Had GA been used, the empty-slot test would have been *undefined on IncP* — an
unrecognised `trbJ` upstream of VirB6 is indistinguishable from an unrecognised
occupant, which is exactly the failure §4.1 names. `TIGR02762` was tested as an
alternative and rejected: it is TraL, not VirB5, and hits none of the four.

**Finding 3 — VirB7 is a blind spot.** No Pfam family hits VirB7 in any control
(pKM101 `traN`, 48 aa — too small and divergent). Recorded per step 5; VirB7 must
not be used for slot definition.

**Slot architecture, confirmed identical in all four:**

```
[VirB5-ish 226-258 aa]  [SLOT: 69-76 aa exclusion gene]  [VirB6 342-572 aa]
```

The 3x length gap between the VirB5 position and the slot occupant is itself a
usable discriminator, independent of any HMM.

### 1.2 Non-Pfam anchors

`IPR027628` (DotA_TraY), `IPR000515` (MetI-like) and `cd06261` **are not HMMs**
and must first be resolved to member-database signatures. `cd06261` is a CDD
PSSM, not in Pfam-A, and needs a separate fetch. **Do this now, not on discovering
it mid-run.** (Both are needed only when extending to MPF_I / MPF_FA; the MPF_T
main line does not use them.)

### 1.3 🔴 RED LINE: exclusion gene families must never enter the anchor set

**`TIGR04359` (TrbK_RP4), `NF033894` (Eex_IncN), `NF041429` (EexR),
`NF033891` (surf_exc_IncI1), `PF10624` (TraS), `PF14729` (DUF4467) must not be
used as anchors or as filter criteria.**

They may be used only for **post-hoc classification** ("is this candidate already
covered by an existing family?"). Once they enter the anchors it is label
leakage — the scan would "discover" the thing used to define the set.

### 1.4 🔴 RED LINE: anchors are defined by VirB number / Pfam family, never by IncP gene name

* `trbA` is an IncP-specific regulator with no VirB homolog and does not exist on
  pKM101 or R388 → using it as an anchor biases the set toward IncP
* `trbK` **is the target gene itself**
* the naming trap has already bitten once: TraG is VirB6 in F/IncC/SXT but the
  T4CP (VirD4) in RP4/R27

---

## 2. Step 0 — CDS calling

**`Pyrodigal -p anon`, whole plasmid, single pass.** Implementation:
`scripts/lib/orf_caller.py`; decision record: `docs/orf_caller_decision.md`.

Basis: all four positive controls hit (69/75/75/76 aa, every start exact). The
exclusion genes are not hidden in intergenic gaps — they are ordinary genes in a
distinctive position.

### 2.1 Hard constraints (assertions, not comments)

| # | assertion | why |
|---|---|---|
| A1 | must be `-p anon`; `normal` forbidden | normal has a hard floor, `sequence must be at least 20000 characters`; U09868 (12.4 kb) and J01566 (6.6 kb) raise. A batch run would **fail silently** on a subset |
| A2 | any new caller adapter must reproduce one known CDS exactly before use | orfipy's BED excludes the stop codon (off-by-3) and nearly invalidated the whole benchmark |
| A3 | any config value outside its enum → hard fail | historical bug: notes written into `restrict_inc_group` |
| A4 | any regex/query matching 0 records → hard fail, not empty return | historical bug: an undefined gene was silently skipped and reported as no-baseline |
| A5 | any `str_replace` yielding 0 replacements → hard fail | historical bug: the uniparc/crc64 edit did not match, leaving the layer with no provenance |
| A6 | one value only per `(gene_a, gene_b, metric)` triple | historical bug: TraS identity existed as both 17.7% and 27.9% |

**All six are now implemented and self-tested against the real historical bugs:**

| # | implementation | self-test |
|---|---|---|
| A1 | `scripts/lib/orf_caller.py` `assert_anon_only()` | trained mode unreachable |
| A2 | `scripts/lib/adapters.py` | catches a simulated orfipy off-by-3; refuses registration without a golden test |
| A3–A5 | `scripts/lib/assertions.py` + `scripts/74_validate.py` | all four historical bugs re-injected and caught |
| A6 | `scripts/lib/metrics.py` | rejects a duplicate triple and an unregistered metric name |
| **A7** | `scripts/lib/assertions.py` `is_blank` / `assert_clean_field` / `assert_no_cr`, wired into `74_validate.py` | catches `"\r"` reading as truthy; found 8 further CR-bearing tables on first run |

### A7 is the most valuable of the seven, because its bug produced a FALSE signal

A1–A6 all guard "matched nothing, returned empty" — an *absence*, which has some
chance of being noticed. A7's bug did the opposite. `csv.writer` defaults to
`\r\n`, so the last field of every row held a stray `\r`; an **empty** anchor
field was therefore truthy and **146/146 plasmids scored as anchor-bearing**.

**100% is precisely the value at which §6.1 says to abandon the direction.** Had
the smoke test not checked correctness as well as existence, a line-terminator
character would have produced "slot occupancy = 100%" across the library and the
project would have been killed by its own decision rule.

Rule: never apply bare truthiness to a field parsed from a file. Use
`is_blank()`, or compare explicitly against `""` / `None`.

A6 turned out to matter beyond bookkeeping. The TraS 17.7% / 27.9% pair was not a
contradiction — 36 identities over 203 alignment columns (gaps counted) versus 129
ungapped columns. Both were right under different denominators, stored under one
unqualified name. `metrics.py` makes the denominator part of the metric name.
**Consequence:** `humbert2019`'s "17%" was matched to the 17.7% figure; if that
paper used the ungapped denominator its comparator is 27.9% and the match was
coincidence. Carried as a `provenance_caveat` in `data/metrics.tsv`.

### 2.2 Six-frame enumeration: off the main line

Scope compressed to two situations, **neither on the MPF_T main line**:

| situation | why six-frame is needed | on the MPF_T line? |
|---|---|---|
| overlapping pair sharing a stop (R64 excA/excB) | one-per-stop cannot express it | no (MPF_I) |
| different-frame nesting (ColE1 mob3/mob4) | same | no |
| boundary genes Pyrodigal misses | recall (boundary case 2/4) | none of the four missed |

**The decision is analytic:** within the same window and length floor, six-frame
all-starts output is a **strict superset** of any one-per-stop caller's output.
So `recall(six-frame) ≥ recall(any caller)` always holds and no further
experiment can change the ordering. Experiment can only measure FP cost, and that
is done (stage 1 ≈ 3/window, stage 2 ≈ 7/window, GC correlation not significant
at n=10).

Reattach it when extending to MPF_I, via the **two-stage architecture**:

* **stage 1, occupancy** — unit = stop-to-stop locus, one per stop, FP ≈ 3/window.
  **Cross-plasmid enrichment testing happens only at this layer.**
* **stage 2, sequence identity** — enumerate all in-frame starts, only for loci
  that passed stage 1

The 2.3× all-starts cost then never enters the enrichment denominator.
excA/excB resolves naturally: for occupancy it is one occupied slot; the two
products are a stage-2 question.

---

## 3. Step 1 — set construction

```
1a  pyhmmer, all anchors (the set determined empirically in §1), over all 72k
    ENTRY CRITERION: VirB4 (PF03135) hit        <- the ONLY filter

1b  CONJscan (MacSyFinder v2, T4SS_typeT) on the VirB4-positive subset
    purpose: attach an MPF_T label that maps onto the literature
    VirB4-positive but CONJscan-negative -> review separately, do not discard

1c  report: fraction of VirB4-positives in which VirB6 (PF04610) is detectable
    <- the anchor's own recall. This must be reported.
```

### 3.1 Why the entry criterion is VirB4, not VirB6

**Using VirB6 as the entry criterion makes it permanently impossible to know how
often VirB6 detection failed — the denominator is defined away by the anchor
itself.**

And that failure is not random: what `PF04610` misses at GA is precisely the
**divergent** VirB6, hence divergent conjugation systems, hence the plasmids most
likely to carry novel exclusion genes. **Filtering on VirB6 systematically
removes exactly what the project is looking for.**

Three reasons for VirB4:

1. VirB4 is the most conserved T4SS component and remains detectable in
   divergent systems
2. literature standard: Guglielmini et al. 2014 (NAR 42:5715) use VirB4 ±20-gene
   windows for system detection
3. **it provides a denominator**: the VirB4⁺ / VirB6⁻ subset is re-checked at a
   lowered threshold

**Report BOTH numbers from the start, not on discovering a problem.** §1.1a
Finding 2 already showed GA dropping divergent members of a real anchor (VirB5
missed in both IncP controls at GA, strongly significant by E-value). Treat that
as the expected behaviour of GA rather than a surprise:

| stage | threshold | reported as |
|---|---|---|
| primary | GA (`--cut_ga`) | the headline count |
| sensitivity | E ≤ 1e-5 | the divergent-member count |

Both go in the output. The difference between them **is** the estimate of how
much VirB6 detection is missing, which is the whole reason VirB4 is the entry
criterion.

### 3.2 The entry criterion is a conjunction, not a union

"any of 12 anchors hits" would let a single incidental GA hit drag an unrelated
plasmid into the set — every HMM has its own FP rate at GA. Entry uses VirB4
alone; the other anchors are for annotation and slot definition, never for
filtering.

### 3.3 pyhmmer rather than the hmmsearch CLI

Multithreaded, no intermediate files, results identical to the CLI. The GA
threshold workflow is already validated on 8 models × 20 seed proteins; keep it
as is.

CONJscan runs only on the subset because **it is the real bottleneck**
(MacSyFinder is slow). Running pyhmmer with 1 versus 15 models over 72k costs
about the same; the bottleneck is Pyrodigal and I/O.

### 3.4 Why CONJscan is still needed

Its role is to **attach a label that maps onto the literature**, not to discover
anchors. Saying "MPF_T" in a paper must mean what Guglielmini/Cury mean.

* **not MOB-suite**: MPF typing is marker-based, coarser than CONJscan, and gives
  no component coordinates
* **not COPLA for filtering**: too expensive, and it does PTU classification
  rather than component localisation

Keep the union under review: CONJscan calls systems by quorum and accessory
components are frequently absent (GitHub issue #81 reports misses), while most
virB skeleton members are accessory. **Do not rely on CONJscan component hits
alone to define the skeleton.**

---

## 4. Step 2 — slot definition

### 4.1 Two required anchors

| anchor | purpose | status |
|---|---|---|
| VirB6 = `PF04610` | downstream slot boundary | validated, GA |
| VirB5 / TrbJ = `PF07996` | deciding that a slot is EMPTY | **resolved, but requires E ≤ 1e-5 — GA misses both IncP controls (§1.1a Finding 2)** |

**The second is critical and easy to miss.** The clean definition of "slot empty"
is "the gene immediately upstream of VirB6 *is* VirB5/TrbJ". Without a VirB5
anchor you cannot distinguish **an empty slot** from **a slot occupied by
something unrecognised** — and the latter is exactly the target.

### 4.2 Formalising "immediately upstream": compute both

| definition | |
|---|---|
| D1 | Pyrodigal gene index −1, relative to VirB6's own strand orientation |
| D2 | nearest same-strand gene with gap < X bp |

The two agree on all four controls (gaps 6–15 bp). **They will not agree
library-wide.**

**Decision: compute both and record the disagreement rate.** Low → adopt D1
(simpler). High → the disagreement rate is itself a reportable observation
(skeleton rearrangement frequency).

### 4.3 Gap threshold: do not fix it yet

The four controls give 6–15 bp, very tight, but n=4. **Record the distribution
first, then set a threshold.**

### 4.4 Full virB1–11 ordering: demoted to secondary analysis

Retained for two purposes: (a) detecting skeleton rearrangement, (b) extending to
other slots later. **Not on the main line.**

If done, the ordering must be fault-tolerant: Guglielmini 2013 notes that in
MPF_T, virB5/virB6 sometimes sit after virB10. **Name slots by the locally
observed adjacent pair (X|Y), never by forcing a global VirB order**, or
rearranged plasmids get mis-scored as "containing an insertion".

---

## 5. Step 3 — capture the slot occupant

For each VirB4⁺ plasmid:

```
1. locate the VirB6 anchor (PF04610, GA)
2. take the gene immediately upstream of VirB6 (compute D1 and D2 separately)
3. record:
     aa length
     strand / same strand as VirB6?
     gap to VirB6 (bp)
     start codon
     lipobox (precursor length, cleavage site, mature length)
     all HMM hits (including exclusion families -- POST-HOC CLASSIFICATION ONLY)
     is it VirB5?  (-> slot is empty)
```

### 5.1 The lipobox is worth using as a secondary feature

All four controls carry one, with cleavage sites at different positions
(C15/C17/C22/C23). **TrbK and Eex sit at background identity yet share an
architecture**: MKK/MKN start + hydrophobic stretch + Cys at 15–23. That is a
feature homology cannot detect, usable as independent confirmation of a candidate.

⚠️ **UniProt's Signal annotation is systematically wrong on lipoproteins:**
RP4 TrbK off by one (Signal 1..21, Cys at 23), ICEBs1 YddJ off by eight, pKM101
Eex has no signal annotation at all. **Compute the lipobox yourself; do not use
UniProt's feature.**

---

## 6. Step 4 — report three numbers first

**Before hunting new candidates, answer whether the criterion carries signal.**

| # | number | reads as |
|---|---|---|
| 1 | plasmids that are VirB4⁺ and CONJscan-MPF_T (one count per record, one per cluster) | set size |
| 2 | **fraction of upstream slots that are occupied** | ← decides whether the criterion has signal |
| 3 | **fraction of occupants hit by no existing Eex HMM**, split by lipobox | ← candidate pool, **not** a de novo count (see §6.3) |

### 6.1 Decision thresholds

**Number 2:**

| value | conclusion |
|---|---|
| near 100% | slot occupancy is not a variable; presence/absence carries no signal → change the criterion (sequence divergence rate, or within-slot family diversity) |
| 60–90% | signal present, continue |
| very low | check the anchor and strand logic for a bug; do not treat as a result |

**Number 3:**

| value | conclusion |
|---|---|
| near 0 | the eight existing Eex HMMs already cover this slot → a new slot (the traN one) is required for anything novel |
| clearly > 0 | that set is the de novo candidate list; proceed to downstream validation |

### 6.3 Number 3 must be reported split, and is not a candidate count

**Do not report number 3 as "de novo candidates".** §1.1a Finding 3 found VirB7
invisible to Pfam — pKM101 `traN` is 48 aa and hits nothing in the whole
database. Generalise it: **any 40–80 aa T4SS component may be invisible to Pfam**,
and the slot occupants are 69–76 aa, squarely in that band.

So an occupant that no Eex HMM hits is not necessarily a novel exclusion gene. It
may be an unannotated *structural* component. That is the dominant contamination
of number 3 and there is no anchor available to remove it — VirB7 has no family
to search with.

Mitigation is reporting, not filtering. Split number 3 by:

| stratum | interpretation |
|---|---|
| lipobox present | enriched for exclusion genes — all four controls have one (`LAG↓C23`, `VAG↓C22`, `LVA↓C15`, `LSA↓C17`) |
| lipobox absent | weaker; report separately |

⚠️ **The lipobox alone is not sufficient**, because VirB7 is also a lipoprotein.
Combine with length and strand to narrow it. State this when the number is
reported; the honest label is "candidate pool", and a de novo count only exists
after the structural-component contamination is assessed some other way.

### 6.2 The baseline that must be beaten

The old heuristic ("in-operon + VirB6 ± window") has **precision = 11%**
(53 in-window decoys; `docs/heuristic_precision.md`). The new method's precision
must clearly exceed that or there is no improvement.

Precision proxy: among occupants, the fraction hit by an existing Eex HMM serves
as an approximate true-positive rate, since those are near-certainly exclusion
genes.

---

## 7. The sfx benchmark (a scope decision — state it now)

The project's self-declared final benchmark is recovering IncC `sfx` from scratch
(TraDIS-discovered, no HMM family at all).

`sfx` lies between `traN` and `acaB`, and `traN` is not one of VirB1–11.

**Earlier text here said this put `sfx` "structurally outside coverage" so that
"no amount of statistical power will retrieve it". That was wrong, and the
correction matters more than the original claim.** It rested on the same mistake
as the retracted VirB7 blind spot (`docs/virb7_blindspot_RETRACTED.md`): assuming
that because a component has no VirB number, it has no family model. It does —
`PF06986` (F_T4SS_TraN) and `PF20898` (P_T4SS_TraN), both now in
`data/anchors/anchor_set.tsv`. Anchors are defined by Pfam family, never by VirB
number and never by IncP gene name, so a TraN anchor is expressible today.

**Decision: the MPF_T main line does not include the sfx benchmark — deferred,
not infeasible.** What it actually requires:

1. extension to MPF_F (`sfx` is IncC; the MPF_T line has no IncC systems)
2. widening the skeleton definition to CONJscan mandatory + accessory components,
   which `traN` already satisfies via `PF06986`

Both are work, neither is a barrier. Deferred because the MPF_T main line is not
finished, not because the target is unreachable.

**Do not restate this as a method limit.** The failure mode to avoid is the one
already observed twice: an absence produced by our own tooling being written down
as a property of the biology. `sfx` is out of scope for the current run and in
scope for the method.

Related: the entry-only scope removed the hardest targets — the TraDIS-class,
family-less cases are all surface exclusion. De novo targets are currently only
F `traS` and R27 `eexA`; **n=2 cannot evaluate a discovery method.** The primary
purpose of extending to surface exclusion is to bring the de novo target count to
a usable size, not taxonomic completeness.

---

## 8. Tools deleted from v1

Recorded explicitly so nobody executes v1 by mistake.

| tool | v1 role | now |
|---|---|---|
| Corekaburra | blueprint for steps 4–5 | **deleted**; no pan-genome core definition needed |
| PPanGGOLiN / panRGP | same | **deleted**; RGP works at genomic-island scale, irrelevant to 69 aa, and needs ≥15 genomes |
| MCScanX / i-ADHoRe / SibeliaZ | synteny blocks | excluded earlier |
| six-frame all-starts | primary discovery layer | **demoted to conditional supplement** (§2.2) |
| StORF-Reporter | overlap recovery | demoted. Note it **cannot see nested genes in principle** — mob4 sits inside annotated mob3, so that region is never a UR |
| full virB1–11 ordering | step 3 | demoted to secondary analysis (§4.4) |
| smORFinder | small-protein discovery | not used; biased to ≤50 aa and inherits Prodigal's misses |
| Bakta sORF module | same | not used; only <30 aa, and discards anything without a homology hit — the opposite of the false-negative priority |
| Balrog | caller | not used; 11% fewer extra predictions than Prodigal (wrong direction) and slower |
| clinker / pyGenomeViz | visualisation | **kept**, use directly |

---

## 9. Three things to fix before running

| # | open item | why it blocks |
|---|---|---|
| 1 | **Mash clustering threshold** (95% ANI suggested) | sets the statistical denominator; changing it re-runs every statistic |
| 2 | ~~VirB5/TrbJ Pfam accession~~ | **DONE** — `PF07996` at E ≤ 1e-5 (§1.1a) |
| 3 | ~~empirical anchor set~~ | **DONE** — `data/anchors/anchor_set.tsv` (§1.1a) |

Items 2 and 3 are resolved. **Only item 1, the Mash clustering threshold, still
blocks.** It sets the statistical denominator and changing it re-runs everything,
so it must be fixed before step 1 is launched, not after.

---

## 10. Sentinel rule (already fixed policy)

**Every benchmark must include a baseline expected to score perfectly. If it does
not, stop and fix before reading any other result.**

Origin: the GenBank baseline failed its own annotation in the ORF-caller
benchmark, and what it caught was a test-set bug (R621a excA is `−` strand and
the test set said `+`), not a caller bug. That was the **first of five silent
failures caught by design rather than after the fact.**

---

## Appendix: outstanding provenance items (not blocking this pipeline)

| item | status |
|---|---|
| origin of the F TraG 610–673 window | unresolved. The 17% in the same protocol sentence is disproven (it is TraS-vs-TraS), so the whole sentence is unreliable. Empirical null gives p = 0.022, with traB at 6.75× vs TraG 6.88× — essentially tied → **not significant, cannot be a primary claim** |
| F TraS aa_len | DISPUTED. UniProt conflict 144–173 / 1987 record 149 aa / AP001918 173 aa / audette2007 gives no length |
| IncC eexC coordinates | NOT_RECOVERED; needs humbert2019 full text (PMC6482922 reCAPTCHA; ASM 403) |
| R27 eexA / trhG | E3 (cross-Inc-group homology from pAPEC-O1-R, 78.1%/76.1%, plus positional corroboration); needs gunton2008 |
| pKPC 151 vs 152 aa | closed. Only one in-frame start at that stop (152 aa ATG); 151 = 152 minus the initiator Met |
| CP017937.1 | classified `accession_related` (gene calls differ from CP009465.1); **scripts must never use it for coordinate validation** |

### Domain-noise warning (six independent instances)

The literature carries systematic **attribute-mismatch** noise. Verify the
attribution of any number against the original before adopting it.

1. 17% (TraS-vs-TraS) welded onto a TraG window
2. kamruzzaman2022 states RP4's lipobox + Asp@+2 as a property of pKPC (no such
   ORF exists in the sequence)
3. NCBIfam family name `surf_exc_IncI1` hits ExcA, which is entry exclusion
4. Furuya & Komano 1994's title says "surface exclusion"; the gene is entry
   exclusion
5. UniProt `gene:eexA`/`eexB` → RP4 TrbJ/TrbK, colliding with R27's eexA/eexB
6. within one sentence, 69 aa includes the initiator Met and 151 aa does not

**Consequences:** gene names cannot be join keys (must be
`(element_id, gene_name)`); family names cannot imply mechanism; paper-title
keywords cannot drive automatic classification.

---

## Appendix B: control roster for the production run (added 2026-09-03)

| control | Inc | in PLSDB? | role in step 1a |
|---|---|---|---|
| RP4 | IncP-α | **no** (TPA record `BN000925.1`) | injected into the sketch explicitly |
| R751 | IncP-β | yes, `NC_001735.4` | native |
| R388 | IncW | yes, `NC_028464.1` | native |
| pKM101 | IncN | **no** — `U09868.1` is a 12.4 kb *tra* fragment | **structurally cannot enter** |
| **R46** | **IncN** | **yes, `NC_003292.1`** (50,969 bp) | **substitute for pKM101** |

pKM101 has no complete plasmid sequence, so it cannot be sketched and cannot be
an input to step 1a. Without a substitute, IncN — one of the two exclusion
families — would have no validated representative in the occupancy statistics.

**R46 verified as the substitute** (`AY046276.1` = `NC_003292.1`, 50,969 bp,
*S.* Typhimurium). Running the slot logic on it blind reproduces the architecture
exactly:

```
traC  7471..8184  237 aa  VirB5 (PF07996, E<=1e-5, aa>150)
eex   8192..8419   75 aa  <- SLOT, lipobox LVA|C15, mature 61 aa, "entry exclusion protein"
traD  8435..9475  346 aa  VirB6 (PF04610 GA, score 145)
```

Gaps 7 bp and 15 bp, matching the other four. The occupant hits `Eex_IncN`
(NF033894) at 70.5 — post-hoc classification only.

**R46's value is PROCESS validation, not biological evidence — record the two
separately.** The slot logic was previously demonstrated on `U09868.1`, a 12.4 kb
*tra* fragment with a local coordinate frame. R46 shows it works blind on a
**50,969 bp complete plasmid** — which is the kind of input step 1a actually
receives. That is a real result about the pipeline and it is not diminished by
what follows.

**But R46 is NOT an independent observation.** pKM101 is a deletion derivative of
R46 and their `eex` proteins are byte-identical (`fd988b8cd6b2e4c3`, 75 aa,
`AAA86454.1` vs `AAL13388.1`). R46 restores IncN's structural presence in the
production run; it adds no new evidence for the family. The family-level claim
remains **n=2** (`TrbK_RP4` vs `Eex_IncN`), as `docs/clustering_is_deduplication.md`
states.
