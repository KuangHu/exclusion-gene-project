# ORF-caller benchmark — results

10 target genes · 9 methods · 4 tests. Coordinates from frozen records; the
harness is `testset.py`, `sixframe.py`, `methods.py`, `test1..test4`.

**Three adapter bugs had to be fixed before any number here was meaningful**, and
the first run was measuring them rather than the callers:
* R621a `excA` was recorded as `+` strand; it is `−`. The GenBank baseline
  correctly failed its own annotation, which is what exposed it.
* orfipy's BED **excludes the stop codon** (ColE1 mob2: orfipy 1842..2187 len 345
  vs GenBank 1843..2190, 348 nt). Uncorrected this gave orfipy 0/10.
* StORF-Reporter needs `-anno Pyrodigal Single_FASTA` (two values) and emits a
  **combined** pyrodigal+StORF GFF with **stop-to-stop** coordinates.

---

## TEST 1 — recall

| method | 3′ recall | 5′ exact | verdict |
|---|---|---|---|
| 1 GenBank baseline | **10/10** | 10/10 | — |
| 2 Pyrodigal anon default | 9/10 | 8/10 | fails go/no-go |
| 3 Pyrodigal anon min_gene=60 | 9/10 | 8/10 | fails go/no-go |
| 4 Pyrodigal normal (trained) | 7/10 + 2 errors | 6/10 | worst |
| 5 orfipy | **10/10** | 3/10 | passes 3′ |
| 6 six-frame all-starts | **10/10** | **10/10** | **passes both** |
| 6b six-frame first-start | **10/10** | 3/10 | passes 3′ |
| 7 StORF default | 9/10 | 7/10 | fails go/no-go |
| 8 StORF con_storfs | 9/10 | 7/10 | fails go/no-go |

**Only method 6 achieves 10/10 on both.** Retaining every in-frame start is what
makes the 5′ column full: the true start is always among the candidates.

### Where the design's predictions did not hold

* **Prodigal did NOT miss the short genes.** `trbK` (210 nt) and `eex` (228 nt)
  were both found, with **exactly correct starts**, by methods 2 and 3. The
  length-penalty hypothesis predicted misses here and it did not happen on these
  records.
* **The single miss for every statistical caller is ColE1 `mob4`** — 186 nt,
  fully nested in a different frame inside `mob3`. The failure mode is
  **overlap, not length**.
* **StORF also misses `mob4`**, and structurally so: it searches only
  *unannotated* regions, and `mob4` lies inside pyrodigal-annotated `mob3`, so
  that region is not a UR at all. `con_storfs` does not help.
* **Method 4 was predicted worst and is.** Its two errors are explicit:
  `sequence must be at least 20000 characters` on U09868.1 (12,416 bp) and
  J01566.1 (6,646 bp). The training-length limit is a hard floor, not a
  degradation.

---

## TEST 2 — false-positive cost per slot window

Windows defined by the flanking Pyrodigal skeleton calls, as the architecture
would: **829–2,354 bp**.

| method | decoys @40 aa | @60 aa | @80 aa |
|---|---|---|---|
| Pyrodigal anon | 1 med / 3 max | 1 / 3 | 0 / 3 |
| StORF | 2 / 3 | 2 / 3 | 2 / 3 |
| orfipy | 3 / 8 | 3 / 8 | 1 / 6 |
| six-frame first-start | 3 / 9 | 3 / 9 | 1 / 6 |
| **six-frame all-starts** | **7 / 39** | **7 / 39** | 2 / 34 |

The geometric estimate (1–3 per slot at a 60 aa floor with ATG/GTG/TTG) **matches
the one-ORF-per-stop enumerators exactly** — median 3. All-starts costs about
2.3× that, median 7, which is the price of the 5′ column being full.

**GC dependence is real: Pearson r(GC, decoys/kb) = 0.613.** Range 1.74/kb
(R100, 39.8% GC) to 16.57/kb (ColE1, 53.6% GC). RP4 at 63.1% GC gives 8.44/kb.
The ColE1 window is also the widest (2,354 bp) and is a dense mobilization
region, so it inflates the correlation; the trend holds without it but is weaker.

Every non-target candidate in these windows is a negative-set member.

---

## TEST 3 — overlap retention

| method | R64 excA/excB (same stop) | ColE1 mob3/mob4 (nested, diff frame) |
|---|---|---|
| 1 GenBank | **PASS** | **PASS** |
| 2/3/4 Pyrodigal | excA exact, excB 3′-only | **FAIL** — mob4 not proposed at all |
| 5 orfipy | excA exact, excB 3′-only | 3′-only for both |
| **6 six-frame all-starts** | **PASS** | **PASS** |
| 6b six-frame first-start | excA exact, excB 3′-only | 3′-only for both |
| 7/8 StORF | excA exact, excB 3′-only | **FAIL** |

**The prediction that only methods 5–8 would pass is wrong.** Only method 6
passes. Methods 5, 6b, 7 and 8 are all one-ORF-per-stop, so on the same-stop pair
they emit `excA` and reach `excB`'s stop but never its start — 3′-only, not a
pass. This is a model-expressiveness limit, not a threshold setting: no amount of
parameter tuning makes a one-per-stop caller emit two proteins from one stop.

---

## TEST 4 — reproducibility on identical sequence

`CP009465.1` vs `CP017937.1` — byte-identical, 43,621 bp both.

| | result |
|---|---|
| **1 GenBank annotation** | **DIFFERS: +3 / −3 of 47 CDS (13%)** |
| all 8 computational methods | **IDENTICAL** |

The differing calls:

```
only in CP009465.1                          only in CP017937.1
 25984..26175 +  63aa AJE44587.1              43..279   -  78aa APB53843.1
 26859..27209 + 116aa AJE44589.1 trbM family  27251..27463 + 0aa
 43209..43466 -  85aa AJE44601.1              39614..40513 - 0aa Tn3 transposase
```

**The caller is not what manufactures presence/absence signal — the deposited
annotation is.** Every caller tested, statistical ones included, is deterministic
on identical input. `AJE44589.1` is a *trbM family protein* — a real conjugation
gene present in one deposit and absent from the other. `AJE44587.1` is the 63-aa
lipobox ORF previously considered as a pKPC `trbK-like` candidate.

This inverts the test's rationale in a useful direction: the argument for
re-calling ORFs uniformly is not that callers are unreliable, but that
**deposited annotations are**, and a presence/absence study built on GenBank CDS
features inherits ~13% annotation-pipeline noise that looks like biology.

---

## Conclusion

The architectural recommendation holds, but for a different reason than proposed.

* Not length. Prodigal handled 210 and 228 nt genes correctly, starts included.
* **Overlap.** Every non-enumerative method fails on a nested different-frame
  gene, and every one-per-stop method fails on a same-stop second start.
* Uniform re-calling is justified by TEST 4, not by caller unreliability.

**Recommended layering, supported by these numbers:**
Pyrodigal for the skeleton and slot definition (9/10 3′, 8/10 5′, ~1 decoy per
window), plus six-frame all-starts **inside slot windows only** to recover
overlaps and multi-start cases (10/10 both, ~7 decoys per window of 829–2354 bp).
Running all-starts genome-wide is unnecessary; running it inside slots costs a
median of 7 candidates and is the only configuration that expresses excA+excB.

---

# Round 2 — corrections and the hard-case set

## TEST 6 — hard-case recall. This replaces Test 1 as the headline number.

Test 1's targets are all GenBank-annotated CDS, and GenBank annotation comes from
PGAP / GeneMarkS-2+ — statistical callers. **The test set was pre-filtered to
genes a statistical caller can already find.** Pyrodigal's 9/10 there is an
optimistic upper bound, not an unbiased estimate.

`CP009465.1` and `CP017937.1` are the same sequence annotated twice by production
pipelines that disagree. Every gene in one and not the other was **missed by a
real annotation pipeline on the exact sequence the other pipeline saw** — a hard
set with no circularity, since the label comes from one pipeline and the
difficulty is demonstrated by another. (Two of the six differing entries are
`/pseudo`, "incomplete; partial on complete genome; missing stop", and are
excluded: those are annotation-status differences, not gene-call differences.
The real difference is **4 of 47 CDS ≈ 8.5%**, not the 13% reported earlier.)

| method | AJE44587 (63aa) | AJE44589 (116aa, trbM) | AJE44601 (85aa) | APB53843 (78aa) | 3′ |
|---|---|---|---|---|---|
| GenBank CP009465 | 5′ | 5′ | 5′ | **MISS** | 3/4 |
| GenBank CP017937 | **MISS** | **MISS** | **MISS** | 5′ | 1/4 |
| Pyrodigal anon (def / m60) | 5′ | 5′ | **MISS** | **MISS** | **2/4** |
| Pyrodigal normal | 5′ | **MISS** | 5′ | 3′ | 3/4 |
| orfipy | 5′ | 5′ | 3′ | 3′ | 4/4 |
| **six-frame all-starts** | **5′** | **5′** | **5′** | **5′** | **4/4** |
| six-frame first-start | 5′ | 5′ | 3′ | 3′ | 4/4 |
| StORF (default / con) | 5′ | 5′ | 5′ | **MISS** | 3/4 |

**Pyrodigal falls from 9/10 (90%) on the pre-filtered set to 2/4 (50%) here.**
Six-frame all-starts holds at 100% and is the only method with exact starts on
all four.

This materially changes the recommendation. On the pre-filtered set, Pyrodigal
looked adequate as a sole primary layer and six-frame was needed only for
overlap. On the unbiased set, Pyrodigal's recall is not adequate on its own —
which supports the original recall-based argument after all, on evidence Test 1
could not provide.

**Caveat: n = 4.** 2/4 vs 4/4 is directionally clear but not statistically
strong. It is the only near-unbiased set currently available, and it should be
expanded — every duplicate-deposit pair in PLSDB is a free hard-case generator.

## TEST 5 — start ambiguity, resolved

pKPC `trbK-like`: **only one in-frame start ≥40 aa exists at that stop**, the
152-aa ATG (`AJE44581.1`, 21029..21487). All eight methods agree; a 151-aa form
would need a start 3 nt downstream and there is none.

So the 151-vs-152 discrepancy is **not** start ambiguity. Given the paper's own
loose arithmetic ("twice" 69 aa = 138 ≠ 151), 151 is almost certainly 152 counted
without the initiator methionine. No sequence question remains.

## GC correlation — withdrawn

r = +0.613 (t = 2.20, n = 10). Removing ColE1: r = +0.622 (t = 2.10, n = 9).
**Not significant either way, and not single-point-driven — simply underpowered.**
Do not set GC-stratified thresholds on this.

## Two-stage architecture — FP inflation leaves the statistics

The 7-median / 39-max all-starts cost is misleading because **those candidates
are not independent hypotheses**: every extra start shares a stop codon with one
already counted. 39 candidates is roughly 8 stop-groups × ~5 starts.

**Stage 1 — is this slot occupied?**
Unit = the stop-to-stop locus, counted one-per-stop. FP ≈ 3 per window.
**Cross-plasmid enrichment testing happens only here.**

**Stage 2 — given occupancy, what is the protein?**
Enumerate all in-frame starts, but only within loci that passed stage 1.

The 2.3× all-starts cost then never enters the statistical test; it applies only
to positions already judged occupied, so the enrichment denominator is unaffected.
R64 excA/excB falls out naturally: for occupancy it is **one** occupied slot, and
the two products are a stage-2 annotation question.

## Two engineering constraints for the schema

* **Pyrodigal normal mode has a hard floor, not graceful degradation:**
  `sequence must be at least 20000 characters`. U09868.1 (12.4 kb) and J01566.1
  (6.6 kb) raise. Fragment records and small plasmids **must** use anon mode.
  This needs an assertion or a batch run will fail silently on a subset.
* **Every new caller adapter must first reproduce one known CDS exactly** before
  use. orfipy's stop-codon-excluding BED (off-by-3) nearly invalidated the whole
  benchmark.

## The canary worked

The GenBank baseline failing its own annotation is what exposed the R621a strand
error — a test-set bug, not a caller bug. Unlike the four earlier silent
failures, this one was caught by a sentinel that was in the design.

**Rule: every benchmark must include a baseline that is expected to score
perfectly. If it does not, stop and fix before reading any other result.**
