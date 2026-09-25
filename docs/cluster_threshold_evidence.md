# ANI clustering threshold — the evidence, before choosing

Pipeline v2 §0.1 blocker #1. Job `25502424` (64 cores, 1h04), all-vs-all Mash
over PLSDB + RP4 injected. **56,192,194 pairs** at d ≤ 0.10.
Scripts: `slurm/81_mash_dist.sh`, `scripts/82_cluster_threshold.py`.

## 1. There is a valley, and it is nowhere near 95%

The distribution **rises monotonically** toward the d=0.10 cutoff, so 95% ANI
(d=0.05) sits mid-slope with no natural break — as anticipated. The real
structure is at the top end:

```
d 0.000-0.001   99.9-100%    1,028,344   <- oversampling spike
d 0.001-0.002   99.8-99.9      324,784
d 0.002-0.003   99.7-99.8      199,034   <- 5x below the spike
d 0.003-0.004   99.6-99.7      197,284
d 0.005-0.006   99.4-99.5      206,394
d 0.006-0.007   99.3-99.4      190,740   <- TROUGH MINIMUM
d 0.007-0.008   99.2-99.3      250,098
...rising...
d 0.029-0.030   97.0-97.1      426,520
```

The 1.03 M pairs above 99.9% ANI are the repeated-sequencing regime §0.1 warns
about. The trough separating that from everything else spans
**d ≈ 0.002–0.007 (99.3–99.8% ANI)**, minimum at d = 0.006–0.007.

## 2. The control-separation test does NOT discriminate

Proposed test: if RP4 (IncP-α) and R751 (IncP-β) merge, the threshold is too
loose — their TrbK proteins are only 37.8% identical, so they are independent
events.

**Result: all three testable controls stay in separate clusters at every
threshold from 99.9% down to 95%.**

| d | ANI | controls |
|---|---|---|
| 0.001 … 0.050 | 99.9 … 95.0 | separate (3/3) at every value |

Sound test, but it does not constrain the choice here — RP4 and R751 are far
enough apart that even 95% keeps them apart. (pKM101 is untestable: no complete
plasmid sequence, so it cannot be sketched. RP4 is a TPA record, absent from
PLSDB, and was injected explicitly.)

## 3. What DOES disqualify 95% is percolation

| d | ANI | clusters | largest | largest % of DB |
|---|---|---|---|---|
| 0.007 | 99.3 | 34,436 | 860 | 1.2% |
| 0.010 | 99.0 | 31,535 | 864 | 1.2% |
| 0.016 | 98.4 | 26,999 | 1,893 | 2.6% |
| **0.018** | **98.2** | 25,830 | **6,523** | **9.0%** ← 3.4× jump |
| 0.020 | 98.0 | 24,826 | 12,582 | 17.3% |
| 0.030 | 97.0 | 20,788 | 24,039 | 33.1% |
| **0.050** | **95.0** | 15,856 | **38,712** | **53.4%** |

**At 95% ANI a single cluster holds 38,712 of 72,557 plasmids — 53% of the
database as one "independent observation"**, with the top 10 clusters holding
56.8%. That is precisely the denominator failure §0.1 exists to prevent, and it
is far worse than the 3,000-member cluster the review posited.

The graph percolates between **d = 0.016 and d = 0.018 (ANI 98.2–98.4%)** — a
3.4× jump in one step. Any threshold looser than that is unusable regardless of
what the controls do.

## 4. Recommendation

**d = 0.005 (99.5% ANI)**, with d = 0.007 (99.3%) an equally defensible
alternative.

| criterion | at d=0.005 |
|---|---|
| inside the distribution trough | yes (206 k/bin vs 1.03 M at the spike) |
| below percolation | yes, by 3.6× margin (0.005 vs 0.018) |
| largest cluster | 847 = 1.2% of the DB |
| clusters | 36,767 from 72,557 records (2.0 records/cluster) |
| singletons | 29,830 |
| controls separate | yes |

Three independent lines agree on the top end of the range, and none supports 95%.

**Not adopted yet.** §0.1 says fixing this re-runs every statistic, so the number
is presented rather than chosen. Once set, record it here and in
`config/` before step 1a launches.

## 5. Reporting consequence

Occupancy must still be reported per record **and** per cluster (§0.1). At
d=0.005 the two differ by roughly 2× on average, but the ratio is very uneven —
29,830 singletons against a largest cluster of 847 — so cluster-weighted and
record-weighted occupancy can diverge sharply for any Inc group that happens to
be heavily resequenced. Report both, always, and never a single blended figure.
