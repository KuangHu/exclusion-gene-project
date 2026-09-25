# The benchmark has to be redefined: 7 of 9 exclusion genes are ALREADY findable by HMM

`scripts/73_hmm_benchmark.py` · models in `/global/scratch/users/kh36969/exclusion_gene/hmm/`
· raw hits `data/level2_uniprot/hmm_seed_hits.domtbl`

## What was tested

Eight HMMs, all pre-existing public families, run with **Pfam/NCBIfam GA
thresholds** (`hmmsearch --cut_ga`) against the 20 frozen seed proteins:

| model | accession | description | source |
|---|---|---|---|
| TraS | `PF10624.15` | **"Plasmid conjugative transfer entry exclusion protein TraS"** | local Pfam-A |
| DUF4467 | `PF14729.13` | ICEBs1 YddJ fold | local Pfam-A |
| TrbL | `PF04610.20` | TrbL/VirB6 conjugal transfer | local Pfam-A |
| TraG_N | `PF07916.17` | TraG-like, N-terminal | local Pfam-A |
| EexR | `NF041429.1` | **"EexR/EexS family entry exclusion protein"** | NCBI FTP |
| Eex_IncN | `NF033894.1` | **"EexN family lipoprotein"** | NCBI FTP |
| surf_exc_IncI1 | `NF033891.1` | "plasmid IncI1-type surface exclusion" | NCBI FTP |
| TrbK_RP4 | `TIGR04359.1` | **"entry exclusion lipoprotein TrbK"** | NCBI FTP |

## Result: hmm_recoverable

**Exclusion side — 7 of 9 with sequence are recovered:**

| element | gene | model | score |
|---|---|---|---|
| R100 | traS | `TraS` PF10624 | 312.5 |
| R64 | excA | `surf_exc_IncI1` NF033891 | 289.7 |
| SXT | eex | `EexR` NF041429 | 274.5 |
| R621a | excA | `surf_exc_IncI1` NF033891 | 216.4 |
| ICEBs1 | yddJ | `DUF4467` PF14729 | 98.9 |
| RP4 | trbK | `TrbK_RP4` TIGR04359 | 75.6 |
| pKM101 | eex | `Eex_IncN` NF033894 | 70.5 |
| **F** | **traS** | **no hit** | — |
| **R27** | **eexA** | **no hit** | — |

**Partner side — 8 of 11.** `TraG_N` covers F/R100/IncC/SXT/R27; `TrbL` covers
RP4/pKM101/pKPC. Misses: R64 and R621a TraY (MPF_I, needs DotA_TraY) and ICEBs1
ConG (MPF_FA). This empirically confirms the four-anchor-by-MPF-class split.

## What this does to the project's premise

The stated rationale — that exclusion genes cannot be found by compositional or
homology features and require a relational signal — is **substantially weakened**.
Most of the seed is already covered by public HMMs, several of them purpose-built
for entry exclusion.

**Success is therefore NOT "recover the seed".** Recovering the seed largely
re-finds what `hmmsearch --cut_ga` finds in minutes. The correct criterion:

> The synteny pipeline's value is what it finds OUTSIDE `hmm_recoverable`.

Consequences:

* Every seed row needs an `hmm_hit` column. A row hit by `TIGR04359` is a
  **positive control**, not a de novo discovery target.
* The only genuine de novo targets in the current seed are **F traS** and
  **R27 eexA** — 2 of 9, not 9 of 9.
* F traS failing `PF10624` while R100 traS scores 312.5 is consistent with their
  measured 17.7% identity: the family is built around the IncFII type. F is a
  real HMM blind spot, and the best single de novo test case available.
* The real benchmark targets are TraDIS-class discoveries with no family at all
  (IncC `sfx`) — but `sfx` is surface exclusion and out of scope under the
  entry-only cut. **The entry-only scope removed the hardest cases.** That
  tension needs resolving before the scan is evaluated.

## Naming trap, fifth instance — now baked into a family name

`NF033891` is named **`surf_exc_IncI1`** ("surface exclusion") but its hits are
R64 and R621a **ExcA**, which are **entry** exclusion (`furuya1994`'s title makes
the same error). Any downstream step that infers mechanism from a family name
inherits it. Mechanism must come from the element record, never from a model name.

## Not yet done

`hmm_recoverable` over all of PLSDB (72,556 plasmids) is the denominator this
argument ultimately needs. It requires prodigal over PLSDB then hmmsearch —
a SLURM array job, not yet written. The seed-level result above already forces
the redefinition; the PLSDB run quantifies how much.
