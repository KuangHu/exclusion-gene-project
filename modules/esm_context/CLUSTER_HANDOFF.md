# Cluster handoff: PLSDB / IMGPR short-protein extraction

Start from this module's README. The first operational task is **CPU extraction**;
do not submit a full GPU run before inspecting the extraction manifest and input
volume. All filesystem paths, account and partition names are site-specific.

1. Locate the local PLSDB / IMGPR files. Prefer already available full CDS caches
   when their format and provenance match; otherwise use raw nucleotide FASTA and
   the repository's unchanged frozen ORF caller. No new database download is needed.
2. Use the existing Python 3.11 conda environment, install `requirements.txt`, and
   run the 12 tests. The adapter coordinate sentinel must pass.
3. Write a tab-delimited manifest with `dataset`, `path`, `format`. Never supply a
   protein file that has already been filtered by size/function as a full CDS set.
4. On a compute node, run a hash-selected smoke-test shard. Inspect counts of
   records, centers, unique required proteins and complete ±3 windows. Empty input,
   duplicate conflicting accessions and changed resume parameters are hard errors.
5. Submit the full CPU extraction array with identical settings. Reuse each shard's
   own directory on interruption. Merge only a complete shard set, or process
   extraction shards separately if their combined ESM input is too large for RAM.
6. GPU inference uses the local ESM2-3B checkpoint and the bundled reference/SDPA
   checks. It writes full-sequence vectors and validation files. No silent truncation.
7. Scoring is CPU-only and retains center, context, and fusion outputs separately.
   Incomplete contexts receive no context/fusion score. Functional labels remain
   unassigned; annotation/architecture review and anchor-qualified shortlisting
   are separate later tasks.

The packaged small heads are a full-cohort exploratory refit. They are distinct
from the earlier evaluation-fold heads. The validation figures in the README
describe the prior procedure, not independent validation of this refit. Numerical
reference proteins are compatibility checks, not positive functional sentinels.

Return the extraction/scoring manifests and summary counts first. Keep full
genomes, intermediate databases, and large embedding matrices on the cluster.
