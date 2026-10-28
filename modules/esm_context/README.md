# Short proteins from PLSDB / IMGPR: center and ±3 ESM scores

This independent module extracts short-protein centers from whole plasmid records,
preserves **±7 CDS** for audit, and supplies **center + six ±3 neighbors** to two
separate ESM-based classifier heads. It does not modify the curated catalogues,
anchor criteria, or source annotations in the parent repository.

The default center range is **40–300 aa**, configurable. All CDS, including
partial or invalid translations, remain in the neighborhood order. Missing slots
are never compressed or filled with zeros. Offsets follow the center gene's
orientation. Circular wrap-around is disabled.

## Cluster setup

Use Python 3.11. An existing conda environment is fine; a system `python3.11`
binary is not required. From this module directory:

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -p 'test_*.py' -v
```

For GPU inference use a separate environment with the cluster-supported
CUDA-enabled PyTorch, then `pip install -r requirements-gpu.txt`. Stage the
official `facebook/esm2_t36_3B_UR50D` checkpoint on disk. No runtime model download
occurs. Do not mix ESM C, other ESM2 sizes, different layers, or truncated pooling.

## 1. Specify local input files

Copy `inputs.example.tsv` to a working location and replace its paths. Remove any
dataset row you are not using. Columns are `dataset`, `path`, `format`, separated
by **tabs**. Relative paths resolve against the manifest's directory. Dataset
names distinguish PLSDB and IMGPR records; accessions are retained verbatim.

Supported explicit formats:

| Format | Input |
|---|---|
| `genome-fasta` | Complete nucleotide records in `.fna` / `.fasta`, optionally `.gz` |
| `project-cds-faa` | Existing full CDS cache from `scripts/109_full_cds_cache.py`: `accession\|ordinal\|start1\|end1\|strand\|aa_length` |
| `prodigal-faa` | Full per-contig protein set with `# start # end # strand # ...partial=...` headers |
| `genbank` | CDS features with existing translations; source gene/product fields retained |

Protein FASTAs must contain the **full ordered CDS set**, not just selected short
proteins. Filtered/gapped ordinals and non-contiguous reappearance of a contig are
rejected. GenBank compound CDS are retained but not used as vector slots. Existing
annotations and frozen ORF calls can differ; their source formats are recorded.

For raw nucleotide input, the module imports the parent repository's frozen
`scripts/lib/orf_caller.py` unchanged (Pyrodigal meta, min_gene=90 nt,
closed=True, translation table 11). The repository adapter's exact-coordinate
golden test must pass before extraction. This sentinel does not set model thresholds.
Closed linear calling can miss origin-spanning or nested CDS; absence is not
evidence of absent function.

## 2. Extract on a CPU compute node

One complete run:

```bash
python extract.py --inputs /path/to/inputs.tsv --out /path/to/work/extracted
```

A deterministic 1/100 **smoke-test shard**, not the full dataset:

```bash
python extract.py --inputs /path/to/inputs.tsv --out /path/to/work/pilot/shard_0 \
  --num-shards 100 --shard-index 0
```

For a full 16-shard Slurm array, submit **from this module directory** (or set
`EEX_MODULE_DIR` to its absolute path):

```bash
export EEX_PYTHON=/path/to/cpu_env/bin/python
export INPUT_MANIFEST=/path/to/inputs.tsv
export EXTRACT_ROOT=/path/to/work/extraction_shards
export EEX_NUM_SHARDS=16
sbatch --account=YOUR_ACCOUNT --partition=YOUR_CPU_PARTITION --array=0-15 extract.sbatch
```

Sharding is SHA256(dataset + record ID), not file order. Each worker scans the
manifest's input streams but calls CDS only for its assigned records. For very
large IMGPR input, use separate manifests for record-complete input blocks to
reduce repeated shared-filesystem reads. Do not split in the middle of a record.

Each output has an exclusive writer lock and a transactional SQLite checkpoint.
Rerun the **same** command to resume. Changed file checksums, caller configuration,
length bounds, or shard parameters require a new output directory. Equal IDs with
conflicting content fail; duplicate identical source records are audited.

After every shard finishes, optionally merge that complete shard set:

```bash
python merge.py --shards-root /path/to/work/extraction_shards \
  --out /path/to/work/extracted
```

The merger requires exactly the expected shard indices, rejects conflicting
configurations, and deduplicates required protein sequences across shards.
An incomplete merge has no completion manifest; use a new directory to retry.
For large datasets, embedding/scoring individual extraction shards is also valid;
deduplicate by sequence/window hashes when combining statistics afterward.

### Extraction outputs

- `centers.tsv`: every eligible short center, original annotation, coordinates,
  full-sequence hash, and ±3/±7 completeness flags.
- `slots.tsv`, `windows.jsonl`: all ±7 positions, including missing/excluded CDS.
- `embedding_slots.tsv`: only the positions requested for ESM (default center ±3).
- `all_proteins.faa`: required valid sequences, deduplicated by full SHA256.
- `extraction.sqlite`, `extraction_manifest.json`: resume state and provenance.

`functional_label` stays `unassigned`. Center length is not a functional label.
Source occurrence counts are database occurrences, not independent biological
samples. Identical ordered ±7 sequence windows have a shared `window_sequence_id`.
No source records parsed is an error for an unsharded run. A legitimate empty
shard is recorded explicitly; skip GPU inference if there are no required sequences.

## 3. Calculate only the required ESM vectors on a GPU node

```bash
export ESM_PYTHON=/path/to/gpu_env/bin/python
export EXTRACTION=/path/to/work/extracted
export EMBED_OUT=/path/to/work/embeddings
export ESM_MODEL_DIR=/path/to/esm2_t36_3B_UR50D
export ESM_MODEL_FORMAT=esm2
sbatch --account=YOUR_ACCOUNT --partition=YOUR_GPU_PARTITION embed.sbatch
```

`esmfold-stem` is also supported for a compatible ESMFold checkpoint. The HF
ESM2 vocabulary is read from `vocab.txt` when `config.vocab_list` is null.
Inference is fp16 with float32 full-residue mean pooling, no BOS/EOS/padding,
no sequence truncation, and no chunk averaging. Overlength sequences are recorded
as extended context; a single-sequence OOM fails explicitly instead of truncating.

The inference path is adapted from the project's previously exercised cluster
runner. On each new checkpoint/run it compares eager and SDPA attention, and
checks three bundled **numerical reference vectors** against the classifier's
vector space. These references check arithmetic/checkpoint compatibility, not
functional recall or decision thresholds. Completion requires all vectors and
matrix/cache integrity checks. GPU inference was not rerun during CPU packaging;
the cluster must pass these numerical checks before scoring.

Rerun with the same output to resume its per-sequence cache. Never run two writers
in one directory. The runner keeps input sequence metadata in memory; size GPU
jobs/host RAM to the sequence count, or process extraction shards separately.

## 4. Keep two independent scores, then apply late fusion

```bash
python score.py --extraction /path/to/work/extracted \
  --embeddings /path/to/work/embeddings --out /path/to/results/scores
# Or set EXTRACTION, EMBED_OUT, SCORE_OUT, EEX_PYTHON and submit score.sbatch.
```

The scorer verifies input, matrix, model, and numerical-reference checksums.
`window_scores.tsv` retains the original gene/product and writes distinct
`center_*`, `context_*`, and `fusion_*` model fields. Incomplete ±3 windows can
receive a center score; context/fusion fields stay empty. A repeated exact copy
of the center in ±3 also disables context scoring.

`center_scores.tsv` averages scores over distinct full ±7 sequence-window variants
per center, with no database-copy weighting. Context/fusion averages use only
scorable variants; their denominator is explicit. Window-level outputs are also
retained so variant disagreement is not hidden. Exact matches to training center
hashes are flagged. These are not calibrated probabilities or validated functions.

The bundled small NumPy heads are **full-cohort exploratory refits**, not the
three validation-fold models. They were fitted on 21,421 centers / 26,436 windows:
728 known-label EEX, 741 conjugation comparators, 19,952 provisional GTDB backgrounds.
Training preserved family/source weights and full-length ESM2-3B vectors. Late
fusion uses standardized separate scores with fixed 50:50 weights. Normalization
and thresholds come from grouped out-of-fold calibration before the full refit.
See `model/model_card.json` and `model/refit_validation.json`.

The preceding paired internal comparison detected 725/728 EEX after fusion vs
710/728 with the center alone, with conjugation positive calls 8/741 vs 1/741 and
GTDB positive calls 230/19,952 vs 290/19,952. Those are validation-procedure results,
**not an independent test of the bundled all-data refit**. In particular 77/78 TraS
were in one strict group. PLSDB/IMGPR population shifts, caller differences, lengths
outside the training range, and newly included ±3-only windows require review.

Model calls never overwrite source annotation. This module does not create a
novelty shortlist: that requires a separate existing class-matched conjugation
anchor within ±7 (excluding the candidate), plus reference and architecture review.

## Module provenance

Embedding/cache helpers were adapted from `KuangHu/DL_entry_exclusion_gene` and
retain the tested attention and pooling path. The null-vocabulary fallback and
classifier numerical-reference gate are included here. Curated parent `config/`
and release catalogues are untouched. No genome datasets, full embedding matrices,
CUDA model checkpoints, local environments, or cluster credentials are in this module.
