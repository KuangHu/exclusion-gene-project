#!/usr/bin/env bash
#SBATCH --job-name=context-esm
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
set -euo pipefail
: "${EXTRACTION:?Set EXTRACTION to a completed extraction or merged directory}"
: "${EMBED_OUT:?Set EMBED_OUT to a new or resumable embedding directory}"
: "${ESM_MODEL_DIR:?Set ESM_MODEL_DIR to the staged local checkpoint}"
MODULE_DIR="${EEX_MODULE_DIR:-${SLURM_SUBMIT_DIR:?Submit from the module directory}}"
PYTHON="${ESM_PYTHON:-python}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
exec "$PYTHON" "$MODULE_DIR/embed_neighborhoods.py" --fasta "$EXTRACTION/all_proteins.faa" \
  --slots "$EXTRACTION/embedding_slots.tsv" --out "$EMBED_OUT" \
  --model-dir "$ESM_MODEL_DIR" --model-format "${ESM_MODEL_FORMAT:-esm2}" \
  --batch-size "${ESM_BATCH_SIZE:-16}" --token-budget "${ESM_TOKEN_BUDGET:-4096}" "$@"
