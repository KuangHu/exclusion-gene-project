#!/usr/bin/env bash
#SBATCH --job-name=short-context
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
set -euo pipefail
: "${INPUT_MANIFEST:?Set INPUT_MANIFEST to the input TSV}"
: "${EXTRACT_ROOT:?Set EXTRACT_ROOT to the output root}"
MODULE_DIR="${EEX_MODULE_DIR:-${SLURM_SUBMIT_DIR:?Submit from the module directory}}"
PYTHON="${EEX_PYTHON:-python}"
INDEX="${SLURM_ARRAY_TASK_ID:-0}"
SHARDS="${EEX_NUM_SHARDS:-1}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec "$PYTHON" "$MODULE_DIR/extract.py" --inputs "$INPUT_MANIFEST" \
  --out "$EXTRACT_ROOT/shard_$INDEX" --num-shards "$SHARDS" --shard-index "$INDEX" \
  --min-aa "${EEX_MIN_AA:-40}" --max-aa "${EEX_MAX_AA:-300}" "$@"
