#!/usr/bin/env bash
#SBATCH --job-name=context-score
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
set -euo pipefail
: "${EXTRACTION:?Set EXTRACTION}"
: "${EMBED_OUT:?Set EMBED_OUT}"
: "${SCORE_OUT:?Set SCORE_OUT to a new output directory}"
MODULE_DIR="${EEX_MODULE_DIR:-${SLURM_SUBMIT_DIR:?Submit from the module directory}}"
PYTHON="${EEX_PYTHON:-python}"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
exec "$PYTHON" "$MODULE_DIR/score.py" --extraction "$EXTRACTION" --embeddings "$EMBED_OUT" --out "$SCORE_OUT" "$@"
