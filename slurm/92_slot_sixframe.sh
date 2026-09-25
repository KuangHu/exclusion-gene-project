#!/bin/bash
#SBATCH --job-name=slot6f
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --array=0-15
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/92_slot6f_%A_%a.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/92_slot6f_%A_%a.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
P=/global/home/users/kh36969/exclusion_gene_project
python3 "$P/scripts/lib/orf_caller.py" >/dev/null || { echo "caller sentinel FAILED"; exit 1; }
python3 "$P/scripts/92_slot_sixframe_recheck.py" \
    --src "${SRC:-$P/data/positional}" --tag "${TAG:-v2}" \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}"
