#!/bin/bash
#SBATCH --job-name=fullcache
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --array=0-31
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/fullcache_%A_%a.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/fullcache_%A_%a.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
P=/global/home/users/kh36969/exclusion_gene_project
python3 "$P/scripts/lib/orf_caller.py" >/dev/null || { echo "caller sentinel FAILED"; exit 1; }
python3 -u "$P/scripts/109_full_cds_cache.py" \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}"
