#!/bin/bash
#SBATCH --job-name=up6
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --array=0-7
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/95_up6_%A_%a.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/95_up6_%A_%a.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
P=/global/home/users/kh36969/exclusion_gene_project
python3 "$P/scripts/95_upstream_of_virb6.py" --cpus 8 \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}"
