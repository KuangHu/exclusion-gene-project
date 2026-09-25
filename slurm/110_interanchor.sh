#!/bin/bash
#SBATCH --job-name=ia6f
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --array=0-15
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/ia6f_%A_%a.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/ia6f_%A_%a.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
P=/global/home/users/kh36969/exclusion_gene_project
python3 -u "$P/scripts/110_sixframe_interanchor.py" --cpus 8 \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}"
