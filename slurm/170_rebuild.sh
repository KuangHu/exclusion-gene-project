#!/bin/bash
#SBATCH --job-name=shardfix
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/shardfix_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/shardfix_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
export TMPDIR=/global/home/users/kh36969/tmp
cd /global/home/users/kh36969/exclusion_gene_project

python3 "scripts/lib/orf_caller.py" >/dev/null || { echo "caller sentinel FAILED"; exit 1; }

# shard 000 is readable and is rebuilt first as the perfect-score baseline.
# If it is not byte-identical the script stops and rebuilds nothing.
python3 -u scripts/170_rebuild_blocked_shards.py \
    --plsdb /global/scratch/users/kh36969/exclusion_gene/recovery/sequences_r2.fasta \
    --control 0 \
    --shards 2,5,6,9,12,17,18,23,30
