#!/bin/bash
#SBATCH --job-name=round0
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/round0_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/round0_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
cd /global/home/users/kh36969/exclusion_gene_project
python3 -u scripts/140_round0_cluster_v11.py --tag "${TAG:-VirB5+1}" --cpus 16
