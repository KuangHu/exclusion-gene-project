#!/bin/bash
#SBATCH --job-name=v11cat
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/v11cat_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/v11cat_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
python3 -u /global/home/users/kh36969/exclusion_gene_project/scripts/130_catalogue_v1_1.py --cpus 16
