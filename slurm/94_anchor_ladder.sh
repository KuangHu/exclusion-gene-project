#!/bin/bash
#SBATCH --job-name=ladder
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/94_ladder_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/94_ladder_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
P=/global/home/users/kh36969/exclusion_gene_project
python3 "$P/scripts/74_validate.py" || { echo "A9/validate FAILED"; exit 1; }
python3 "$P/scripts/94_anchor_ladder.py" --cpus 32
