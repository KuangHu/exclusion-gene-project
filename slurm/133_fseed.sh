#!/bin/bash
#SBATCH --job-name=fseed
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/fseed_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/fseed_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
python3 -u /global/home/users/kh36969/exclusion_gene_project/scripts/133_mpffa_seed_readout.py
