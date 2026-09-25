#!/bin/bash
#SBATCH --job-name=mpfphage
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=20:00:00
#SBATCH --array=0-31%8
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/mpfphage_%A_%a.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/mpfphage_%A_%a.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
cd /global/home/users/kh36969/exclusion_gene_project
python3 -u scripts/150_mpf_type_phage.py
