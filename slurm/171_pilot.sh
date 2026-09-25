#!/bin/bash
#SBATCH --job-name=mpfpilot
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_debug
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:20:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/mpfpilot_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/mpfpilot_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
export TMPDIR=/global/home/users/kh36969/tmp
cd /global/home/users/kh36969/exclusion_gene_project
python3 -u scripts/171_pilot_gbc_fata.py
