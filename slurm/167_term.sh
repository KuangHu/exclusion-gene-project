#!/bin/bash
#SBATCH --job-name=term
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/term_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/term_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
export HF_HOME=/global/home/users/kh36969/.hfcache
export TMPDIR=/global/home/users/kh36969/tmp
cd /global/home/users/kh36969/exclusion_gene_project
python3 -u scripts/167_terminator_pilot.py
