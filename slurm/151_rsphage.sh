#!/bin/bash
#SBATCH --job-name=rsphage
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/rsphage_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/rsphage_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
cd /global/home/users/kh36969/exclusion_gene_project
python3 -u scripts/151_mpf_type_refseq_phage.py
