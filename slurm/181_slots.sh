#!/bin/bash
#SBATCH --job-name=slots181
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=10:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/slots181_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/slots181_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
export TMPDIR=/global/home/users/kh36969/tmp
cd /global/home/users/kh36969/exclusion_gene_project
python3 -u scripts/181_slot_pipeline.py --dataset "${1:-plsdb}" --classes "${2:-T,F,I,FA}"
