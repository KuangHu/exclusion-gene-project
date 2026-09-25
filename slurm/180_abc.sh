#!/bin/bash
#SBATCH --job-name=abc180
#SBATCH --account=pc_rubinlab
#SBATCH --partition=es1
#SBATCH --qos=es_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/abc180_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/abc180_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
export HF_HOME=/global/home/users/kh36969/.hfcache
export TMPDIR=/global/home/users/kh36969/tmp
cd /global/home/users/kh36969/exclusion_gene_project
python3 -u scripts/180_core_rules_abc.py --dataset "${DS}"
