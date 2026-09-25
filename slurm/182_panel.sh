#!/bin/bash
#SBATCH --job-name=panel182
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=3:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/panel182_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/panel182_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
export TMPDIR=/global/home/users/kh36969/tmp
# TMbed raced on the shared HF cache and died in rmtree cleanup
# (OSError 39, Directory not empty). Offline mode stops it trying to
# re-resolve a model that is already on disk.
export HF_HOME=/global/home/users/kh36969/.hfcache
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
cd /global/home/users/kh36969/exclusion_gene_project
python3 -u scripts/182_abc_panel_recall.py
