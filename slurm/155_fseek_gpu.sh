#!/bin/bash
#SBATCH --job-name=fseekgpu
#SBATCH --account=pc_rubinlab
#SBATCH --partition=es1
#SBATCH --qos=es_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/fseekgpu_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/fseekgpu_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
cd /global/home/users/kh36969/exclusion_gene_project
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "no GPU visible"
export FOLDSEEK_BIN=/global/scratch/users/kh36969/bin/gpu/foldseek/bin/foldseek
python3 -u scripts/155_foldseek_remainder.py
