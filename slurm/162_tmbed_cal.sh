#!/bin/bash
#SBATCH --job-name=tmbedcal
#SBATCH --account=pc_rubinlab
#SBATCH --partition=es1
#SBATCH --qos=es_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/tmbedcal_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/tmbedcal_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
export HF_HOME=/global/home/users/kh36969/.hfcache
export TMPDIR=/global/home/users/kh36969/tmp
cd /global/home/users/kh36969/exclusion_gene_project
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# CALIBRATION FIRST. Nothing downstream is readable until the i/o label
# convention is pinned on proteins with known topology:
#   CrcB  (P37002, 127 aa) -- residue 127 is PERIPLASMIC -> anchors the 'o' label
#   RstR  (O85264, 111 aa) -- cytoplasmic repressor, ZERO TM -> anchors 'i' and
#                             the zero-TM false-positive rate
cat data/anchors/calibration.faa data/anchors/panel.faa > /global/home/users/kh36969/tmp/tmbed_in.faa
grep -c '^>' /global/home/users/kh36969/tmp/tmbed_in.faa

tmbed predict \
  -f /global/home/users/kh36969/tmp/tmbed_in.faa \
  -p data/anchors/tmbed_panel.pred \
  --out-format 0 \
  --model-dir /global/home/users/kh36969/.hfcache/prott5 \
  --use-gpu

echo "--- raw predictions ---"
cat data/anchors/tmbed_panel.pred
