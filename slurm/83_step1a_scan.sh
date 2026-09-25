#!/bin/bash
#SBATCH --job-name=step1a_scan
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --array=0-31
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/83_step1a_%A_%a.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/83_step1a_%A_%a.err
#
# Pipeline v2 step 0 + 1a, fused: CDS calling and the anchor scan in one pass, so
# 72k protein FASTAs never hit disk.
#
#   step 0  Pyrodigal -p anon  (scripts/lib/orf_caller.py, the only permitted caller)
#   step 1a pyhmmer, all anchors from data/anchors/anchor_set.tsv
#           ENTRY CRITERION: VirB4 = PF03135 at GA -- the sole filter (v2 3.1/3.2)
#
# Anchor thresholds are per-anchor, from the anchor_set file:
#   PF03135 GA        VirB4, entry
#   PF04610 GA        VirB6, slot downstream boundary
#   PF07996 E<=1e-5 AND aa>150   VirB5, empty-slot test
#                     (GA misses BOTH IncP controls -- v2 1.1a Finding 2; the
#                      length conjunct is what makes the loose E-value safe,
#                      because a too-loose VirB5 call turns an OCCUPIED slot into
#                      an "empty" one, which is the unfilterable direction)
#
# Exclusion-gene families are NOT scanned here. Post-hoc classification only (v2 1.3).
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"

PROJ=/global/home/users/kh36969/exclusion_gene_project
OUT=/global/scratch/users/kh36969/exclusion_gene/step1a
mkdir -p "$OUT"

# sentinels first -- v2 §10: a baseline expected to score perfectly must pass
# before any result is read
python3 "$PROJ/scripts/lib/orf_caller.py" > /dev/null || { echo "FATAL: caller sentinel failed"; exit 1; }
python3 "$PROJ/scripts/lib/adapters.py"   > /dev/null || { echo "FATAL: adapter golden test failed"; exit 1; }
echo "$(date)  sentinels OK  task ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT}"

python3 "$PROJ/scripts/83_step1a_scan.py" \
    --shard "${SLURM_ARRAY_TASK_ID}" \
    --nshards "${SLURM_ARRAY_TASK_COUNT}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --out "$OUT"

echo "$(date)  task ${SLURM_ARRAY_TASK_ID} done"
