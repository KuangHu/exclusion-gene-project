#!/bin/bash
#SBATCH --job-name=mash_plsdb
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/81_mash_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/81_mash_%j.err
#
# Pipeline v2 §0.1 -- all-vs-all Mash distances over PLSDB + the RP4 control.
#
# Purpose: pick the ANI clustering threshold from the DATA, not from convention.
# 95% is the usual default but plasmid ANI structure differs from chromosomal
# (shared backbone + variable accessory), so the inter-cluster gap may not fall
# where it does for genomes. Getting this wrong re-runs every downstream statistic.
#
# Two things this must answer:
#   1. does the pairwise ANI distribution have a valley? -> put the threshold there
#   2. at a candidate threshold, do the controls land in SEPARATE clusters?
#      RP4 (IncP-alpha) and R751 (IncP-beta) merging would be disqualifying:
#      their TrbK proteins are only 37.8% identical, so they are plainly
#      independent evolutionary events and must not count as one observation.
#
# -d 0.10 keeps every threshold from 90% ANI upward available without a re-run.
set -euo pipefail
export PATH="/global/home/groups/pc_rubinlab/databases/kuang/system/.conda/envs/mob_suite/bin:$PATH"

WORK=/global/scratch/users/kh36969/exclusion_gene/mash
cd "$WORK"

[[ -s plsdb_plus.msh ]] || { echo "FATAL: plsdb_plus.msh missing"; exit 1; }
N=$(mash info -t plsdb_plus.msh | tail -n +2 | wc -l)
echo "$(date)  sketches: ${N}  threads: ${SLURM_CPUS_PER_TASK}"
[[ "$N" -eq 72557 ]] || { echo "FATAL: expected 72557 sketches, got $N"; exit 1; }

mash dist -p "${SLURM_CPUS_PER_TASK}" -d 0.10 \
    plsdb_plus.msh plsdb_plus.msh \
    | awk -F'\t' '$1!=$2' \
    | gzip -1 > pairs_d010.tsv.gz

echo "$(date)  done"
echo "pairs (excluding self): $(zcat pairs_d010.tsv.gz | wc -l)"
ls -la pairs_d010.tsv.gz
