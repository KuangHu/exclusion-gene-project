#!/bin/bash
# Shared paths for the exclusion-gene database project.
# Source from every stage:  source "$(dirname "$0")/lib/common.sh"

export PROJ_DIR="/global/home/users/kh36969/exclusion_gene_project"
export CONFIG_DIR="${PROJ_DIR}/config"
export DATA_DIR="${PROJ_DIR}/data"

# The seed set is FROZEN IN THE REPO, not on scratch: §7 requires a stable
# accession set per release, and scratch is not a durable medium.
export SEED_DIR="${DATA_DIR}/seed"
export GB_DIR="${SEED_DIR}/genbank"
export EXTRACT_DIR="${SEED_DIR}/extracted"

# Scratch is for bulk expansion corpora only (ICEberg, PLSDB work, step 2+).
export SCRATCH_ROOT="/global/scratch/users/kh36969/exclusion_gene"
export ICEBERG_DIR="${SCRATCH_ROOT}/iceberg"

# Local resources already on this cluster
export PLSDB_DIR="/global/scratch/users/kh36969/plsdb"
export MACSY_MODELS="/global/scratch/users/kh36969/funcannot_dbs/macsy_models"

export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"

# SLURM defaults (cf1 is the partition with current free capacity)
export SLURM_ACCOUNT="pc_rubinlab"
export SLURM_PARTITION="cf1"
export SLURM_QOS="cf_normal"
