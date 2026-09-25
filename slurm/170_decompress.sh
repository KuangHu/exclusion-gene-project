#!/bin/bash
#SBATCH --job-name=plsdb_bz2
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/global/home/users/kh36969/exclusion_gene_project/logs/plsdb_bz2_%j.out
#SBATCH --error=/global/home/users/kh36969/exclusion_gene_project/logs/plsdb_bz2_%j.err
set -euo pipefail
export PATH="/global/home/users/kh36969/.conda/envs/claude-env/bin:$PATH"
export TMPDIR=/global/home/users/kh36969/tmp

SRC=/global/scratch/users/kh36969/plsdb/sequences.fasta.bz2
DST=/global/scratch/users/kh36969/exclusion_gene/recovery/sequences_r2.fasta

# sequences.fasta is unreadable (Lustre OST 106). The .bz2 is intact. The output
# MUST live in recovery/, which is pinned to the ddn_hdd pool: the default layout puts
# a new file's first PFL component in ddn_nvme2, and OST 106 is 1 of that pool's 4
# members, so ~25% of new files are born unreadable. The first attempt (job 26150157)
# was cancelled for exactly that reason.
echo "decompressing $SRC -> $DST"
bunzip2 -c "$SRC" > "$DST"

N=$(grep -c '^>' "$DST")
echo "records in rebuilt fasta: $N"
EXP=$(wc -l < /global/scratch/users/kh36969/plsdb/sequences.fasta.fai)
echo "records expected (.fai):  $EXP"
if [ "$N" -ne "$EXP" ]; then
    echo "RECORD COUNT MISMATCH -- do not use this file"
    exit 1
fi
echo "OK: rebuilt PLSDB matches the .fai record count"
ls -la "$DST"
