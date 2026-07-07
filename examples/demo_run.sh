#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Demo script for pooled CRISPR screen analysis pipeline
# ============================================================
#
# Before running:
#   1. Edit DIR to point to your project directory
#   2. Prepare raw FASTQ files under ${DIR}/rawdata/<sample_number>/
#   3. Edit config/meta.data.csv to match your sample metadata
#   4. Edit config/TF_library.csv and config/TF_library_control_id.txt
#      with your actual library sequences and gene annotations
#   5. Install conda environment: conda env create -f environment.yml

DIR=/path/to/project
THREADS=32
PARALLEL=8

mkdir -p ${DIR}/cutadapt

# Step 1: Trim adapters and extract sgRNA spacers from raw FASTQ (R1 reads)
bash scripts/run_cutadapt.sh ${DIR} ${THREADS}

# Step 2: Generate pairwise MAGeCK comparison commands and Slurm script
#   All pairwise comparisons among groups in meta.data.csv will be generated.
#   Comparisons that already have VolcanoView_ZF.pdf are skipped.
#   Add --submit to automatically submit the Slurm job.
python scripts/screening_run_all.py \
    -d ${DIR} \
    -l config/TF_library.csv \
    -c config/TF_library_control_id.txt \
    -p ${PARALLEL} \
    -n TF

# After Slurm finishes, results are in:
#   ${DIR}/TF_<treatment>.vs.<control>/
