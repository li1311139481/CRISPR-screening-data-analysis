#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Demo script for pooled CRISPR screen analysis pipeline
# ============================================================
#
# Usage:
#   Run from the repository root directory:
#
#     cd /path/to/CRISPR-screening-data-analysis
#     bash examples/demo_run.sh /path/to/your/project
#
#   where <project> is a directory containing:
#     - rawdata/<sample_number>/<sample_number>_1.fq.gz
#     - meta.data.csv
#     (edit config/TF_library.csv and config/TF_library_control_id.txt
#      with your actual library sequences and gene annotations)
#
#   Example:
#     DIR=$(pwd)                        # project = the repo itself
#     bash examples/demo_run.sh ${DIR}
#
# Optional knobs: THREADS / PARALLEL can be overridden by setting
#   THREADS=32 PARALLEL=8 bash examples/demo_run.sh ${DIR}
# ============================================================

# Project directory (data) is passed as the first argument.
if [ $# -lt 1 ]; then
    echo "Error: missing project directory argument." >&2
    echo "Usage: bash examples/demo_run.sh /path/to/your/project" >&2
    exit 1
fi
DIR="${1}"

THREADS="${THREADS:-32}"
PARALLEL="${PARALLEL:-8}"

# The scripts below reference ./scripts, ./config and ./rawdata,
# so the repository root must be the current working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${SCRIPT_DIR}"

mkdir -p "${DIR}/cutadapt"

# Step 1: Trim adapters and extract sgRNA spacers from raw FASTQ (R1 reads)
bash scripts/run_cutadapt.sh "${DIR}" "${THREADS}"

# Step 2: Generate pairwise MAGeCK comparison commands and Slurm script
#   All pairwise comparisons among groups in meta.data.csv will be generated.
#   Comparisons that already have VolcanoView.pdf are skipped.
#   Add --submit to automatically submit the Slurm job.
#   --with-postprocess makes run.slurm additionally auto-generate the merged
#   count table, per-sample library QC, and per-comparison MAGeCKFlute reports
#   after count/test finish (the full one-command end-to-end run).
python scripts/screening_run_all.py \
    -d "${DIR}" \
    -l config/TF_library.csv \
    -c config/TF_library_control_id.txt \
    -p "${PARALLEL}" \
    -n TF \
    --with-postprocess

# After Slurm finishes, results are in:
#   ${DIR}/TF_<treatment>.vs.<control>/
