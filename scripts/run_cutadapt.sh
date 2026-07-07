#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: bash scripts/run_cutadapt.sh <PROJECT_DIR> [THREADS]"
    exit 1
fi

DIR="$1"
THREADS="${2:-32}"

mkdir -p "${DIR}/cutadapt"

ADAPTER='TTGTGGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATAGCAA'

ls ${DIR}/rawdata/*/*_1.f*q.gz | sort --version-sort | while read -r id; do
    echo "Processing ${id}"
    RAWDATA=${id%.f*q.gz}
    BASE=${RAWDATA##*/}
    cutadapt -j "${THREADS}" \
        -g "${ADAPTER}" \
        --discard-untrimmed \
        --minimum-length 10 \
        -e 0.1 \
        -q 20 \
        -o "${DIR}/cutadapt/${BASE}.trm.fq.gz" \
        "${id}" > "${DIR}/cutadapt/${BASE}_trimlog.txt"
done
