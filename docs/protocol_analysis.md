# Computational analysis of pooled CRISPR screen sequencing data

This workflow extracts sgRNA spacer sequences from raw FASTQ files, quantifies sgRNA abundance, identifies genes enriched or depleted between experimental groups, and generates quality-control and visualization outputs. The workflow uses `meta.data.csv` to define sample groups and automatically generates all directional pairwise MAGeCK comparisons.

## Main steps

1. Prepare project directory and metadata.
2. Extract sgRNA spacers using cutadapt.
3. Generate all pairwise MAGeCK comparisons with `screening_run_all.py`.
4. Run `mageck count` and `mageck test` with `mageck.py`.
5. Generate QC plots and summary reports with `mageck_flute.R`.

## Notes

- The default downstream wrapper expects trimmed FASTQ files named as `{sample_number}_1.trm.fq.gz`.
- The script generates directional comparisons. For `n` groups, it generates `n × (n - 1)` comparisons.
- The visualization step generates mapping rate, Gini index, missed sgRNA ratio plots, volcano plots, and an integrated MAGeCKFlute report.
