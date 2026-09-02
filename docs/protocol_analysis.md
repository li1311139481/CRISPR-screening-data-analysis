# Computational analysis of pooled CRISPR screen sequencing data

This workflow extracts sgRNA spacer sequences from raw FASTQ files, quantifies sgRNA abundance, identifies genes enriched or depleted between experimental groups, and generates quality-control and visualization outputs. The workflow uses `meta.data.csv` to define sample groups and automatically generates all directional pairwise MAGeCK comparisons.

## Main steps

1. Prepare project directory and metadata (`meta.data.csv`, sgRNA library and control files).
2. Extract sgRNA spacers with cutadapt (`run_cutadapt.sh`).
3. Generate and run all pairwise MAGeCK comparisons with `screening_run_all.py` (calls `mageck.py`, which runs `mageck count` + `mageck test`).
4. (Optional, one-command end-to-end) Pass `--with-postprocess` so that after count/test the `postprocess.py` step merges count tables, generates per-sample library QC, and produces MAGeCKFlute reports for every comparison.
5. Inspect outputs: per-comparison count/test results, `all_samples.count.txt`, per-sample `library_qc_*.pdf`, and per-comparison MAGeCKFlute reports.

## Notes

- The default downstream wrapper expects trimmed FASTQ files named as `{sample_number}_1.trm.fq.gz`.
- The script generates directional comparisons. For `n` groups, it generates `n × (n - 1)` comparisons.
- The visualization step generates mapping rate, Gini index, missed sgRNA ratio plots, volcano plots, and an integrated MAGeCKFlute report.
- Re-generate all figures on demand with `plot_library_qc.py --all` (library QC) and `postprocess.py --flute-only` (MAGeCKFlute reports).
