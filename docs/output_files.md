# Output files

Each comparison is written to a folder named:

```text
{prefix}_{treatment}.vs.{control}/
```

Major output files:

| File | Description |
|---|---|
| `*_count.count.txt` | Raw sgRNA count matrix |
| `*_count.count_normalized.txt` | Normalized sgRNA count matrix |
| `*_count.countsummary.txt` | Mapping rate, total sgRNAs, zero-count sgRNAs, and Gini index |
| `*_test.gene_summary.txt` | MAGeCK gene-level enrichment/depletion statistics |
| `*_test.sgrna_summary.txt` | MAGeCK sgRNA-level statistics |
| `*_mix_gene_summary.csv` | Unified positive/negative gene-level summary table |
| `MapRatesView.pdf` | sgRNA mapping-rate plot |
| `GiniIndexView.pdf` | sgRNA read-distribution evenness plot |
| `MissedsgRNAView.pdf` | Missed sgRNA ratio plot |
| `VolcanoView.pdf` | MAGeCKFlute volcano plot (also serves as completion marker) |
| `*_MAGeCKFlute_report.pdf` | Integrated report |

The complete output produced by the demo dataset is committed under [`expected_results/`](../expected_results/), with one subfolder per comparison plus the per-sample library QC plots. Use it as a reference when verifying a new run.
