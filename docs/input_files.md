# Input files

## Raw FASTQ files

Raw FASTQ files should be organized under the project `rawdata` directory:

```text
rawdata/S1/S1_1.fq.gz
rawdata/S1/S1_2.fq.gz
```

## meta.data.csv

Required columns:

| Column | Description |
|---|---|
| `sample_name` | Biological replicate label used by MAGeCK |
| `sample_number` | Sample identifier used to locate FASTQ files |
| `group_name` | Biological group used for pairwise comparisons |

## sgRNA library file

Three comma-separated columns:

1. sgRNA ID
2. sgRNA spacer sequence
3. Target gene symbol

## Control sgRNA file

A text file with header `Name`. The sgRNA identifiers must match entries in the first column of the sgRNA library file.
