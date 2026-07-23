# CRISPR screening data analysis

A computational workflow for pooled CRISPR screen sequencing analysis. The workflow extracts sgRNA spacer sequences from FASTQ files, quantifies sgRNA abundance with [MAGeCK](https://sourceforge.net/p/mageck/), performs all pairwise enrichment/depletion analyses, and generates quality-control plots and summary reports.

## Quick start

```bash
# 1. Install conda environment
conda env create -f environment.yml
conda activate crispr-screen-analysis

# If MAGeCKFlute fails to install via conda, install from Bioconductor:
# Rscript -e 'BiocManager::install("MAGeCKFlute")'

# 2. Prepare your project (see Input files below)
# 3. Run the analysis
bash examples/demo_run.sh
```

## Workflow

```
  Raw FASTQ                    meta.data.csv
       │                             │
       ▼                             │
  ┌──────────┐                       │
  │ cutadapt  │  extract sgRNA       │
  │   (R1)    │  spacer sequences    │
  └────┬─────┘                       │
       │  *.trm.fq.gz                │
       ▼                             ▼
  ┌──────────────────────────────────────┐
  │        screening_run_all.py          │
  │  - all pairwise group comparisons    │
  │  - generates commands.sh + run.slurm │
  └────────────────┬─────────────────────┘
                   │
                   ▼
  ┌──────────────────────────────────────┐
  │           mageck.py                  │
  │  ┌──────────┐  ┌──────────┐          │
  │  │ count    │→ │   test   │          │
  │  └──────────┘  └────┬─────┘          │
  │                     ▼                │
  │  ┌──────────────────────────────┐    │
  │  │     mageck_flute.R           │    │
  │  │  QC plots + volcano + report │    │
  │  └──────────────────────────────┘    │
  └──────────────────────────────────────┘
```

## Installation

### Prerequisites

- Linux server with **Miniconda** or **Anaconda** installed. If you don't have it:

  ```bash
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  # restart shell, then configure channel priority:
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda config --set channel_priority flexible
  ```

- **Slurm** for job scheduling (usually pre-installed on HPC clusters).

### Create environment

```bash
cd /path/to/CRISPR-screening-data-analysis
conda env create -f environment.yml
conda activate crispr-screen-analysis
```

### Verify

```bash
cutadapt --version             # 4.8
mageck --version               # 0.5.9.5
parallel --version             # any
python -c "import pandas; print('pandas', pandas.__version__)"
Rscript -e 'library(MAGeCKFlute); message("OK")'
```

## Repository structure

```text
CRISPR-screening-data-analysis/
├── README.md
├── LICENSE
├── .gitignore
├── meta.data.csv
├── environment.yml
├── requirements.txt
├── scripts/
│   ├── run_cutadapt.sh             # Step 1: sgRNA spacer extraction
│   ├── screening_run_all.py        # Step 2: pairwise comparison generator
│   ├── mageck.py                   # Step 2 (worker): MAGeCK count + test
│   ├── mageck_flute.R              # Step 2 (worker): QC & visualization
│   └── plot_library_qc.py          # Library QC: KDE / Skew ratio / AUC plots
├── config/
│   ├── TF_library.csv
│   └── TF_library_control_id.txt
├── rawdata/
│   ├── S1/
│   │   └── S1_1.fq.gz
│   ├── S2/
│   │   └── S2_1.fq.gz
│   └── ...
├── examples/
│   ├── demo_project_structure.txt  # Expected project directory layout
│   ├── demo_run.sh                 # Example run script
│   └── expected_outputs.txt        # Expected output files per comparison
├── docs/
│   ├── protocol_analysis.md
│   ├── input_files.md
│   └── output_files.md
└── results/
    └── README.md
```

> **Note:** The `rawdata/` files in this repository are **downsampled demo data** (80万 reads per sample) for testing the workflow. Full-size data should be prepared following the layout described in [Input files](#input-files).

## Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| cutadapt | >= 4.0 | sgRNA spacer extraction |
| MAGeCK | >= 0.5.9 | sgRNA count and statistical test |
| Python | >= 3.10 | Workflow orchestration |
| R | >= 4.4 | MAGeCKFlute QC and visualization |
| GNU parallel | any | Parallel job execution |
| Slurm | any | HPC job scheduling |

### Python packages

| Package | Purpose |
|---------|---------|
| pandas | Data manipulation |
| numpy | Numerical computation |
| matplotlib | Volcano plot generation |

### R packages

| Package | Purpose |
|---------|---------|
| ggplot2 | Plotting |
| org.Mm.eg.db | Mouse gene annotation |
| clusterProfiler | Gene ID conversion |
| dplyr | Data manipulation |
| data.table | Fast data I/O |
| scales | Axis label formatting |
| MAGeCKFlute | QC plots and volcano plots |

## Input files

### Project directory layout

```text
project/
├── rawdata/
│   ├── S1/
│   │   └── S1_1.fq.gz
│   ├── S2/
│   │   └── S2_1.fq.gz
│   └── ...
└── meta.data.csv
```

### meta.data.csv

Three required columns:

| Column | Description | Example |
|--------|-------------|---------|
| `sample_name` | Biological replicate label for MAGeCK | `Input_1` |
| `sample_number` | Sample ID matching FASTQ subdirectory | `S1` |
| `group_name` | Biological group for pairwise comparisons | `Input` |

Example:

```csv
sample_name,sample_number,group_name
Input_1,S1,Input
Input_2,S2,Input
Day7_Total_1,S3,Day7_Total
Day7_Total_2,S4,Day7_Total
```

### sgRNA library file

Three comma-separated columns (no header):

```csv
sgRNA_ID,spacer_sequence,target_gene
MEM.Aebp1.g1101,TATGCCAGGAGTACCGCGAT,Aebp1
MEM.Aebp1.g1102,GATGTACACCAATGGCTACG,Aebp1
```

### Control sgRNA file

Header `Name`, followed by sgRNA IDs matching the library file:

```text
Name
Neg_control_1
Neg_control_2
```

## Usage

### Step 1: Extract sgRNA spacers

The adapter sequence targets the vector backbone used in this protocol (`TTGTGGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATAGCAA`).

```bash
DIR=/path/to/project
bash scripts/run_cutadapt.sh ${DIR} 32
```

The script processes `rawdata/*/*_1.f*q.gz` (R1 reads only). Trimmed FASTQ files are named `cutadapt/{sample_number}_1.trm.fq.gz`.

### Step 2: Generate and submit MAGeCK comparisons

```bash
python scripts/screening_run_all.py \
    -d /path/to/project \
    -l config/TF_library.csv \
    -c config/TF_library_control_id.txt \
    -p 8 \
    -n TF \
    --submit
```

**Options:**

| Flag | Description |
|------|-------------|
| `-d` | Project directory containing `cutadapt/` |
| `-l` | sgRNA library annotation file (CSV) |
| `-c` | Control sgRNA ID list |
| `-m` | Metadata CSV path (default: `DIR/meta.data.csv`) |
| `-n` | Output prefix for result directories |
| `-p` | Number of parallel MAGeCK jobs |
| `--submit` | Automatically submit the Slurm job |
| `--nodelist` | Slurm nodelist (e.g. node2), omit if not needed |
| `--python` | Python executable path (default: current interpreter) |
| `--mageck-wrapper` | Path to `mageck.py` (default: auto-detect) |

The script generates `n × (n - 1)` directional pairwise comparisons from all groups in `meta.data.csv`. Comparisons that already have `VolcanoView.pdf` are automatically skipped, allowing safe re-runs on incomplete analyses.

Without `--submit`, only `commands.sh` and `run.slurm` are generated. Submit manually with:

```bash
sbatch run.slurm
```

### Step 3: Inspect outputs

Each comparison is written to a directory named:

```text
{prefix}_{treatment}.vs.{control}/
```

Example: `TF_Day7_Total.vs.Input/`

**Output files per comparison:**

| File | Description |
|------|-------------|
| `*_count.count.txt` | Raw sgRNA count matrix |
| `*_count.count_normalized.txt` | Normalized count matrix |
| `*_count.countsummary.txt` | Mapping rate, Gini index, zero-count stats |
| `*_test.gene_summary.txt` | Gene-level RRA statistics |
| `*_test.sgrna_summary.txt` | sgRNA-level statistics |
| `*_gene_summary.csv` | Annotated full gene summary with ENTREZID |
| `*_mix_gene_summary.csv` | Unified positive/negative gene summary |
| `MapRatesView.pdf` | sgRNA mapping rate plot |
| `GiniIndexView.pdf` | Read distribution evenness plot |
| `MissedsgRNAView.pdf` | Missed sgRNA ratio plot |
| `VolcanoView.pdf` | MAGeCKFlute volcano plot (also serves as completion marker) |
| `*_MAGeCKFlute_report.pdf` | Combined QC and volcano report |

## Notes

- If using a different vector backbone, adjust the adapter sequence in `run_cutadapt.sh` accordingly.
- `VolcanoView.pdf` serves as a completion marker — if it exists, that comparison is skipped on re-run.
- Replace `config/TF_library.csv` and `config/TF_library_control_id.txt` with your own library files, and edit `meta.data.csv` to match your sample metadata before running.
- This workflow was developed for mouse CRISPR screens using `org.Mm.eg.db`. For human screens, update `mageck_flute.R` to use `org.Hs.eg.db`.

## Citation

If you use this workflow, please cite:

- MAGeCK: Li W, et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. *Genome Biology*, 2014.
- MAGeCKFlute: Wang B, et al. Integrative analysis of pooled CRISPR genetic screens using MAGeCKFlute. *Nature Protocols*, 2019.

## License

This repository is distributed under the MIT License. See `LICENSE` for details.
