# CRISPR screening data analysis

A computational workflow for pooled CRISPR screen sequencing analysis. The workflow extracts sgRNA spacer sequences from FASTQ files, quantifies sgRNA abundance with [MAGeCK](https://sourceforge.net/p/mageck/), performs all pairwise enrichment/depletion analyses, and generates quality-control plots and summary reports.

## Quick start

```bash
# 1. Clone and enter the repository
git clone https://github.com/li1311139481/CRISPR-screening-data-analysis.git
cd CRISPR-screening-data-analysis
DIR=$(pwd)

# 2. Install conda environment
conda env create -f environment.yml
conda activate crispr-screen-analysis

# If MAGeCKFlute fails to install via conda, install from Bioconductor:
# Rscript -e 'BiocManager::install("MAGeCKFlute")'

# 3. Prepare your own project data (see Input files below):
#    - rawdata/<sample_number>/<sample_number>_1.fq.gz
#    - meta.data.csv

# 4. Run the full analysis. ${DIR} is your project directory containing
#    rawdata/ + meta.data.csv.
bash examples/demo_run.sh ${DIR}
sbatch run.slurm

# After the Slurm job finishes, these outputs are produced automatically:
#   - per-comparison MAGeCK count/test results
#   - per-comparison MAGeCKFlute reports (VolcanoView.pdf, QC plots, ...)
#   - per-sample library QC plots at the project root (library_qc_<sample>.pdf)
#   - a merged all-samples count table at all_samples.count.txt
#
# 5. (Optional) Re-generate figures on demand. See "Step 4: Batch plotting" below.
python scripts/plot_library_qc.py -i ${DIR}/all_samples.count.txt --all -o ${DIR}/
python scripts/postprocess.py -d ${DIR} -n TF --flute-only
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
                   │  (parallel, per comparison)
                   ▼
  ┌──────────────────────────────────────┐
  │           mageck.py                  │
  │  ┌──────────┐   ┌──────────┐         │
  │  │  count   │ → │   test   │         │
  │  └──────────┘   └────┬─────┘         │
  │                      │ writes run_flute.sh per comparison │
  └──────────────────────┼────────────────┘
                         ▼
  ┌──────────────────────────────────────┐
  │           postprocess.py             │
  │  - merge all count tables            │
  │    -> all_samples.count.txt          │
  │  - per-sample library QC             │
  │    (plot_library_qc.py, overwrites)  │
  │  - per-comparison MAGeCKFlute report │
  │    (mageck_flute.R via run_flute.sh) │
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
cutadapt --version             # any (no specific version pinned)
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
│   ├── postprocess.py              # Step 2 (worker): merge counts + library QC + MAGeCKFlute
│   ├── mageck_flute.R              # MAGeCKFlute QC & visualization (R)
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
| cutadapt | any | sgRNA spacer extraction |
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
| scipy | Library QC density estimation |

All Python packages above are already declared in `environment.yml`, so `conda env create -f environment.yml` sets everything up. `requirements.txt` is provided as a lightweight alternative (e.g. `pip install -r requirements.txt`) when you only need the Python components and will manage MAGeCK/R separately.

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
Day7.5_Total_1,S3,Day7.5_Total
Day7.5_Total_2,S4,Day7.5_Total
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
| `--partition` | Slurm partition name (omit to use the cluster default) |
| `--python` | Python executable to embed in generated commands (default `python`, resolved from PATH after `conda activate`) |
| `--mageck-wrapper` | Path to `mageck.py` (default: auto-detect) |
| `--with-postprocess` | Also merge counts + generate library QC + MAGeCKFlute reports after count/test (used by `demo_run.sh`) |

The script generates `n × (n - 1)` directional pairwise comparisons from all groups in `meta.data.csv`. Comparisons that already have `VolcanoView.pdf` are automatically skipped, allowing safe re-runs on incomplete analyses.

Without `--submit`, only `commands.sh` and `run.slurm` are generated. Submit manually with:

```bash
sbatch run.slurm
```

By default, `run.slurm` only runs the parallel MAGeCK **count/test** jobs. To also
auto-generate the merged count table, per-sample library QC, and per-comparison
MAGeCKFlute reports after count/test finish, pass `--with-postprocess` (this is what
`examples/demo_run.sh` does):

```bash
python scripts/screening_run_all.py \
    -d ${DIR} -l config/TF_library.csv -c config/TF_library_control_id.txt \
    -p 8 -n TF --submit --with-postprocess
```

With `--with-postprocess`, `run.slurm` runs `scripts/postprocess.py` after `parallel`, which:

1. merges every comparison's `*_count.count.txt` into a single project-level table `all_samples.count.txt`;
2. generates a per-sample library QC plot `library_qc_<sample>.pdf` in the project root;
3. runs the per-comparison MAGeCKFlute report for every comparison directory (`bash {dir}/run_flute.sh`).

So you only need to submit the job once; all outputs are produced end-to-end.

### Step 3: Inspect outputs

Each comparison is written to a directory named:

```text
{prefix}_{treatment}.vs.{control}/
```

Example: `TF_Day7.5_Total.vs.Input/`

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

**Project-level outputs (written to the project root by the post-processing step):**

| File | Description |
|------|-------------|
| `all_samples.count.txt` | Merged sgRNA count table across all samples |
| `library_qc_<sample>.pdf` | Per-sample library QC (KDE / skew ratio / AUC) |

### Step 4: Batch plotting (re-generate all figures in one command)

The automated `run.slurm` already produces these figures. If you want to re-generate
them **on demand** after the analysis, two commands cover everything:

**a) All-samples library QC** (one PDF per sample, project root):

```bash
python scripts/plot_library_qc.py \
  -i ${DIR}/all_samples.count.txt \
  --all \
  -o ${DIR}/
```

**b) MAGeCKFlute report for every comparison** (VolcanoView + QC plots +
integrated report per comparison directory):

```bash
python scripts/postprocess.py \
  -d ${DIR} \
  -n TF \
  --flute-only
```

Both commands **overwrite** any existing PDFs, so you can re-run them at any time
without having to delete old output files first.

For **one** comparison as a quick check, run that comparison's own script
(generated by `mageck.py`):

```bash
Rscript scripts/mageck_flute.R ${DIR}/TF_Day7.5_Total.vs.Input/ \
  Input_1,Input_2,Day7.5_Total_1,Day7.5_Total_2,Day7.5_Total_3 Day7.5_Total.vs.Input
```

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
