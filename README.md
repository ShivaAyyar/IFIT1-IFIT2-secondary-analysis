# IFIT2-IFIT3 eCLIP Analysis Scripts

Command-line Python scripts for analyzing IFIT2-IFIT3 eCLIP data from GSE284636 (Glasner et al. 2025).

## Overview

This pipeline analyzes eCLIP (enhanced CLIP-seq) data to identify human transcripts bound by IFIT2-IFIT3 complexes and analyzes their 5' UTR length preferences.

**Key Features:**
- Processes single-end 75nt eCLIP reads
- Uses Yeo lab's official eCLIP analysis pipeline (CLIPper + normalization scripts)
- Distinguishes IFIT2-only, IFIT3-only, and complex binding
- Analyzes 5' UTR length distributions
- Generates publication-quality figures

## Pipeline Workflow

### Complete Analysis Steps

```
┌─────────────────────────────────────────────────────────────┐
│ STEP 1: One-Time Setup (Run Once)                          │
├─────────────────────────────────────────────────────────────┤
│ sbatch slurm/setup_reference.sbatch                        │
│   → Downloads hg19 genome                                   │
│   → Downloads GENCODE v19 annotations                       │
│   → Builds STAR index (~64 GB RAM, ~2 hours)               │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ STEP 2: Download FASTQ Data (Required)                     │
├─────────────────────────────────────────────────────────────┤
│ sbatch slurm/download_data.sbatch   (OR run on login node) │
│   → Downloads all 12 samples from SRA                       │
│   → Saves to data/raw_fastq/                               │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ STEP 3: Sample Processing (Choose A or B)                  │
├─────────────────────────────────────────────────────────────┤
│ Option A: Parallel (Recommended)                           │
│   sbatch slurm/process_sample.sbatch                       │
│     → Processes 12 samples simultaneously                   │
│     → Each: QC → Trim → Align → Deduplicate               │
│                                                             │
│ Option B: Sequential (Simpler)                             │
│   sbatch slurm/run_all_samples.sbatch                      │
│     → Runs entire pipeline in one job                       │
│     → Automatically proceeds to Steps 4-6                   │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ STEP 4: Peak Calling (After sample processing completes)   │
├─────────────────────────────────────────────────────────────┤
│ python scripts/04_call_peaks.py                            │
│   → Calls peaks with CLIPper (Poisson statistics)          │
│   → Normalizes IP vs input (Fisher's exact test)           │
│   → Filters by fold-change and p-value                     │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ STEP 5: UTR Analysis                                       │
├─────────────────────────────────────────────────────────────┤
│ python scripts/05_analyze_utrs.py                          │
│   → Parses GENCODE GTF for 5' UTRs                         │
│   → Identifies IFIT2/IFIT3/complex bound transcripts       │
│   → Analyzes UTR length distributions                      │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ STEP 6: Generate Figures                                   │
├─────────────────────────────────────────────────────────────┤
│ python scripts/06_visualize.py                             │
│   → Creates comparison plots                                │
│   → Generates UTR length histograms                         │
│   → Publication-quality PDFs                                │
└─────────────────────────────────────────────────────────────┘
```

## Quick Start

### On CWRU HPC

#### **Step 1: One-Time Setup (First Time Only)**

```bash
# Clone repository
cd /home/ssa163
git clone <repository_url> IFIT1-IFIT2-secondary-analysis
cd IFIT1-IFIT2-secondary-analysis

# Create conda environment
conda env create -f environment.yml
conda activate eclip

# Set up reference genome (requires 64 GB RAM)
sbatch slurm/setup_reference.sbatch

# Monitor job completion
squeue -u ssa163
```

#### **Step 2: Download FASTQ Data (Required)**

**Option A: Submit as SLURM job (Recommended for HPC)**

```bash
# Submit download job
sbatch slurm/download_data.sbatch

# Monitor progress
squeue -u ssa163
```

**Option B: Run on login node**

```bash
# Download all 12 samples from SRA
python scripts/02_download_data.py --all

# OR download specific samples
python scripts/02_download_data.py --sample IFIT2_IFIT3_FLAG_IP
```

**Note:** This step is **required** before processing. The FASTQ files must be downloaded to `data/raw_fastq/` before running the SLURM jobs.

#### **Step 3: Main Analysis (Choose Option A or B)**

**Option A: Parallel Processing (Recommended - Fastest)**

```bash
# 1. Process all 12 samples in parallel
sbatch slurm/process_sample.sbatch

# 2. Wait for ALL jobs to complete (check with: squeue -u ssa163)

# 3. After all samples finish, run analysis steps (Steps 4-6)
python scripts/04_call_peaks.py
python scripts/05_analyze_utrs.py
python scripts/06_visualize.py
```

**Option B: Sequential Processing (Simpler - Slower)**

```bash
# Run entire pipeline as one job
sbatch slurm/run_all_samples.sbatch

# All analysis runs automatically in this job
```

### On Local Machine

```bash
# Run complete pipeline
python scripts/run_pipeline.py
```

## Script Descriptions

### Core Modules

- **config.py** - Configuration, sample information, and parameters
- **utils.py** - Utility functions (download, QC, alignment, etc.)
- **analysis.py** - Analysis functions (peaks, UTR analysis)
- **visualization.py** - Plotting functions

### Executable Scripts

1. **01_setup_reference.py** - Download hg19 and build STAR index
   ```bash
   python scripts/01_setup_reference.py --threads 8
   ```

2. **02_download_data.py** - Download FASTQ files from SRA
   ```bash
   # Download all samples
   python scripts/02_download_data.py --all

   # Download specific sample
   python scripts/02_download_data.py --sample IFIT2_IFIT3_FLAG_IP
   ```

3. **03_process_samples.py** - QC, trim, align, and deduplicate
   ```bash
   # Process specific sample
   python scripts/03_process_samples.py --sample IFIT2_IFIT3_FLAG_IP --threads 8

   # Process by index (for SLURM arrays)
   python scripts/03_process_samples.py --sample-index 0
   ```

4. **04_call_peaks.py** - Call peaks for all IP samples
   ```bash
   python scripts/04_call_peaks.py
   ```

5. **05_analyze_utrs.py** - Analyze 5' UTR binding patterns
   ```bash
   python scripts/05_analyze_utrs.py
   ```

6. **06_visualize.py** - Generate publication figures
   ```bash
   python scripts/06_visualize.py
   ```

7. **run_pipeline.py** - Run complete pipeline
   ```bash
   # Full pipeline
   python scripts/run_pipeline.py

   # Skip reference setup
   python scripts/run_pipeline.py --skip-reference

   # Skip download
   python scripts/run_pipeline.py --skip-download
   ```

## SLURM Batch Scripts

### setup_reference.sbatch
Submits reference genome setup as batch job (requires 64 GB RAM).

```bash
sbatch slurm/setup_reference.sbatch
```

### process_sample.sbatch
Processes all 12 samples in parallel using array jobs.

```bash
sbatch slurm/process_sample.sbatch
# This creates 12 parallel jobs (array indices 0-11)
```

### run_all_samples.sbatch
Runs the complete pipeline sequentially.

```bash
sbatch slurm/run_all_samples.sbatch
```

## Sample Information

The pipeline processes 12 samples from GSE284636:

### IFIT2+IFIT3 Co-expression
- `IFIT2_IFIT3_FLAG_IP` (SRR31773513) - IFIT2 binding
- `IFIT2_IFIT3_HA_IP` (SRR31773515) - IFIT3 binding
- `IFIT2_IFIT3_FLAG_SMInput` (SRR31773512)
- `IFIT2_IFIT3_HA_SMInput` (SRR31773514)

### IFIT2 Alone
- `IFIT2_alone_FLAG_IP` (SRR31773519)
- `IFIT2_alone_SMInput` (SRR31773518)

### IFIT3 Alone
- `IFIT3_alone_HA_IP` (SRR31773517)
- `IFIT3_alone_SMInput` (SRR31773516)

### Wild-Type Controls
- `WT_FLAG_IP` (SRR31773539)
- `WT_HA_IP` (SRR31773541)
- `WT_FLAG_SMInput` (SRR31773538)
- `WT_HA_SMInput` (SRR31773540)

## Output Files

### Directory Structure

```
IFIT1-IFIT2-secondary-analysis/
├── data/
│   ├── reference/
│   │   ├── hg19.fa.gz
│   │   ├── gencode.v19.annotation.gtf.gz
│   │   ├── STAR_index_hg19/
│   │   └── 5utr_regions.bed
│   ├── raw_fastq/          # Downloaded FASTQ files
│   ├── trimmed/            # UMI-extracted, adapter-trimmed reads
│   ├── aligned/            # STAR-aligned, deduplicated BAMs
│   ├── peaks/              # Called peaks (BED files)
│   └── results/
│       ├── ifit2_alone_bound_transcripts.txt
│       ├── ifit3_alone_bound_transcripts.txt
│       ├── complex_bound_transcripts.txt
│       ├── utr_length_distribution.csv
│       └── figures/
│           ├── ifit_comparison.pdf
│           └── utr_length_histogram.pdf
└── logs/                   # Log files for all steps
```

### Key Output Files

- **bound transcripts lists** - Transcripts bound by IFIT2, IFIT3, or complex
- **utr_length_distribution.csv** - Complete UTR analysis data
- **ifit_comparison.pdf** - 6-panel comparison figure
- **utr_length_histogram.pdf** - UTR length distribution plots

## System Requirements

### Memory
- **Reference indexing**: 64 GB RAM (one-time)
- **Sample processing**: 32 GB RAM per sample
- **Analysis steps**: 8-16 GB RAM

### Disk Space
- **Reference genome**: ~30 GB
- **Per sample**: ~2-3 GB
- **Total recommended**: 100 GB free space

### Software Dependencies

**Recommended:** Use the provided environment file:

```bash
# From the project root directory
conda env create -f environment.yml
conda activate eclip
```

**Or install manually:**

```bash
# Install conda packages
conda install -c bioconda -c conda-forge \
    sra-tools \
    fastqc \
    cutadapt \
    star \
    samtools \
    bedtools \
    umi_tools \
    perl \
    perl-statistics-basic \
    perl-statistics-distributions \
    perl-statistics-r \
    cython \
    -y

# Install Python packages
pip install pysam pybedtools pandas matplotlib seaborn scipy

# Install CLIPper from GitHub (specific commit used by eCLIP pipeline)
# Note: --no-build-isolation allows CLIPper to use conda-installed numpy/cython
pip install --no-build-isolation git+https://github.com/YeoLab/clipper.git@5d865bb17b2bc6787b4c382bc857119ae917ad59
```

**Key Tools:**

This pipeline uses the official Yeo lab eCLIP analysis tools from https://github.com/yeolab/eclip:

- **CLIPper** ([YeoLab/clipper](https://github.com/YeoLab/clipper)) - Peak caller using Poisson-based statistics
- **Normalization scripts** - Perl scripts (`overlap_peakfi_with_bam.pl`) for IP vs input normalization using Fisher's exact test
- **pybedtools** - Genomic overlap analysis
- **scipy** - Statistical testing (Mann-Whitney U test for UTR analysis)

## Troubleshooting

### "STAR: command not found"
Ensure conda environment is activated:
```bash
conda activate eclip
```

### "Permission denied" when downloading from SRA
Configure SRA Toolkit cache:
```bash
vdb-config --interactive
```

### Out of memory during STAR indexing
Use HPC or cloud instance with 64+ GB RAM. This step cannot run on machines with less memory.

### Slow downloads
Use prefetch before fasterq-dump:
```bash
prefetch SRR31773513
fasterq-dump --outdir ./data/raw_fastq SRR31773513
```

## Citation

If you use this pipeline, please cite:

- **Data source**: Glasner et al. (2025) - GSE284636
- **STAR aligner**: Dobin et al. (2013)
- **GENCODE**: Frankish et al. (2019)

## Peak Calling Method

This pipeline uses the Yeo lab's standard eCLIP analysis approach:

1. **Peak Calling**: CLIPper identifies enriched regions in IP samples using a Poisson-based statistical model
2. **Normalization**: Yeo lab Perl scripts (`overlap_peakfi_with_bam.pl`) normalize IP peaks against size-matched input controls using Fisher's exact test or chi-square test
3. **Filtering**: Peaks are filtered by fold-change (≥2.0) and p-value (≤0.001) thresholds
4. **Output**: Normalized peaks with -log10(p-value) and log2(fold-change) values

This approach matches the published eCLIP methods and provides statistically rigorous peak identification.

## Additional Documentation

For detailed installation and usage instructions, see:
- [Installation Guide](INSTALLATION_GUIDE.md) - Complete setup instructions for HPC and local machines
- [environment.yml](environment.yml) - Conda environment specification
