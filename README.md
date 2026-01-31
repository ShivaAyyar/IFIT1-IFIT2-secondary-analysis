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

## Quick Start

### On CWRU HPC

```bash
# 1. Activate conda environment
conda activate eclip

# 2. Set up reference genome (one-time, submit as SLURM job)
sbatch slurm/setup_reference.sbatch

# 3. Process all samples in parallel
sbatch slurm/process_sample.sbatch

# 4. Run analysis and generate figures
python scripts/04_call_peaks.py
python scripts/05_analyze_utrs.py
python scripts/06_visualize.py
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

# Install Python packages and CLIPper from GitHub
pip install pysam pybedtools pandas numpy matplotlib seaborn scipy
pip install git+https://github.com/YeoLab/clipper.git@master
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
