# Installation Guide for IFIT2-IFIT3 eCLIP Analysis Pipeline

## Overview

This pipeline requires bioinformatics tools for processing eCLIP sequencing data. The pipeline uses command-line Python scripts and can be run on HPC systems (SLURM) or local machines.

## Quick Start (Recommended)

### Step 1: Install Conda/Mamba

If you don't already have conda installed, download Miniforge:
- **Linux/HPC**: https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
- **Windows**: https://github.com/conda-forge/miniforge/releases/latest
- **macOS**: https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh

```bash
# Linux/HPC installation
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b
~/miniforge3/bin/conda init bash
source ~/.bashrc
```

### Step 2: Create Analysis Environment

**Option A: Install from environment file (recommended)**

```bash
# Clone the repository
cd /path/to/IFIT1-IFIT2-secondary-analysis

# Create environment from file
conda env create -f environment.yml
conda activate eclip
```

**Option B: Manual installation**

```bash
# Create environment
conda create -n eclip python=3.10 -y
conda activate eclip

# Install tools from bioconda
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
pip install pysam pybedtools pyBigWig pandas numpy matplotlib seaborn scipy
pip install git+https://github.com/YeoLab/clipper.git@master
```

### Step 3: Verify Installation

```bash
# Check that tools are installed
which fasterq-dump
which fastqc
which cutadapt
which STAR
which samtools
which bedtools
which umi_tools
which clipper
which perl

# Test Python packages
python -c "import pysam, pybedtools, pandas, numpy, matplotlib, seaborn, scipy; print('All packages installed')"
```

### Step 4: Run the Pipeline

```bash
# On local machine - run complete pipeline
python scripts/run_pipeline.py

# On HPC - submit SLURM jobs
sbatch slurm/setup_reference.sbatch    # One-time setup
sbatch slurm/process_sample.sbatch    # Process all samples
```

---

## Key Tools

### Yeo Lab eCLIP Pipeline

This pipeline uses the official Yeo lab eCLIP analysis tools:

- **CLIPper** ([YeoLab/clipper](https://github.com/YeoLab/clipper)) - Peak caller using Poisson statistics
- **Normalization scripts** - Perl scripts for IP vs input normalization using Fisher's exact test

These are the same tools used in published eCLIP studies and provide statistically rigorous peak identification.

### Required Perl Modules

The Yeo lab normalization scripts require:
- `Statistics::Basic`
- `Statistics::Distributions`
- `Statistics::R`

These are automatically installed when using the provided `environment.yml` file.

---

## System Requirements

### Memory

- **Reference indexing**: 64 GB RAM (one-time)
- **Sample processing**: 32 GB RAM per sample
- **Analysis steps**: 8-16 GB RAM

### Disk Space

- **Reference genome**: ~30 GB
- **Per sample**: ~2-3 GB
- **Total recommended**: 100 GB free space

### Software Versions

Expected versions (as of January 2025):

- **Python**: 3.10+
- **sra-tools**: 3.0+
- **FastQC**: 0.12+
- **cutadapt**: 4.0+
- **STAR**: 2.7+
- **samtools**: 1.17+
- **bedtools**: 2.30+
- **umi_tools**: 1.1+
- **clipper**: 2.0+
- **perl**: 5.30+

---

## HPC-Specific Setup

### CWRU HPC

```bash
# Load modules (if needed)
module load python/3.10
module load conda

# Activate environment
conda activate eclip

# Set up reference genome (requires 64 GB RAM)
sbatch slurm/setup_reference.sbatch

# Process all samples in parallel
sbatch slurm/process_sample.sbatch
```

### SLURM Job Arrays

The pipeline supports parallel processing using SLURM array jobs:

```bash
# Process all 12 samples in parallel
sbatch slurm/process_sample.sbatch

# This creates jobs for array indices 0-11
# Each sample runs independently
```

---

## Troubleshooting

### "STAR: command not found"

Make sure conda environment is activated:
```bash
conda activate eclip
```

### "Permission denied" when downloading from SRA

Configure SRA Toolkit cache:
```bash
vdb-config --interactive
# Set a cache directory with write permissions
```

### Out of memory during STAR indexing

Use HPC or cloud instance with 64+ GB RAM. Alternatively:
```bash
# Reduce threads in script
python scripts/01_setup_reference.py --threads 4
```

### Slow downloads from SRA

Use prefetch before fasterq-dump:
```bash
prefetch SRR31773513
fasterq-dump --outdir ./data/raw_fastq SRR31773513
```

### "Perl module not found"

Ensure Perl modules are installed:
```bash
conda install -c conda-forge \
    perl-statistics-basic \
    perl-statistics-distributions \
    perl-statistics-r
```

---

## Testing Installation

Run this test script to verify all tools:

```bash
#!/bin/bash
# test_installation.sh

echo "Testing tool installation..."
echo "=============================="

echo -n "fasterq-dump: "
fasterq-dump --version 2>&1 | head -1

echo -n "fastqc: "
fastqc --version

echo -n "cutadapt: "
cutadapt --version

echo -n "STAR: "
STAR --version

echo -n "samtools: "
samtools --version | head -1

echo -n "bedtools: "
bedtools --version

echo -n "umi_tools: "
umi_tools --version

echo -n "clipper: "
clipper --version 2>&1 | head -1

echo -n "perl: "
perl --version | head -2 | tail -1

echo -n "Python packages: "
python -c "import pysam, pybedtools, pandas, numpy, matplotlib, seaborn, scipy; print('All installed')"

echo "=============================="
echo "Installation test complete!"
```

Run the test:
```bash
bash test_installation.sh
```

---

## Getting Help

If you encounter issues:

1. **Yeo lab eCLIP**: https://github.com/YeoLab/eclip
2. **CLIPper documentation**: https://github.com/YeoLab/clipper
3. **Bioconda issues**: https://github.com/bioconda/bioconda-recipes/issues
4. **SRA Toolkit**: https://github.com/ncbi/sra-tools/wiki
5. **STAR aligner**: https://github.com/alexdobin/STAR

---

## Next Steps

Once all tools are installed:

1. Set up reference genome: `sbatch slurm/setup_reference.sbatch` (HPC) or `python scripts/01_setup_reference.py` (local)
2. Download samples: `python scripts/02_download_data.py --all`
3. Process samples: `sbatch slurm/process_sample.sbatch` (HPC) or `python scripts/03_process_samples.py --sample SAMPLE_NAME` (local)
4. Call peaks: `python scripts/04_call_peaks.py`
5. Analyze UTRs: `python scripts/05_analyze_utrs.py`
6. Generate figures: `python scripts/06_visualize.py`

See [scripts/README.md](scripts/README.md) for detailed usage instructions.
