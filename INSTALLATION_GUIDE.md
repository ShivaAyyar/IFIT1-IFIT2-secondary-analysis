# Installation Guide for IFIT2-IFIT3 eCLIP Analysis Pipeline

## Overview

This pipeline requires several bioinformatics tools for processing eCLIP sequencing data. Since you're on Windows, the recommended approach is to use **Windows Subsystem for Linux (WSL)** or **conda/mamba** environments.

## Option 1: Using Conda/Mamba (Recommended for Windows)

### Step 1: Install Miniforge (includes conda and mamba)

If you don't already have conda installed:

1. Download Miniforge for Windows: https://github.com/conda-forge/miniforge/releases/latest
2. Run the installer and follow the prompts
3. Open "Miniforge Prompt" from Start Menu

### Step 2: Create Analysis Environment

```bash
# Create a new conda environment for eCLIP analysis
conda create -n eclip python=3.10 -y
conda activate eclip

# Install bioinformatics tools from bioconda
conda install -c bioconda -c conda-forge \
    sra-tools \
    fastqc \
    cutadapt \
    star \
    samtools \
    bedtools \
    umi_tools \
    -y

# Install Python packages
pip install pysam pybedtools pyBigWig pandas numpy matplotlib seaborn scipy jupyter
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
```

### Step 4: Launch Jupyter Notebook

```bash
# Navigate to your project directory
cd "C:\Users\s3ayy\OneDrive\Documents\GitHub\IFIT1-IFIT2-secondary-analysis"

# Launch Jupyter
jupyter notebook
```

---

## Option 2: Using WSL (Windows Subsystem for Linux)

### Step 1: Install WSL

Open PowerShell as Administrator and run:

```powershell
wsl --install -d Ubuntu
```

Restart your computer when prompted.

### Step 2: Set Up Ubuntu Environment

```bash
# Update package list
sudo apt update
sudo apt upgrade -y

# Install basic dependencies
sudo apt install -y wget curl git build-essential

# Install Miniconda in WSL
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
~/miniconda3/bin/conda init bash
source ~/.bashrc
```

### Step 3: Create Analysis Environment (same as Option 1)

```bash
conda create -n eclip python=3.10 -y
conda activate eclip

conda install -c bioconda -c conda-forge \
    sra-tools \
    fastqc \
    cutadapt \
    star \
    samtools \
    bedtools \
    umi_tools \
    -y

pip install pysam pybedtools pyBigWig pandas numpy matplotlib seaborn scipy jupyter
```

### Step 4: Access Windows Files from WSL

Your Windows C: drive is mounted at `/mnt/c/` in WSL:

```bash
cd /mnt/c/Users/s3ayy/OneDrive/Documents/GitHub/IFIT1-IFIT2-secondary-analysis
jupyter notebook --no-browser
```

---

## Tool Versions and Compatibility

Expected versions (as of January 2025):

- **Python**: 3.10+
- **sra-tools**: 3.0+
- **FastQC**: 0.12+
- **cutadapt**: 4.0+
- **STAR**: 2.7+
- **samtools**: 1.17+
- **bedtools**: 2.30+
- **umi_tools**: 1.1+

---

## Disk Space Requirements

Before starting the analysis, ensure you have adequate disk space:

- **Reference genome (hg19)**: ~3 GB
- **STAR index**: ~30 GB
- **GENCODE annotation**: ~50 MB
- **Raw FASTQ files** (per sample): ~400-800 MB
- **Aligned BAM files** (per sample): ~200-400 MB
- **Working space**: ~10 GB

**Total recommended**: At least **100 GB free space**

---

## Memory Requirements

- **STAR genome indexing**: 30-35 GB RAM
- **STAR alignment**: 30-35 GB RAM
- **Peak calling and analysis**: 8-16 GB RAM

If your system has less than 32 GB RAM, you may need to:
1. Use a compute cluster or cloud instance
2. Process samples one at a time
3. Use STAR's `--limitGenomeGenerateRAM` parameter

---

## Testing Your Installation

Create a test script to verify all tools:

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

echo -n "Python packages: "
python -c "import pysam, pybedtools, pyBigWig, pandas, numpy, matplotlib, seaborn, scipy; print('All packages installed')"

echo "=============================="
echo "Installation test complete!"
```

Run the test:

```bash
bash test_installation.sh
```

---

## Troubleshooting

### Issue: "STAR: command not found"

**Solution**: Make sure conda environment is activated:
```bash
conda activate eclip
```

### Issue: "Permission denied" when downloading from SRA

**Solution**: Configure SRA Toolkit:
```bash
vdb-config --interactive
# Navigate to "cache" and set a cache directory with write permissions
```

### Issue: Out of memory during STAR indexing

**Solution**: Use a smaller sjdbOverhang or process on a machine with more RAM:
```bash
# In the notebook, when building STAR index:
star_index = build_star_index(
    REFERENCE_DIR / 'hg19.fa.gz',
    REFERENCE_DIR / 'gencode.v19.annotation.gtf.gz',
    REFERENCE_DIR,
    threads=4,  # Reduce threads
    sjdb_overhang=74
)
```

### Issue: Slow downloads from SRA

**Solution**: Use prefetch to download first, then extract:
```bash
prefetch SRR31773513
fasterq-dump --outdir ./data/raw_fastq SRR31773513
```

---

## Alternative: Using Docker (Advanced)

If you're familiar with Docker, you can use a pre-built bioinformatics container:

```bash
# Pull a bioinformatics container
docker pull biocontainers/samtools:v1.9-4-deb_cv1

# Or use a complete environment like Jupyter + bioinformatics tools
docker pull jupyter/datascience-notebook
```

---

## Getting Help

If you encounter issues:

1. **Check tool documentation**: Most bioinformatics tools have extensive docs
2. **Bioconda issues**: https://github.com/bioconda/bioconda-recipes/issues
3. **SRA Toolkit**: https://github.com/ncbi/sra-tools/wiki
4. **STAR aligner**: https://github.com/alexdobin/STAR

---

## Next Steps

Once all tools are installed:

1. ✅ Open the Jupyter notebook
2. ✅ Run Section 1 (Environment Setup)
3. ✅ Run Section 2 (Sample Information) - verify the sample table displays
4. ✅ Proceed to download reference genome (Section 4)
5. ✅ Download and analyze eCLIP data
