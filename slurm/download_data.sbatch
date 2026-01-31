#!/bin/bash
#SBATCH --job-name=eclip_download
#SBATCH --output=logs/slurm_download_%j.out
#SBATCH --error=logs/slurm_download_%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --mail-user=ssa163@case.edu
#SBATCH --mail-type=end


# SLURM batch script for downloading FASTQ files from SRA
# Downloads all 12 samples (~400-800 MB each)

echo "Starting FASTQ download from SRA..."
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# Activate conda environment
source ~/.bashrc
conda activate eclip

# Download all samples
python scripts/02_download_data.py --all

echo ""
echo "Download completed at: $(date)"
