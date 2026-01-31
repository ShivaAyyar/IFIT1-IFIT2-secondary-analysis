#!/bin/bash
#SBATCH --job-name=eclip_reference
#SBATCH --output=logs/slurm_reference_%j.out
#SBATCH --error=logs/slurm_reference_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00

# SLURM batch script for reference genome setup
# This job requires ~64 GB RAM for STAR indexing

echo "Starting reference genome setup..."
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# Activate conda environment
source ~/.bashrc
conda activate eclip

# Run reference setup script
python scripts/01_setup_reference.py --threads $SLURM_CPUS_PER_TASK

echo ""
echo "Job completed at: $(date)"
