#!/bin/bash
#SBATCH --job-name=eclip_full
#SBATCH --output=logs/slurm_full_%j.out
#SBATCH --error=logs/slurm_full_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mail-user=ssa163@case.edu
#SBATCH --mail-type=end


# SLURM batch script for running the complete pipeline sequentially
# This processes all samples one after another

echo "Starting complete eCLIP analysis pipeline..."
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# Activate conda environment
source ~/.bashrc
conda activate eclip

# Run complete pipeline
python scripts/run_pipeline.py --threads $SLURM_CPUS_PER_TASK

echo ""
echo "Pipeline completed at: $(date)"
