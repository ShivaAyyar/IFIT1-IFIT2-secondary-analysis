#!/bin/bash
#SBATCH --job-name=eclip_sample
#SBATCH --output=logs/slurm_sample_%A_%a.out
#SBATCH --error=logs/slurm_sample_%A_%a.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --array=0-11
#SBATCH --mail-user=ssa163@case.edu
#SBATCH --mail-type=end


# SLURM array job script for processing all 12 samples in parallel
# Each array task processes one sample

echo "Starting sample processing..."
echo "Job ID: $SLURM_JOB_ID"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# Activate conda environment
source ~/.bashrc
conda activate eclip

# Run sample processing script with array index
python scripts/03_process_samples.py \
    --sample-index $SLURM_ARRAY_TASK_ID \
    --threads $SLURM_CPUS_PER_TASK

echo ""
echo "Job completed at: $(date)"
