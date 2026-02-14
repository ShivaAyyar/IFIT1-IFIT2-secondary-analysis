#!/bin/bash
#SBATCH --job-name=fix_dedup
#SBATCH --output=logs/slurm_fix_dedup_%j.out
#SBATCH --error=logs/slurm_fix_dedup_%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mail-user=ssa163@case.edu
#SBATCH --mail-type=END,FAIL

# SLURM job to fix empty deduplicated BAM files
# This re-runs umi_tools dedup on the STAR output BAM files

echo "=============================================="
echo "Fix Empty Deduplicated BAM Files"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# Change to project directory
cd /mnt/vstor/SOM_CCCC_JGS25/sayyar/IFIT1-IFIT2-secondary-analysis

# Activate conda environment
source ~/.bashrc
conda activate eclip

# Run the fix script
chmod +x scripts/fix_dedup_bams.sh
bash scripts/fix_dedup_bams.sh

echo ""
echo "Job completed at: $(date)"
