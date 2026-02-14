#!/bin/bash
#SBATCH --job-name=clipper_peaks
#SBATCH --output=logs/slurm_peaks_%A_%a.out
#SBATCH --error=logs/slurm_peaks_%A_%a.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --array=0-5
#SBATCH --mail-user=ssa163@case.edu
#SBATCH --mail-type=END,FAIL

# SLURM array job script for peak calling on all 6 IP samples in parallel
# Each array task processes one IP sample with CLIPper
#
# Usage:
#   sbatch slurm/call_peaks.sbatch
#
# Array indices map to IP samples:
#   0: IFIT2_IFIT3_FLAG_IP
#   1: IFIT2_IFIT3_HA_IP
#   2: IFIT2_alone_FLAG_IP
#   3: IFIT3_alone_HA_IP
#   4: WT_FLAG_IP
#   5: WT_HA_IP

echo "=============================================="
echo "CLIPper Peak Calling - SLURM Array Job"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Start time: $(date)"
echo ""

# Change to project directory
cd /mnt/vstor/SOM_CCCC_JGS25/sayyar/IFIT1-IFIT2-secondary-analysis

# Activate conda environment
source ~/.bashrc
conda activate eclip

# Verify Singularity is available
if ! command -v singularity &> /dev/null; then
    echo "ERROR: Singularity not found. Loading module..."
    module load singularity 2>/dev/null || true
fi

# Verify container exists
if [ ! -f "clipper_5d865bb.sif" ]; then
    echo "ERROR: Singularity container not found: clipper_5d865bb.sif"
    echo "Please run: singularity pull docker://brianyee/clipper:5d865bb"
    exit 1
fi

# Run peak calling for this sample index
echo "Running peak calling for sample index: $SLURM_ARRAY_TASK_ID"
echo ""

python scripts/04_call_peaks.py --sample-index $SLURM_ARRAY_TASK_ID

EXIT_CODE=$?

echo ""
echo "=============================================="
echo "Job completed at: $(date)"
echo "Exit code: $EXIT_CODE"
echo "=============================================="

exit $EXIT_CODE
