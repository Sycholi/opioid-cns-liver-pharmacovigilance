#!/bin/bash
#SBATCH --job-name=faers_gate
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
cd $SLURM_SUBMIT_DIR
/usr/bin/Rscript gate_faers.R
