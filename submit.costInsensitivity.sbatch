#!/bin/bash
#SBATCH --job-name=costInsensitivity
#SBATCH --output=logs/costInsensitivity.out 
#SBATCH --error=logs/costInsensitivity.err 
#SBATCH --time=36:00:00 
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14 
#SBATCH --mem-per-cpu=4000

module load python

echo "SLURM_JOBID="$SLURM_JOBID

cat prev_snakemake.py

snakemake --snakefile prev_snakemake.py --cores all
