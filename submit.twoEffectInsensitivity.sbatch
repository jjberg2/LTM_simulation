#!/bin/bash
#SBATCH --job-name=twoEffectInsensitivity 
#SBATCH --output=logs/twoEffectInsensitivity.out 
#SBATCH --error=logs/twoEffectInsensitivity.err 
#SBATCH --time=36:00:00 
#SBATCH --partition=broadwl
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=14 
#SBATCH --mem-per-cpu=4000

module load python

echo "SLURM_JOBID="$SLURM_JOBID

cat twoEffectInsensitivity_snakemake.py

snakemake --snakefile twoEffectInsensitivity_snakemake.py --cores all
