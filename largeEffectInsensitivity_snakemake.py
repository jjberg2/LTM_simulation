#!/bin/bash
#SBATCH --job-name=largeEffectInsensitivity
#SBATCH --output=logs/largeEffectInsensitivity.out 
#SBATCH --error=logs/largeEffectInsensitivity.err 
#SBATCH --time=36:00:00 
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14 
#SBATCH --mem-per-cpu=4000

module load python

echo "SLURM_JOBID="$SLURM_JOBID

cat largeEffectInsensitivity_snakefile.py

snakemake --snakefile largeEffectInsensitivity_snakemake.py --cores all
