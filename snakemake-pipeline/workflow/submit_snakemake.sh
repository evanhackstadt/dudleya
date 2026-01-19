#!/bin/bash
#
#SBATCH --job-name=Snakemake_Master
#SBATCH --output=out-snakemake-master-%j.out 
#
#SBATCH --partition=cmp
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1   # Cores for the Controller
#SBATCH --mem-per-cpu=4G    # Memory for the Controller
#SBATCH --time=24:00:00     # Max time for the Controller
#
#SBATCH --mail-user=ehackstadt@scu.edu
#SBATCH --mail-type=END
#
set -e
echo "Snakemake master wrapper script started at $(date)"
#
# 1. Load and activate conda environment where Snakemake is installed
module load Anaconda3
conda activate /WAVE/projects/whittalllab/conda_envs/genomics
#
# 2. Run the main Snakemake command
# All args specified in config/profile.yaml (including executor: slurm)
snakemake --profile "/WAVE/projects/whittalllab/dudleya/snakemake-pipeline/config/"