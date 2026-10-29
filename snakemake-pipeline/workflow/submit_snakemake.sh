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
#SBATCH --time=96:00:00     # Max time for the Controller
#
#SBATCH --mail-user=ehackstadt@scu.edu
#SBATCH --mail-type=END
#
set -e
echo "Snakemake master wrapper script started at $(date)"
#
# Store variables for path to specific config file, or additional flags, if passed
config=${1:-''}
flags=${2:-''}
# Store variable for profile (static)
profile="/WAVE/projects/whittalllab/dudleya/snakemake-pipeline/config/"
base_config="/WAVE/projects/whittalllab/dudleya/snakemake-pipeline/config/project.yaml"
echo "Using profile: $profile"
echo "Using config file: $config"
echo "Adding flags: $flags"
#
# 1. Load and activate conda environment where Snakemake is installed
module load Anaconda3
conda activate /WAVE/projects/whittalllab/conda_envs/genomics
#
# 2. Run the main Snakemake command
# Rest of the args are specified in config/profile.yaml (including executor: slurm)
snakemake --profile $profile --configfile $base_config $config $flags