# Dudleya Genomics Pipeline

Dr. Justen Whittall, Evan Hackstadt, Karina Martinez, Dante Cable

## Instructions for Current Snakemake Pipeline
Last updated: 1/6/26

### The Architecture

* `Snakefile` — contains the bulk of the **pipeline**. It is a list of steps (terminal commands / scripts) to process the data.
* The Snakefile depends on two **config files** that tell it what to do:
  * `profiles/config.yaml` — contains parameters telling Snakemake to use SLURM and specifying default resources. We do NOT normally change this.
  * `samples_config/config.yaml` — paths to the ref genome, anc genome, and a dictionary mapping sample names to the paths to each raw file (R1 and R2). We MAY change this.
  * `update_sample_config.py` — utility script providing an easy way to update samples config file.
* `submit_snakemake.sh` — batch file that allows us to run the pipeline in the background (as a master SLURM job).

### Running the Pipeline

#### Update sample config
Let's say we have **new samples** in a folder `data/` that we want to process.
First we need to udpate `samples_config/config.yaml` to point to our samples.

You can use the utility script, `update_sample_config.py`, to do this easily without editing the file directly.

The script takes a few arguments and has optional flags. To see options, run:
```
python scripts/update_sample_config.py --help
```

Example usage:
```
python scripts/update_sample_config.py data/ snakemake/samples_config/
```

The script should print what it's doing to the terminal. Note that this script CANNOT change paths to the ref and anc genomes, so these must be changed manually if needed.

#### Run Snakemake
Now we can run the pipeline.

To run Snakemake in the background (as a SLURM job):
```
sbatch snakemake/submit_snakemake.sh
```

And that's it! To check status, use `squeue` or `tail <output_file>` or `ls <dirs_being_created>`.
The pipeline creates a master `results/` directory containing subdirectories for relevant outputs from each step.

To run Snakemake directly (not recommended since this requires you to stay logged in until it finishes):
```
snakemake --profile profiles/config.yaml
```

#### Currently: Visualization is Manual
We hope to integrate visualization (PCA plotting) into the pipeline, but currently it must be done manually.

Once Snakemake finishes, you should have `results/pca/` containing `population.cov` and `population.info` files.
We need to pass these to `pcangsd_visualize.py` to create a plot. This script requires the visualization conda env.

Example usage:
```
conda activate /WAVE/projects/whittalllab/conda_envs/visualization
python scripts/pcangsd_visualize.py snakemake/results/pca/
```
The PCA plot should be saved to the directory provided. You can then `scp` it onto your personal computer to view.
