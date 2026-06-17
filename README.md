# Dudleya Genomics Pipeline

Dr. Justen Whittall, Evan Hackstadt, Karina Martinez, Dante Cable
Santa Clara University (SCU)

README last updated: 6/16/2026

## About

This repo contains the bioinformatics pipeline used to study the conservation genomics of *Dudleya setchellii*, an endangered succulent endemic to Santa Clara county, California. The research project is led by Dr. Justen Whittall (Biology Department, Santa Clara University). Preliminary results were presented at the 2026 West Coast Biological Sciences Undergraduate Research Conference. Journal publication is anticipated to be forthcoming.
The results should be reproducible, but note that the pipeline is specifically designed for use on SCU's WAVE High Performance Cluster (HPC) and our specific conda environments and datasets.


## Pipeline Overview

This pipeline takes a collection of whole-genome DNA sequences, cleans them, maps them to a reference genome, performs population analyses, and produces three primary outputs: a PCA plot, a genetic distance matrix, and a phylogenetic tree. It is robust to low coverage sequencing data and preserves uncertainty for the PCA and distance matrix.

First, each sample (forward and reverse read) is processed as follows:
1. fastp - filtering & QC
2. BWA MEM - alignment to reference genome --> .bam
3. Samtools - sort, index, and get QC report

Then, all samples (.bam files) are aggregated for analysis:
4. ANGSD - population analysis --> .beagle, .mafs, .bcf (parallel chromosomal scatter-gather optimizes runtime)
5. bcftools - convert .bcf --> .vcf --> .phy
6. IQTREE - produce phylogenetic tree
7. PCAngsd - produce PCA plot
8. ngsDist - produce genetic distance matrix
9. Custom aggregate QC report


## Pipeline Architecture

Pipeline Files:
* `Snakefile` — contains the bulk of the **pipeline**. It is a list of steps (terminal commands / scripts) to process and analyze the data.
* Two **config files** needed to run the pipeline.
  * The profile (`config/config.yaml`). Contains parameters telling Snakemake to use SLURM and specifying default resources. It MUST be named config.yaml and shouldn't need to be modified.
  * A config file (e.g. `config/samples.yaml`). Paths to the raw input files and other files used by the pipeline. We can create different configs (e.g. `config/biol173-multi-pop.yaml`) to specify different input files or datasets.

Key Scripts:
* `update_sample_config.py` — utility script providing an easy way to create or modify config files.
* `submit_snakemake.sh` — batch file that allows us to run the pipeline in the background (as a master SLURM job).


## Running the Pipeline

### Snakemake Dry Run

Snakemake allows you to perform a **dry run** using the `-n` flag. This simply prints the steps/jobs without actually running them. It's useful for a sanity check before a real run and should be done directly on the command line as follows.

First, activate your conda environment containing Snakemake. Then run the command:
```
snakemake -np --profile {path/to/profile_directory/} --configfile {path/to/congfigfile.yaml}
```
Example usage:
```
snakemake -np --profile config/ --configfile config/samples.yaml
```

Removing the `-n` flag would run the pipeline directly on the command line, but this is inconvenient for long jobs. Instead, use SLURM.

### Snakemake Run on SLURM

To run the Snakemake pipeline in the background (as a SLURM job), use the `submit_snakemake.sh` script:
```
sbatch snakemake/submit_snakemake.sh {path/to/configfile}
```
Example usage:
```
sbatch snakemake/submit_snakemake.sh config/biol173-multi-pop.yaml
```
And that's it! To check status, use `squeue` to look for your jobs, `tail <out_file>` to check logs, or `ls <dirs_being_created>` to check outputs.
The pipeline stores outputs in a subdirectory (named per the config file) within `results/`.
View plots/results by copying them to your machine using the linux `scp` command, or view them through a virtual desktop.


## Creating a New Config
Let's say we have **new samples** in a folder `data/` that we want to process.
First we need to create a new `.yaml` to point to our samples.

We can use the utility script, `update_sample_config.py`, to do this easily without editing the file directly.

The script takes a few arguments and has optional flags. To see options, run:
```
python scripts/update_sample_config.py --help
```
Example usage:
```
python scripts/update_sample_config.py data/ snakemake/samples_config/ --filename new_config
```
The script should print what it's doing to the terminal. Note that this script CANNOT change paths to the ref and anc genomes, so these must be changed manually if needed.

Once the config is successfully created, run the pipeline on it as specified above.
