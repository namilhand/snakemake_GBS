# GBS pipeline

# 1. Setting environment

## 1) install conda

Conda provides a way to setup environment suited for GBS or other pipelines.

See “Quickstart install instructions part” of https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions

```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh

# inactivate conda being started by default
conda config --set auto_activate_base false
```

Or, you can use `mamba` , a faster implementation of conda, instead of conda

See:

- https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
- [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge)

## 2) Clone GBS git repository

```bash
mkdir /path/to/working/directory
cd /path/to/working/directory
git clone https://github.com/namilhand/snakemake_GBS.git .
```

Now all the scripts and data for GBS pipeline is cloned into your working direcotry

## 3) install GBS conda environment

The environment spec is specified in `environment.yaml` file. Use it for creating GBS environment as below.

```bash
# In the cloned directory, you can find "environment.yaml"
conda env create --name GBS --file environment.yaml
conda activate GBS
```

# 2. Run GBS pipeline

## 1) modify config.yaml

Change below accordingly in config.yaml

- THREADS: the number of threads you’re going to use
- sample_size: the number of sample size. Most of the time it is 96.
- lib_name: Name the set of GBS library. e.g., hcr1_set1
- P2_adapter_index: The number of P2 adapter used for library construction

## 2) Prepare input file

Create `raw` directory in your working directory.

Prepare the sequencing data (fastq.gz) in the `raw` directory. You can either copy the files or make symbolic link. Here I explain the simple way—copying the sequencing file into the `raw` directory.

The input file name should be formatted in: `P1_{P1 index}_P2{P2 index}_{read number}.fastq.gz`

Here `{read number}` should be either 1 or 2, which indicates read 1 or read 2 from paired-end sequencing.

There is not standardized way to rename and formatting the raw data. Below I provide a example.

```bash
#!/bin/bash

# The original input file name was formatted by Macrogen as like
# K9M21-[index1]_[1/2].fastq.gz
# P2-3 adapter was used in combination with P1 adapters for library construction.
# To re-format the name of file accordingly I use the code below.

mkdir -p newname

for f in *_1.fastq.gz; do
	read1=$f
	read2=${f%_1.fastq.gz}_2.fastq.gz
	P1=${read1#K9M21-}
	P1=${P1%_1.fastq.gz}
	P2=3

	new_read1=P1_${P1}_P2_${P2}_1.fastq.gz
	new_read2=P1_${P1}_P2_${P2}_2.fastq.gz

# Confirm if the code would work before re-naming it.
	echo "$read1 --- newname/$new_read1"
	echo "$read2 --- newname/$new_read2"

# Change name
# Always move the files with changed name into a different directory.
# Otherwise, the for loop will iterate using the re-named file
# and you will get the files that you won't be able to recognize the origin.

	mv $read1 newname/$new_read1
	mv $read2 newname/$new_read2
done
```

 

## 3) Run Snakefile

Now that all the files are prepared, you can simply run the pipeline by typing below

```bash
snakemake -p --cores [#THREADS]
# --printshellcmds, -p: Print out the shell commands that will be executed. (default: False)
# --cores: specify the number of threads that will be invested to the snakemake process.
```

# 3. Troubleshooting

- If the pipeline is interrupted anyhow, you can re-run from the interrupted point with `--rerun-incomplete` tag: `snakemake -p --cores [#THREADS] --rerun-incomplete`