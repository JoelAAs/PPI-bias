#!/bin/bash
#SBATCH -J Snakemake_main
#SBACTH -o /beegfs/scratch/ieo7513/logs/%j.out
#SBATCH -c 1 
#SBATCH -t 1-00:00:00
#SBATCH --mem=4G

folder_of_work=/beegfs/scratch/ieo7513
mkdir -p ${folder_of_work}/logs
mkdir -p /beegfs/scratch/ieo7513/apptainer

source /hpcnfs/home/ieo7513/miniforge3/etc/profile.d/conda.sh
conda activate snakemake
snakemake \
       	-s SnakeFile.smk \
	--workflow-profile profile/ \
	--use-conda \
	--use-apptainer --apptainer-args="--nv"

