#!/bin/bash
folder_of_work=/beegfs/scratch/ieo7513
#SBATCH -J Snakemake_main
#SBACTH -o ${folder_of_work}/logs/%j.out
#SBATCH -c 1 
#SBATCH -t 1-00:00:00
#SBATCH --mem=4G

mkdir -p ${folder_of_work}/logs
mkdir -p /beegfs/scratch/ieo7513/apptainer

conda activate snakemake
snakemake \
       	-s SnakeFile.smk \
	--workflow-profile profile/ \
	--use-conda \
	--use-apptainer --apptainer-args="--nv" -n

