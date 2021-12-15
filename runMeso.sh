# Script for run analysis on computational cluster managed by SLURM

#!/bin/bash

module load userspace/all
module load python3/3.6.3
source /home/pbrochet/snakemake/bin/activate

snakemake --jobs 50 --use-singularity --singularity-args "-B /scratch/pbrochet:/scratch/pbrochet" --snakefile /scratch/pbrochet/Projets/LensCardio/04.Papier_RNAseq_meth/Snakefile --configfile 02.Config/config.yaml --cluster-config 02.Config/clusterConfig.json --cluster 'sbatch -A {cluster.project} --job-name {cluster.job-name} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cores-number} --mem-per-cpu {cluster.mem-per-cpu} --partition {cluster.partition} --time {cluster.time} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --error {cluster.error} --output {cluster.output}'


# --jobs 1 : nombre jobs en parallèle Snakemake
