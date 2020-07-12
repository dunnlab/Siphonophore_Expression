#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -c 12
#SBATCH --mem=24G
#SBATCH --array=1-3
sleep $((SLURM_ARRAY_TASK_ID*60))

set -e

module load miniconda

source activate agalma

export AGALMA_DB="/home/cm2477/scratch60/siphonophore2016/agalma-siphonophora.sqlite"

export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"

IDS=(
	dbEST_CLYHEM
	NCBI_HYDMAG
	JGI_NEMVEC
)

ID=${IDS[$SLURM_ARRAY_TASK_ID-1]}
echo $ID

agalma import --id $ID
agalma translate --id $ID
