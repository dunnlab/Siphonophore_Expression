#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH -N 8
#SBATCH -c 16
#SBATCH --mem=60G

module load miniconda

source activate agalma

export AGALMA_DB="/home/cm2477/scratch60/siphonophore2016/agalma-siphonophora.sqlite"

export BIOLITE_RESOURCES="threads=${SLURM_CPUS_ON_NODE},memory=${SLURM_MEM_PER_NODE}M"
export BIOLITE_HOSTLIST=$(hostlist -e -s, $SLURM_NODELIST)

ID=SiphonophoraTree

mkdir -p $ID
cd $ID

agalma homologize --id $ID
agalma multalign --id $ID
agalma genetree --id $ID
agalma treeinform --id $ID
agalma homologize --id $ID
agalma multalign --id $ID
agalma genetree --id $ID --bootstrap 100
agalma treeprune --id $ID
agalma multalign --id $ID
agalma supermatrix --id $ID --proportion 0.6
