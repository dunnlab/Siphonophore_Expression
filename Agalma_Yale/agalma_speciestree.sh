#!/bin/bash
#SBATCH -t 6-00:00:00
#SBATCH -c 20
#SBATCH --mem=120G

module load miniconda

source activate agalma

export AGALMA_DB="/home/cm2477/scratch60/siphonophore2016/agalma-siphonophora.sqlite"

ID=SiphonophoraTree
SEED=87167481239

mkdir -p $ID
cd $ID

agalma speciestree --id $ID --outgroup Nematostella_vectensis --bootstrap 100 --seed $SEED
