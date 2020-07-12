#!/bin/bash
#SBATCH -J
#SBATCH -t 92:00:00
#SBATCH -c 8
#SBATCH --mem=60G
#SBATCH -o expression-%j.out
#SBATCH -e expression-%j.out
#SBATCH --array=1-4

IDS=(
OG0000000.fa
OG0000008.fa
OG0000014.fa
OG0000033.fa
)

ID=${IDS[$SLURM_ARRAY_TASK_ID-1]}
echo $ID

mafft --localpair --maxiterate 1000 --anysymbol --thread 8 /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/OrthoFinder/Results_Feb16/WorkingDirectory/Sequences_ids/$ID > /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/OrthoFinder/Results_Feb16/WorkingDirectory/Alignments_ids/$ID
