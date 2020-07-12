#!/bin/bash
#SBATCH -t 12-00:00:00
#SBATCH -c 8
#SBATCH --mem=60G
#SBATCH -o expression-%j.out
#SBATCH -e expression-%j.out
#SBATCH --array=1-10

IDS=(
OG0000001.fa
OG0000002.fa
OG0000005.fa
OG0000006.fa
OG0000007.fa
OG0000016.fa
OG0000019.fa
OG0000021.fa
OG0000023.fa
OG0000026.fa
)

ID=${IDS[$SLURM_ARRAY_TASK_ID-1]}
echo $ID

iqtree -s INPUT -m LG+F+R4 -pre /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/OrthoFinder/Results_Feb16/WorkingDirectory/Alignments_ids/$ID -nt 8
