#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH -N 8
#SBATCH -c 16
#SBATCH --mem=60G

module load IQ-TREE

/home/cm2477/project/OrthoFinder/orthofinder -f /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences -t 16 -s /home/cm2477/scratch60/siphonophore2016/OrthoFinder/speciestree.tree -M msa -T iqtree 
