#!/bin/bash
#SBATCH -A kwcho_lab
#SBATCH -p standard
#SBATCH -J SeqMerge
#SBATCH --nodes=1 ##number of nodes to use
#SBATCH --cpus-per-task=10 ##number of cores to use per node
#SBATCH -error=slurm-%J.err ##Error file


INPUT=/dfs10/bio/clarklh/ATAC/St12/MACS2/ATAC_St12_closestBed_50kb_genes_ForFasta.fasta
OUTPUT=/dfs10/bio/clarklh/ATAC/St12/MACS2/ATAC_St12_closestBed_50kb_genes_ForFasta-PythonCombine.fasta

python SeqMerge.py $INPUT > $OUTPUT



