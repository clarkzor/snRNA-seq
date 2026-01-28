#!/bin/bash
#SBATCH -A kwcho_lab
#SBATCH -p standard
#SBATCH --nodes=1 ##number of nodes to use
#SBATCH --cpus-per-task=35 ##number of cores to use per node
#SBATCH --error=slurm-%J.err ##Error file

cd /dfs10/bio/clarklh/Rcistarget/create_cisTarget_databases/

FASTA=/dfs10/bio/clarklh/ATAC/St12/MACS2/ATAC_St12_closestBed_50kb_genes_ForFasta-PythonCombine.fasta
CBmotifDB=/dfs10/bio/clarklh/Rcistarget/create_cisTarget_databases/CBDB/
CBmotifList=/dfs10/bio/clarklh/Rcistarget/create_cisTarget_databases/St105IraTotalTF_MotifList.csv
OUTPUT=/dfs10/bio/clarklh/Rcistarget/create_cisTarget_databases/Output/ATAC_St12/ATAC_St12_50kbGeneRegions

python create_cistarget_motif_databases.py -f $FASTA -M $CBmotifDB -m $CBmotifList -o $OUTPUT -t 35 
