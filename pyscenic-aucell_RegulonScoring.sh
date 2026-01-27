#!/bin/bash
#SBATCH -A kwcho_lab
#SBATCH -p standard
#SBATCH --nodes=1 ##number of nodes to use
#SBATCH --cpus-per-task=10 ##number of cores to use per node
#SBATCH -error=slurm-%J.err ##Error file
#SBATCH --time=02:00:00

cd /dfs10/bio/clarklh/Rcistarget/pyscenic/

pyscenic aucell /dfs10/bio/clarklh/Rcistarget/pyscenic/ThesisChap3/InnerNNE_St10-13CombTraj_transposed_gene_expression_data.csv /dfs10/bio/clarklh/Rcistarget/pyscenic/ThesisChap3/fuckyou/InnerNNE_Subset_50kb_ATAC-Regions_PPR.motifs.csv --output /dfs10/bio/clarklh/Rcistarget/pyscenic/ThesisChap3/fuckyou/InnerNNE_subset_Adj005_AllGenesMatrix_50kb_ATAC-Regions_aucell_PPR.csv --num_workers 10
