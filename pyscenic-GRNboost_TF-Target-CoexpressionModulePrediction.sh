#!/bin/bash
#SBATCH -A kwcho_lab
#SBATCH -J GRN_pySCENIC
#SBATCH -p standard
#SBATCH --nodes=1 ##number of nodes to use
#SBATCH --cpus-per-task=30 ##number of cores to use per node
#SBATCH -error=slurm-%J.err ##Error file
#SBATCH --time=02:00:00  ##run time 2 min

cd /dfs5/bio/clarklh/Rcistarget/pyscenic

singleNucleiExpressionMatrix=/dfs10/bio/clarklh/Rcistarget/pyscenic/ThesisChap3/fuckyou/NeuralPlate_St10-13CombTraj_transposed_gene_expression_data.csv
TFlist=/dfs10/bio/clarklh/Rcistarget/pyscenic/iraTF_gene_names.csv
OUTPUT=/dfs10/bio/clarklh/Rcistarget/pyscenic/ThesisChap3/fuckyou/Inner_NE_subset_adata_with_stream2_graph_REAL.h5ad

pyscenic grn \
$singleNucleiExpressionMatrix \
$TFlist \
-o $OUTPUT \
--num_workers 30
