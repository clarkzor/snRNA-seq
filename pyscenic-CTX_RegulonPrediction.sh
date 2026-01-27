#!/bin/bash
#SBATCH -A kwcho_lab
#SBATCH -p standard
#SBATCH --nodes=1 ##number of nodes to use
#SBATCH --cpus-per-task=10 ##number of cores to use per node
#SBATCH -error=slurm-%J.err ##Error file
#SBATCH --time=02:00:00

cd /dfs5/bio/clarklh/Rcistarget/pyscenic/

Adjacencies=/dfs10/bio/clarklh/Rcistarget/pyscenic/ThesisChap3/fuckyou/InnerNNE_St10-13CombTraj_GRNboost.csv
RankedMotifRegions=/dfs10/bio/clarklh/Rcistarget/create_cisTarget_databases/Output/ATAC_St12/ATAC_St12_50kbGeneRegions.regions_vs_motifs.rankings.feather
MotifAnnotations=RcisTargetDB2.tbl
SingleNucleiExpressionMatrix=/dfs10/bio/clarklh/Rcistarget/pyscenic/ThesisChap3/fuckyou/InnerNNE_St10-13CombTraj_transposed_gene_expression_data.csv
Output=/dfs10/bio/clarklh/Rcistarget/pyscenic/ThesisChap3/fuckyou/InnerNNE_Subset_50kb_ATAC-Regions.motifs.csv


pyscenic ctx $Adjacencies $RankedMotifRegions --annotations_fname $MotifAnnotations --expression_mtx_fname $SingleNucleiExpressionMatrix  --output $Output --num_workers 15 --mode custom_multiprocessing --mask_dropouts
