#!/bin/bash
#SBATCH -A kwcho_lab
#SBATCH -p standard
#SBATCH --nodes=1 ##number of nodes to use
#SBATCH --cpus-per-task=20 ##number of cores to use per node
#SBATCH --error=slurm-%J.err ##Error file

UPSTREAM=Xtrop_AllGenes_5kbUpstreamRegionsSt105.fasta
DOWNSTREAM=Xtrop_AllGenes_5kbDownstreamRegionsSt105.fasta
OUTPUT=Xtrop_TotalGTFmRNA_Up-Downstream_5kb_GeneRegions.fasta

seqkit concat $UPSTREAM $DOWNSTREAM -o $OUTPUT
