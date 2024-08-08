#!/bin/bash
#SBATCH -J isoseq_algn # Job name
#SBATCH -o job.%j.out # Name of stdout output file (%j expands to %jobId)
#SBATCH -N 1 # Total number of nodes requested
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8


source /home/azureuser/progs/mRNA_Seq_Analysis/modules/aligners/alignLongRead.bash
GENOME_IDX="/home/azureuser/projects/HBV/alignments/2x_genome/88A010_2700-3182_1-3182_1-2100_2x_genome.mmi"

INPUT_FASTA=$1

align_long_reads $INPUT_FASTA "PBMM2" $GENOME_IDX ;