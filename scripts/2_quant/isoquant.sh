#!/bin/bash
#SBATCH --job-name=isoquant-devika
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH -p sapphire


reference=GRCh38_chromosome_only.fa
gff=$1
output=$2

isoquant.py --reference ${reference} \
    --genedb ${gff} \
    --fastq *.fastq \
    --data_type nanopore -o ${output} -t 20 --sqanti_output --count_exons
