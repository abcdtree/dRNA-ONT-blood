#!/bin/bash
#SBATCH --job-name=bash_command
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mem=10G
#SBATCH --time=0-01:00:00
#SBATCH -p sapphire

module load SAMtools

for b in `ls *.bam | cut -f 1 -d .`
do
    samtools view ${b}.bam | awk '/pt:i/{print $1,length($10),$NF}' | sed 's/pt:i://g'  > ${b}.polyA.stats
done
