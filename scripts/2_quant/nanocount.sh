#!/bin/bash
#SBATCH --job-name=minimap2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --account=punim1068
#SBATCH --mail-type=ALL
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH -p sapphire


module load SAMtools/1.16.1
module load minimap2/2.24
module load Python/3.10.4

fastq_home=$1
transcriptome_ref=$2
mmi_ref=trans.mmi

minimap2 -x map-ont -d /${mmi_ref} ${transcriptome_ref}
for b in `ls ${fastq_home} | cut -f 1 -d _`
do
    minimap2 -t 20 -ax splice -uf -k14 -N 10  ${mmi_ref} ${fastq_home}/${b}_fastq_pass.fastq | samtools sort -o ${b}.bam
    samtools sort -o ${b}.bam ${b}.sam
    samtools index ${b}.bam
    NanoCount -i ${b}.bam --extra_tx_info -o ${b}_counts.tsv
done
