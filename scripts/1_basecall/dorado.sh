#!/bin/bash
#SBATCH --job-name=basecalling
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --time=40:00:00

module load SAMtools/1.16.1
module load minimap2/2.24
module load Python/3.10.4
module load dorado/0.5.2-CUDA-11.7.0

dorado_home=$1
model_home=$2

data_home=$3
output_path=$4

for b in `ls ${data_home}`
do
    cd ${output_path}
    pod5 convert fast5 ${data_home}/${b}/*/*/fast5_pass/*.fast5 --output ${b}.pod5
    ${dorado_home}/dorado basecaller ${model_home}/rna002_70bps_hac@v3 ${b}.pod5 --estimate-poly-a  > ${b}.bam
    samtools bam2fq ${b}.bam > ${b}.fq
    #rm ${b}.bam
    rm ${b}.pod5
done
