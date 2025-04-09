#!/bin/bash
#SBATCH --job-name=sqanti3
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mem=200G
#SBATCH --time=24:00:00


module load Apptainer/1.2.3
image_path=sqanti3_latest.error.sif
data_home=$1

my_sequence=${data_home}/$2
anno=${data_home}/gencode.v35.annotation.gtf
gff=${data_home}/gencode.v35.annotation.gff3
output_path=${data_home}/$3
extra_file=${data_home}/SQANTI3-extra/Human
CAGE=human.refTSS_v3.1.hg38.bed
ployA=${extra_file}/human.polyA.txt
tappAS=${extra_file}/tappAS/Homo_sapiens_GRCh38_Ensembl_86.gff3
genome=${data_home}/GRCh38_chromosome_only.fa


apptainer exec --bind ${data_home}:${data_home} ${image_path} sqanti3_qc.py ${my_sequence} ${anno} ${genome}  \
--force_id_ignore -t 20  -o isoquant -d ${output_path} --CAGE_peak ${CAGE} \
--polyA_motif_list ${ployA} --gff3 ${tappAS}

apptainer exec --bind ${data_home}:${data_home} ${image_path} sqanti3_filter.py \
  ml ${output_path}/isoquant_classification.txt -o isoquant_filter -d ${output_path} -i 60 -j 0.7
