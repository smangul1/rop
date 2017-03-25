#!/bin/bash




echo "********************************************************************************
imrep_sra.sh is a script to run ImReP (without QC step of ROP) for sra file with mixture of mapped and unmapped reads.
It requires sratoolkit to be installed. Please add the path to sam-dump to this script manually.

ImReP is written by Igor Mandric and Serghei Mangul.

Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see: https://sergheimangul.wordpress.com/rop/
ROP Tutorial: https://github.com/smangul1/rop/wiki 
********************************************************************************"

if [ $# -lt 3 ]
then
#echo "[1] - file with samples located in the "
echo "[1] bam with mapped and unmapped reads, BAM file needs to be indexed!!!"
echo "[2] dir to save results of ROP"
echo "[3] addditional options for ROP. Should be in double quaotes. In case you don't need addditional options for ROP use \"\""
exit 1
fi

bam=$1
out=$2
options=$3

basename=$(echo  ${bam##*/} | awk -F ".bam" '{print $1}')


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SAMTOOLS=${DIR}/tools/samtools

#--------------------------------------------------------------------------------------
#UNMAPPED READS
mkdir $out




/u/home/h/harryyan//sratoolkit.2.5.0-centos_linux64/bin/sam-dump --unaligned $bam | SAMTOOLS view -bS - | SAMTOOLS bam2fq - > ${out}/${basename}_unmapped.fastq



#--------------------------------------------------------------------------------------
#MAPPED READS
#extract reads from IGH chr14, 106032614..107288051 GRCh37
#extract reads from IGK chr2, 89156874..89630436
#extract reads from IGL chr22, 22380474..23265085
#extract reads from TCRA chr14, 22090057..23021075
#extract reads from TCRB chr7, 141998851..142510972
# reads from TCRD inside TCRA
# #extract reads from TCRG chr7, 38279625..38407656



/u/home/h/harryyan//sratoolkit.2.5.0-centos_linux64/bin/sam-dump --aligned-region 14:106032614-107288051 $bam | SAMTOOLS view -bS - | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq - >${out}/${basename}_mapped_immune.fastq
/u/home/h/harryyan//sratoolkit.2.5.0-centos_linux64/bin/sam-dump --aligned-region 2:89156874-89630436 $bam | SAMTOOLS view -bS - | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq - >>${out}/${basename}_mapped_immune.fastq
/u/home/h/harryyan//sratoolkit.2.5.0-centos_linux64/bin/sam-dump --aligned-region 22:22380474-23265085 $bam | SAMTOOLS view -bS - | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq - >>${out}/${basename}_mapped_immune.fastq
/u/home/h/harryyan//sratoolkit.2.5.0-centos_linux64/bin/sam-dump --aligned-region 14:22090057-23021075 $bam | SAMTOOLS view -bS - | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq - >>${out}/${basename}_mapped_immune.fastq
/u/home/h/harryyan//sratoolkit.2.5.0-centos_linux64/bin/sam-dump --aligned-region 7:141998851-142510972 $bam | SAMTOOLS view -bS - | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq - >>${out}/${basename}_mapped_immune.fastq
/u/home/h/harryyan//sratoolkit.2.5.0-centos_linux64/bin/sam-dump --aligned-region 7:38279625-38407656 $bam | SAMTOOLS view -bS - | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq - >>${out}/${basename}_mapped_immune.fastq


cat ${out}/${basename}_mapped_immune.fastq ${out}/${basename}_unmapped.fastq > ${out}/${basename}_unmapped_plus_immune.fastq

#rm ${out}/${basename}_mapped_immune.fastq


#run ROP-ImReP for unmapped reads plus BCT/TCR reads
#Witpython ${DIR}/rop.py $options --f --immune --f ${out}/${basename}_unmapped_plus_immune.fastq ${out}/rop/

python ${DIR}//tools/imrep/imrep.py --fastq ${out}/${basename}_unmapped_plus_immune.fastq ${out}/${basename}


