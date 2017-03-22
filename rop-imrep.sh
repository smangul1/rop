#!/bin/bash



echo "********************************************************************************
rop-imrep.sh is now available under ROP release v1.0.7
rop-imrep.sh is a script to run ROP-ImReP for bam file with mixture of mapped and unmapped reads.
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
$SAMTOOLS view -f 0x4 -bh  $bam | $SAMTOOLS bam2fq - > ${out}/${basename}_unmapped.fastq


#--------------------------------------------------------------------------------------
#MAPPED READS
#extract reads from IGH chr14, 106032614..107288051 GRCh37
#extract reads from IGK chr2, 89156874..89630436
#extract reads from IGL chr22, 22380474..23265085
#extract reads from TCRA chr14, 22090057..23021075
#extract reads from TCRB chr7, 141998851..142510972
# reads from TCRD inside TCRA
# #extract reads from TCRG chr7, 38279625..38407656



$SAMTOOLS view -bh ${bam} 14:106032614-107288051 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq - >${out}/${basename}_mapped_immune.fastq
$SAMTOOLS view -bh ${bam} 2:89156874-89630436 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq -  >>${out}/${basename}_mapped_immune.fastq
$SAMTOOLS view -bh ${bam} 22:22380474-23265085 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq -  >>${out}/${basename}_mapped_immune.fastq
$SAMTOOLS view -bh ${bam} 14:22090057-23021075 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq -  >>${out}/${basename}_mapped_immune.fastq
$SAMTOOLS view -bh ${bam} 7:141998851-142510972 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq -    >>${out}/${basename}_mapped_immune.fastq
$SAMTOOLS view -bh ${bam} 7:38279625-38407656 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq - >>${out}/${basename}_mapped_immune.fastq

cat ${out}/${basename}_mapped_immune.fastq ${out}/${basename}_unmapped.fastq > ${out}/${basename}_unmapped_plus_immune.fastq



#run ROP-ImReP for unmapped reads plus BCT/TCR reads
python ${DIR}/rop.py $options --f --immune --f ${out}/${basename}_unmapped_plus_immune.fastq ${out}/rop/

