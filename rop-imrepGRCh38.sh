#!/bin/bash




echo "********************************************************************************
rop-imrepGRCh38.sh is now available under ROP release v1.0.7
rop-imrepGRCh38.sh is a script to run ROP-ImReP for bam file with mixture of mapped and unmapped reads.

This is the version for bam file with reads mapped to GRCh38 build of human genome.


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
#extract reads from IGH chr14, 105586437..106879844
#extract reads from IGK chr2, 88857361..90235368
#extract reads from IGL chr22, 22026076..22922913
#extract reads from TCRA chr14, 21621904..22552132
#extract reads from TCRB chr7, 142299011..142813287
# reads from TCRD inside TCRA
# #extract reads from TCRG chr7, 38240024..38368055




$SAMTOOLS view -bh ${bam} 14:105586437-106879844 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq - >${out}/${basename}_mapped_immune.fastq
$SAMTOOLS view -bh ${bam} 2:88857361-90235368 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq -  >>${out}/${basename}_mapped_immune.fastq
$SAMTOOLS view -bh ${bam} 22:22026076-22922913 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq -  >>${out}/${basename}_mapped_immune.fastq
$SAMTOOLS view -bh ${bam} 14:21621904-22552132 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq -  >>${out}/${basename}_mapped_immune.fastq
$SAMTOOLS view -bh ${bam} 7:142299011-1428132872 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq -    >>${out}/${basename}_mapped_immune.fastq
$SAMTOOLS view -bh ${bam} 7:38240024-38368055 | $SAMTOOLS view -bh -F 4 - | $SAMTOOLS bam2fq - >>${out}/${basename}_mapped_immune.fastq

cat ${out}/${basename}_mapped_immune.fastq ${out}/${basename}_unmapped.fastq > ${out}/${basename}_unmapped_plus_immune.fastq

rm ${out}/${basename}_mapped_immune.fastq ${out}/${basename}_unmapped.fastq


#run ROP-ImReP for unmapped reads plus BCT/TCR reads
python ${DIR}/rop.py $options --f --immune --f ${out}/${basename}_unmapped_plus_immune.fastq ${out}/rop/



