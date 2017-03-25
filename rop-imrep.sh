#!/bin/bash


source $(dirname $0)/argparse.bash
argparse "$@" <<EOF || exit 1


parser.add_argument('infile')
parser.add_argument('out')

parser.add_argument('-r', '--release', default="0", type=int,
help='version of human genome release. Default is hg19(GRCh37). Use RELEASE=1 to select GRCh38 ')

parser.add_argument('-o', '--organism', default="0", type=int,
help='Organism. Default is human. Use ORGANISM=1 to select mouse')

parser.add_argument('-i', '--input', default="0", type=int,
help='Input format. Default is bam. Use INPUT=1 to select sra input format. It requires sratoolkit to be installed. Please add the path to sam-dump to this script manually. ')

EOF

echo required infile: "$INFILE"
echo required outfile: "$OUTFILE"
echo optional: "$RELEASE"



exit 1

echo "********************************************************************************
rop-imrep.sh is now available starting from ROP release v1.0.7

rop-imrep.sh is a script to run ROP-ImReP for bam file with mixture of mapped and unmapped reads. It assumes reads are mapped to hg19(GRCh37) version of human genome.
ImReP is written by Igor Mandric and Serghei Mangul.
Released under the terms of the General Public License version 3.0 (GPLv3)
For more details see: https://sergheimangul.wordpress.com/rop/
ROP Tutorial: https://github.com/smangul1/rop/wiki

********************************************************************************"

if [ $# -lt 4 ]
then
#echo "[1] - file with samples located in the "
echo "[1] bam with mapped and unmapped reads, BAM file needs to be indexed!!!"
echo "[2] dir to save results of ROP"
echo "[3] addditional options for ROP. Should be in double quaotes. In case you don't need addditional options for ROP use \"\""
echo "[4] prefix on how chr names are encoded in the BAM file. In case they are encoded as (chr1, chr2, etc), use chr, otherwise use \"\""
exit 1
fi

bam=$1
out=$2
options=$3
chr=$4

basename=$(echo  ${bam##*/} | awk -F ".bam" '{print $1}')


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SAMTOOLS=${DIR}/tools/samtools

#--------------------------------------------------------------------------------------
#UNMAPPED READS
mkdir $out
SAMTOOLS view -f 0x4 -bh  $bam | SAMTOOLS bam2fq - > ${out}/${basename}_unmapped.fastq


#--------------------------------------------------------------------------------------
#MAPPED READS
#extract reads from IGH chr14, 106032614..107288051 GRCh37
#extract reads from IGK chr2, 89156874..89630436
#extract reads from IGL chr22, 22380474..23265085
#extract reads from TCRA chr14, 22090057..23021075
#extract reads from TCRB chr7, 141998851..142510972
# reads from TCRD inside TCRA
# #extract reads from TCRG chr7, 38279625..38407656




SAMTOOLS view -bh ${bam} ${chr}14:106032614-107288051 | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq - >${out}/${basename}_mapped_immune.fastq
SAMTOOLS view -bh ${bam} ${chr}2:89156874-89630436 | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq -  >>${out}/${basename}_mapped_immune.fastq
SAMTOOLS view -bh ${bam} ${chr}22:22380474-23265085 | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq -  >>${out}/${basename}_mapped_immune.fastq
SAMTOOLS view -bh ${bam} ${chr}14:22090057-23021075 | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq -  >>${out}/${basename}_mapped_immune.fastq
SAMTOOLS view -bh ${bam} ${chr}7:141998851-142510972 | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq -    >>${out}/${basename}_mapped_immune.fastq
SAMTOOLS view -bh ${bam} ${chr}7:38279625-38407656 | SAMTOOLS view -bh -F 4 - | SAMTOOLS bam2fq - >>${out}/${basename}_mapped_immune.fastq

cat ${out}/${basename}_mapped_immune.fastq ${out}/${basename}_unmapped.fastq > ${out}/${basename}_unmapped_plus_immune.fastq



#run ROP-ImReP for unmapped reads plus BCT/TCR reads
python ${DIR}/rop.py $options --f --immune --f ${out}/${basename}_unmapped_plus_immune.fastq ${out}/rop/




