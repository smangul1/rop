
dir=$1
ref=$2

while read line 
do
	name=$(echo $line | awk -F '.fastq' '{print $1}')
	echo "bwa mem -S -T 19 $ref ${name}.fq > ${name}_bwa.sam" >> run_${name}.sh 
	echo "perl CIRI.pl -I ${name}_bwa.sam -O ${name}_CIRI -F $ref -S"
done<$3


