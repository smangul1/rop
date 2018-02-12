echo "$PWD/tools/bwa mem  $PWD/db_human/BWA.index/genome.fa" >rop.commands.txt
echo "$PWD/tools/bwa mem  $PWD/db_human/BWA.index/isoforms_GRCh38_Ensembl.fasta">>rop.commands.txt
echo "$PWD/prerequisite_software/tophat-2.1.0.Linux_x86_64/tophat2 -o ./ --fusion-search --keep-fasta-order --no-coverage-search $PWD/db_human/Bowtie2Index/genome">>rop.commands.txt
echo "$PWD/tools/bwa mem $PWD/db_human/bacteria/bacteria.ncbi.february.3.2018.fasta">>rop.commands.txt
echo "$PWD/tools/bwa mem -a $PWD/db_human/viral/viral.ncbi.february.3.2018.fasta">>rop.commands.txt
echo "$PWD/tools/bwa mem -a $PWD/db_human/viral.vipr/NONFLU_All.fastq">>rop.commands.txt
echo "$PWD/tools/bwa mem -a $PWD/db_human/fungi/fungi.ncbi.february.3.2018.fasta">>rop.commands.txt
echo "$PWD/tools/bwa mem -a $PWD/db_human/protozoa/protozoa.ncbi.february.3.2018.fasta">>rop.commands.txt
