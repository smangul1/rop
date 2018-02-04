echo "$PWD/prerequisite_software/bwa/bwa mem  $PWD/db_human/BWA_version0.6.0/genome.fa" >rop_commands.txt
echo "$PWD/prerequisite_software/bwa/bwa mem  $PWD/db_human/BWA_version0.6.0/isoforms_GRCh38_Ensembl.fasta">>rop_commands.txt
echo "$PWD/prerequisite_software/tophat-2.1.0.Linux_x86_64/tophat2 -o ./ --fusion-search --keep-fasta-order --no-coverage-search $PWD/db_human/Bowtie2Index/genome">>rop_commands.txt
echo "$PWD/prerequisite_software/bwa/bwa mem $PWD/db_human/bacteria/bacteria.ncbi.february.3.2018.fasta">>rop_commands.txt
echo "$PWD/prerequisite_software/bwa/bwa mem $PWD/db_human/virus/virus.ncbi.february.3.2018.fasta">>rop_commands.txt
echo "$PWD/prerequisite_software/bwa/bwa mem $PWD/db_human/fungi/fungi.ncbi.february.3.2018.fasta">>rop_commands.txt
echo "$PWD/prerequisite_software/bwa/bwa mem $PWD/db_human/protozoa/protozoa.ncbi.february.3.2018.fasta">>rop_commands.txt
