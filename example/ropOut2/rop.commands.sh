/u/home/galaxy/collaboratory/serghei/code/rop/tools/blastn -task megablast -index_name /u/home/galaxy/collaboratory/serghei/code/rop//db_human/rRNA/rRNA -use_index true -query /u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut2/QC/unmappedExample_lowQ.fa -db /u/home/galaxy/collaboratory/serghei/code/rop//db_human/rRNA/rRNA -outfmt 6 -evalue 1e-05 >/u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut2/QC/unmappedExample_rRNA_blastFormat6.csv 2>log_megablast_rRNA.log
2. Remapping to reference...
/u/home/s/serghei/collab/code/rop-fork1/prerequisite_software/bwa/bwa mem  /u/home/s/serghei/collab/code/rop-fork1/db_human/BWA_version0.6.0/genome.fa /u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut2/QC/unmappedExample_after_rRNA.fasta 2>>/u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut2/lostReads/unmappedExample_bowtieWG.log | /u/home/galaxy/collaboratory/serghei/code/rop/tools/samtools view -SF4 - >/u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut2/lostReads/unmappedExample_genome.sam
/u/home/s/serghei/collab/code/rop-fork1/prerequisite_software/bwa/bwa mem  /u/home/s/serghei/collab/code/rop-fork1/db_human/BWA_version0.6.0/isoforms_GRCh38_Ensembl.fasta /u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut2/QC/unmappedExample_after_rRNA.fasta 2>>/u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut2/lostReads/unmappedExample_bowtieTR.log | /u/home/galaxy/collaboratory/serghei/code/rop/tools/samtools view -SF4 -  >/u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut2/lostReads/unmappedExample_transcriptome.sam
3. Mapping repeat sequences...
