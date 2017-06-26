
For a quick start, a toy example with 2508 unmapped reads is distributed with the ROP package. Please note that selected reads are randomly selected from a normal skin (SRR1146076) RNA-Seq sample and might not represent the typical reads of RNA-Seq experiment. The reads are provided for demonstration purposes and can be accessed under ROP directory. Instruction for the analysis of the full RNA-Seq sample (SRR1146076) are provided [here](https://github.com/smangul1/rop/wiki/ROP-analysis:--one-RNA-Seq-sampe).

```
/rop/example/unmappedExample.fastq
```




Please make sure that the basic unix commands (wget, python, perl) are available on the cluster.  Please install ROP first. The instructions how to install ROP are provided [here](https://github.com/smangul1/rop/wiki/How-to-install-ROP%3F) 

The first operation consists of navigating to ROP directory and creating a subdirectory for storing the training data. 

```
cd rop
mkdir tutorial
cd tutorial
mkdir data
```

Now, download the mapped and unmapped reads from RNA-Seq

```
wget https://googledrive.com/host/0B_NUyiE86yDwaUxoVjhlSjN5SkE/skinExample.tar
tar -xvf skinExample.tar 
```

Now you have skin RNA-Seq sample under the dat directory. Remember that you are located under data directory. Go back to rop directory and run `rop.py` from there. 

ROP is an intensive pipeline requiring substantial amount of computations resources. Thus we don't recommend to run ROP from login nodes. Please check the policy of your cluster, from where to run the ROP pipeline. For hoffman2 (UCLA cluster) read the policy [here] (http://www.hoffman2.idre.ucla.edu/computing/interactive-session/). 
For hoffman2 run this command to get an interactive session : 

```
qrsh -l h_data=16G,h_rt=4:00:00
```


Alternatively, you may consider runing the ROP analysis via qsub. More details are available at [https://github.com/smangul1/rop/wiki/ROP-analysis-via-qsub](https://github.com/smangul1/rop/wiki/ROP-analysis-via-qsub)

ROP requires two mandatory command line arguments, i.e. (1) the unmapped reads and (2) the directory to save the results of ROP.

```
usage: python rop.py [-h] [--qsub] [--qsubArray] [--maui]
                     [--organism ORGANISM] [--b] [--gzip] [--skipLowq]
                     [--skipQC] [--skipPreliminary] [--repeat] [--circRNA]
                     [--immune] [--microbiome] [--outGz] [--rezip] [--clean]
                     [--quiet] [--dev] [--nonReductive] [--f]
                     unmappedReads dir
```

To test ROP for the small training data (toy example) use the following command under the ROP directory, where results will be saved to example/ropOut/ directory. Make

```
python rop.py example/unmappedExample.fastq example/ropOut/
```

You expect the following output of the ROP pipeline:

```
*********************************************
ROP (version 1.0.1) is a computational protocol aimed to discover the source of all reads, originated from complex RNA molecules, recombinant antibodies and microbial communities. Written by Serghei Mangul (smangul@ucla.edu) and Harry Taegyun Yang (harry2416@gmail.com), University of California, Los Angeles (UCLA). (c) 2016. Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see:
https://sergheimangul.wordpress.com/rop/
*********************************************
Processing 2510 unmapped reads
1. Quality Control...
--filtered 2193 low quality reads
--filtered 2 low complexity reads (e.g. ACACACAC...)
--filtered 22 rRNA reads
In toto : 2217 reads failed QC and are filtered out
2. Remaping to human references...
--identified 6 lost human reads from unmapped reads. Among those: 4 reads with 0 mistmathes; 2 reads with 1 mistmath; 0 reads with 2 mistmathes
***Note: Complete list of lost human reads is available from sam files: /u/home/s/serghei/code2/rop/example/test12/lostHumanReads/unmappedExample_genome.sam,/u/home/s/serghei/code2/rop/example/test12/lostHumanReads/unmappedExample_transcriptome.sam
3. Maping to repeat sequences...
-Identify 1 lost repeat sequences from unmapped reads
***Note : Repeat sequences classification into classes (e.g. LINE) and families (e.g. Alu) will be available in next release
4. Non-co-linear RNA profiling
***Note : Trans-spicing and gene fusions  are currently not supported, but will be in the next release.
--identified 2 reads from circRNA
***Note: circRNAs detected by CIRI are available here: unmappedExample_circRNA.csv
5a. B lymphocytes profiling...
--identified 1 reads mapped to immunoglobulin heavy (IGH) locus
--identified 0 reads mapped to immunoglobulin kappa (IGK) locus 
--identified 1 reads mapped to immunoglobulin lambda (IGL) locus
5b. T lymphocytes profiling...
--identified 0 reads mapped to T cell receptor alpha (TCRA) locus
--identified 2 reads mapped to T cell receptor beta (TCRB) locus
--identified 1 reads mapped to T cell receptor delta (TCRD) locus
--identified 0 reads mapped to T cell receptor gamma locus (TCRG) locus
In toto : 5 reads mapped to antibody repertoire loci
***Note : Combinatorial diversity of the antibody repertoire (recombinations of the of VJ gene segments)  will be available in the next release.
6.  Microbiome profiling...
--identified 0 reads mapped bacterial genomes
--identified 0 reads mapped viral genomes
--identified 4 reads mapped ameoba genomes
--identified 1 reads mapped crypto genomes
--identified 0 reads mapped giardia genomes
--identified 0 reads mapped microsporidia genomes
--identified 0 reads mapped piroplasma genomes
--identified 1 reads mapped plasmo genomes
--identified 1 reads mapped toxo genomes
--identified 0 reads mapped trich genomes
--identified 0 reads mapped tritryp genomes
In toto : 7 reads mapped to microbial genomes
Summary:   The ROP protocol is able to account for 2238 reads
***Unaccounted reads (not explained by ROP) are saved to /u/home/s/serghei/code2/rop/example/test12/unmappedExample_unaccountedReads.fasta
```

The ropOut directory now contains the output of ROP. The structure of the output is explained [here](https://github.com/smangul1/rop/wiki/ROP-output-details). For example it contains `/antibodyProfile/` directory with reads spanning antigen receptor gene rearrangement in the variable domain being identified by [IgBLAST](http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/igblast/release/1.4.0/). 

##Genomic profile of RNA-Seq

To get the genomic profile of the mapped reads use `gprofile.py`. To get the genomic profile of the toy bam file (reads from chr22) use the following command:

```
python gprofile.py example/mappedReads_chr22.bam /example/mappedReads_chr22.csv
```

The output of the module is number of reads assigned to each genomic category saved into the `/example/mappedReads_chr22.csv`

```
sampleName,nTotalMapped,nJunction,nCDS,nUTR3,nUTR5,nUTR_,nIntron,nIntergenic,nDeep,nMT,nMultiMapped
mappedReads,397134,129580,101541,96210,7457,22473,19420,3084,649,0,16720
```

You can use `/example/mappedReads_chr22.csv` to create pie chart. The  pie chart corresponding to `/example/mappedReads_chr22.csv` is presented bellow:

![Genomic profile of toy .bam file](https://sergheimangul.files.wordpress.com/2016/05/gprofile.png?w=1280)

Read more about Genomic Profile of RNA-Seq [here](https://github.com/smangul1/rop/wiki/ROP-output-details).


##Profile of repeat elements
```
python rprofile.py example/mappedReads_chr22.bam example/mappedReads_chr22_repeatProfile
```

The ROP provided three levels of repeat profile, i.e class level (e.g SINE, LINE); family level (e.g. Alu, L1); gene level (e.g. L1P4c).  The files contain relative proportions of repeat categories based on the number of reads from the category. More details about the repeat classes used by ROP are [here](https://github.com/smangul1/rop/wiki/What-is-ROP%3F).

Those the output files created for each level



Class level of classification (`mappedReads_chr22_repeatProfilerepeatClass.csv`): 


```
sample,LINE?,LTR,Satellite,Retroposon,DNA,SINE?,RNA,DNA?,RC,LINE,SINE,LTR?
mappedReads_chr22,0,2445,3,6,568,0,0,0,0,2588,5356,0
```

The classes with the `?` are provided by [RepeatMasker](http://www.repeatmasker.org/). For example [MER129](http://www.repeatmasker.org/cgi-bin/ViewRepeat?id=MER129) is classified as `LTR?`. You may merge LTR? with LTR or ignore those.

Those are the repeat classes with sufficient number of reads supporting the class. The repeat profile is presented as a table:
```
sample	LTR	DNA	LINE	SINE
mappedReads_chr22	2445	568	2588	5356
```

Alternatively you can visualize the repeat profile as a pie chart  

![](https://sergheimangul.files.wordpress.com/2016/05/rprofile_class2.png)

Family level of classification (`mappedReads_chr22_repeatProfilerepeatFamily.csv`) is presented bellow. Each family is represented in the following format `Class_Family` allowing to retrieve the particular class the family belongs to. For example `LINE____L1` corresponds to family element L1 from LINE class. 


```
sample,SINE____Alu,DNA____TcMar-Tigger,LTR?____LTR,DNA____PiggyBac,DNA____hAT,Retroposon____SVA_E,LTR____ERV1,LINE____L1,LINE____L2,DNA____MuDR,DNA____TcMar,LINE?____Penelope,LTR____ERVL-MaLR,DNA____DNA,LTR____ERV,DNA____hAT-Charlie,RC____Helitron,DNA____hAT-Blackjack,SINE____MIR,Satellite____centr,DNA____Merlin,LTR____LTR,DNA____hAT-Tip100,RNA____RNA,SINE____Deu,Retroposon____SVA_D,LTR____Gypsy,Retroposon____SVA_F,Retroposon____SVA_A,LINE____CR1,Retroposon____SVA_C,Retroposon____SVA_B,DNA____TcMar-Mariner,LINE____RTE,LINE____RTE-BovB,DNA____TcMar-Tc2,LTR____ERVK,Satellite____acro,Satellite____telo,LTR____ERVL,Satellite____Satellite,SINE?____SINE,DNA?____DNA,SINE____SINE,LINE____Dong-R4
mappedReads_chr22,4502,123,0,1,16,0,859,1476,1009,0,2,0,943,2,0,349,0,4,849,2,0,0,50,0,1,2,3,0,3,93,0,1,21,10,0,0,151,0,0,489,1,0,0,4,0
```

Those are the repeat families with sufficient number of reads supporting the family. The repeat profile is presented as a table:

```
sample	LTR____ERV1	LTR____ERVL-MaLR	LTR____ERVK	DNA____hAT-Charlie	DNA____TcMar-Tigger	DNA____hAT-Tip100	LINE____L1	LINE____L2	SINE____Alu	SINE____MIR
mappedReads_chr22	859	943	151	349	123	50	1476	1009	4502	849
```


Alternatively you can visualize the repeat profile as a pie chart  

![](https://sergheimangul.files.wordpress.com/2016/05/rprofile_family.png)


ROP also provides gene level resolution (`mappedReads_chr22_repeatProfilerepeatGene.csv`) , where the number of reads assigned to each repeat instance is reported. Each instance is represented in the following format `Family__Class__Instance` allowing to retrieve the particular class and family the repeat instance belongs to. For example `L1____LINE____L1MD` corresponds to element L1MD from L1 gamily of LINE class. 


