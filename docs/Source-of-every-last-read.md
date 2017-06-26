ROP is able to detect the source of every last RNA-Seq read. The reads not mapped to the refference genome might originate from complex RNA molecules, recombinant antibodies and microbial communities. The ROP accounts for 98.8% of all reads across poly(A) and ribo-depletion libraries. 

Please refer to the `.log` in the output directory (second mandatory command line argument). Logfile contains the information about the number of reads detected in each step of the ROP pipeline.  

An example of the `.log`

```
*********************************************
ROP is a computational protocol aimed to discover the source of all reads, originated from complex RNA molecules, recombinant antibodies and microbial communities. Written by Serghei Mangul (smangul@ucla.edu) and Harry Taegyun Yang (harry2416@gmail.com), University of California, Los Angeles (UCLA). (c) 2016. Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see:
http://serghei.bioinformatics.ucla.edu/rop/
*********************************************
Processing 2505 unmapped reads
1. Quality Control...
--filtered 2193 low quality reads
--filtered 2 low complexity reads (e.g. ACACACAC...)
--filtered 22 rRNA reads
In toto : 2217 reads failed QC and are filtered out
2. Remaping to human references...
--identified 6 lost human reads from unmapped reads 
3. Mapping to repeat sequences...
-Identify 1 lost repeat sequences from unmapped reads
***Note : Repeat sequences classification into classes (e.g. LINE) and families (e.g. Alu) will be available in next release
3. Non-co-linear RNA profiling
Please use --circRNA options to profile circular RNAs
***Note : Trans-spicing and gene fusions  are currently not supported, but will be in the next release.
4a. B lymphocytes profiling...
--identified 0 reads mapped to immunoglobulin heavy (IGH) locus
--identified 0 reads mapped to immunoglobulin kappa (IGK) locus 
--identified 0 reads mapped to immunoglobulin lambda (IGL) locus
4b. T lymphocytes profiling...
--identified 0 reads mapped to T cell receptor alpha (TCRA) locus
--identified 0 reads mapped to T cell receptor beta (TCRB) locus
--identified 0 reads mapped to T cell receptor delta (TCRD) locus
--identified 0 reads mapped to T cell receptor gamma locus (TCRG) locus
In toto : 0 reads mapped to antibody repertoire loci
***Note : Combinatorial diversity of the antibody repertoire (recombinations of the of VJ gene segments)  will be available in the next release.
5.  Microbiome profiling...
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
Summary:   The ROP protocol is able to account for 2231 reads
```

Additionally ROP creates ``numberReads_<sampleName>.log`` with the number of reads detected in each step in the `.csv` format. An example of ``numberReads_<sampleName>.log`` is provided bellow:

```
sample,totalReads,lowQuality,lowComplexity,rRNA,lostHumanReads,lostRepeatReads,immuneReads,microbiomeReads,unaccountedReads
unmappedExample,2508,2193,2,22,6,1,5,7,272
```

ROP is also able to find the the source of mapped reads. ROP is able to categorize the mapped reads into genomic and repeat features. An example of a file with reads assigned to the genomic features is presented bellow:

```
SRR1146076.3540543,22,junction
SRR1146076.17070789,22,UTR5
SRR1146076.21758998,22,UTR5
SRR1146076.20543213,22,UTR5
SRR1146076.9808865,22,DEEP
SRR1146076.16863731,22,INTERGENIC
SRR1146076.968340,22,UTR3
```


Given ROP analysis for multiple sample, you may consider concatenate the individual logfiles to get the statistics across the multiple samples. Detail instructions are provided [here](https://github.com/smangul1/rop/wiki/ROP-analysis-of-multiple-samples)