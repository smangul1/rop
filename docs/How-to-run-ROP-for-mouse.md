Starting with release 1.0.5, ROP can accept unmapped reads from mouse data. Thanks to Kevin Hsieh(kevin.hsieh@ucla.edu), Linus Chen (u6.30cl@gmail.com) for developing ROP-mouse. 

Make sure to download the mouse database. Please follow [How to install ROP?](https://github.com/smangul1/rop/wiki/How-to-install-ROP%3F)

To run ROP for mouse use `--organism mouse` option. You can use `example/mouse/unmappedExampleMouse.fastqa` to test the ROP-mouse. Toy example distributed with the ROP.  

```
python rop.py --organism mouse example/mouse/unmappedExampleMouse.fastq testMouse/
```

This is expected output of ROP-mouse:

```
Processing 3044 unmapped reads of length 100
1. Quality Control...
--filtered 1163 low quality reads
--filtered 207 low complexity reads (e.g. ACACACAC...)
--filtered 135 rRNA reads
In toto : 1505 reads failed QC and are filtered out
2. Remapping to references...
--identified 59 lost reads from unmapped reads. Among those: 43 reads with 0 mismatches; 13 reads with 1 mismatch; 3 reads with 2 mismatches
***Note: Complete list of lost reads is available from sam files: /u/scratch/b/brigitta/rop/testMouse/lostReads/unmappedExampleMouse_genome.sam,/u/scratch/b/brigitta/rop/testMouse/lostReads/unmappedExampleMouse_transcriptome.sam
3. Mapping to repeat sequences...
-- Identified 61 lost repeat sequences from unmapped reads
***Note : Repeat sequences classification into classes (e.g. LINE) and families (e.g. Alu) will be available in next release
4. Non-co-linear RNA profiling
***Note : Trans-spicing and gene fusions  are currently not supported, but will be in the next release.
--identified 9 reads from circRNA
***Note: circRNAs detected by CIRI are available here: unmappedExampleMouse_circRNA.tsv
5a. B lymphocytes profiling
--identified 7 reads mapped to immunoglobulin heavy (IGH) locus
--identified 3 reads mapped to immunoglobulin kappa (IGK) locus 
--identified 2 reads mapped to immunoglobulin lambda (IGL) locus
5b. T lymphocytes profiling...
--identified 3 reads mapped to T cell receptor alpha (TCRA) locus
--identified 8 reads mapped to T cell receptor beta (TCRB) locus
--identified 4 reads mapped to T cell receptor delta (TCRD) locus
--identified 4 reads mapped to T cell receptor gamma locus (TCRG) locus
In toto : 31 reads mapped to antibody repertoire loci
***Note : Combinatorial diversity of the antibody repertoire (recombinations of the of VJ gene segments)  will be available in the next release.
***Extra step.  Metaphlan profiling...
***Microbiome profiling by Metaphlan2: taxonomic profile of microbial communities detected by Metaphlan2 is available here: /u/scratch/b/brigitta/rop/testMouse/microbiomeProfile//metaphlan/
6.  Microbiome profiling...
--identified 189 reads mapped bacterial genomes
--identified 2 reads mapped viral genomes
--identified 3 reads mapped ameoba genomes
--identified 1 reads mapped crypto genomes
--identified 1 reads mapped giardia genomes
--identified 2 reads mapped microsporidia genomes
--identified 1 reads mapped piroplasma genomes
--identified 2 reads mapped plasmo genomes
--identified 2 reads mapped toxo genomes
--identified 2 reads mapped trich genomes
--identified 1 reads mapped tritryp genomes
In toto : 206 reads mapped to microbial genomes
Summary: The ROP protocol is able to account for 1871 reads
***Unaccounted reads (not explained by ROP) are saved to /u/scratch/b/brigitta/rop/testMouse/unmappedExampleMouse_unaccountedReads.fasta
***Log file with all the commands used is available here: /u/scratch/b/brigitta/rop/testMouse/dev.log
********************
Important: ROP relies on  several open source tools that were developed by other groups. These components are (c) their respective developers and are redistributed with ROP to provide ease-of-use. The list of the tools used by ROP and the parameters/reference databases are provided here: /u/scratch/b/brigitta/rop/testMouse/tools.log 
```


