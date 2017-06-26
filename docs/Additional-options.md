In case you are interested to skip STEP1 (Quality Control) please use option `--skipQC`. Please note that in this case low quality , low complexity and rRNA reads will not be filtered out. The input reads must be in the FASTA format (`.fa`  or `.fasta` extension)

Please run ROP as follows:

```
python rop.py --skipQC example/test.fa example/ropOut60/
```


In case you are interested in unmapped reads with the reads from QC step to be filtered out, please use `--dev` option to save the `afterQC.fastq`.  

##Targeted analysis  

Functionality to run the analysis of interest is supported starting from [ROP  1.0.1](https://sourceforge.net/projects/rop2/files/) (release 05/16/2016). Note that (step 1) QC and  (step 2) Remapping to human references (lost human reads) are mandatory. Please use the following option to run the analysis of interest:

* the --immune option now allows to run the immune profiling only
* the --microbiome option now allows to run the microbiome profiling only
* the --repeat option now allows to run the lost repeat profiling only

For example using `--repeat` option you can run lost repeat profiling only (STEP3) only, as follows:

```
python rop.py --repeat example/unmappedExample.fastq example/ropOut69/
```

The expected output is:

```
*********************************************
ROP (version 1.0.3) is a computational protocol aimed to discover the source of all reads, originated from complex RNA molecules, recombinant antibodies and microbial communities. Written by Serghei Mangul (smangul@ucla.edu) and Harry Taegyun Yang (harry2416@gmail.com), University of California, Los Angeles (UCLA). (c) 2016. Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see:
https://sergheimangul.wordpress.com/rop/
*********************************************
Processing 2511 unmapped reads of length 79
1. Quality Control...
--filtered 2193 low quality reads
--filtered 2 low complexity reads (e.g. ACACACAC...)
--filtered 22 rRNA reads
In toto : 2217 reads failed QC and are filtered out
2. Remapping to human references...
--identified 6 lost human reads from unmapped reads. Among those: 4 reads with 0 mismatches; 2 reads with 1 mismatch; 0 reads with 2 mismatches
***Note: Complete list of lost human reads is available from sam files: /u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut69/lostHumanReads/unmappedExample_genome.sam,/u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut69/lostHumanReads/unmappedExample_transcriptome.sam
*********************************
Non-substractive mode is selected : Low quality, low complexity, rRNA reads and lost human reads are filtered out. Resulting high quality non-human reads are provided as input  for STEP3-STEP6
*********************************
3. Maping to repeat sequences...
-- Identified 1 lost repeat sequences from unmapped reads
***Note : Repeat sequences classification into classes (e.g. LINE) and families (e.g. Alu) will be available in next release
4. Non-co-linear RNA profiling
***Note : Trans-spicing and gene fusions  are currently not supported, but will be in the next release.
4. Non-co-linear RNA profiling is skipped.
5a. B lymphocytes profiling is skipped.
5b. T lymphocytes profiling is skipped.
Extra step.  Metaphlan profiling is skipped.
6.  Microbiome profiling is skipped.
********************
Important: ROP relies on  several open source tools that were developed by other groups. These components are (c) their respective developers and are redistributed with ROP to provide ease-of-use. The list of the tools used by ROP and the parameters/reference databases are provided here: /u/home/galaxy/collaboratory/serghei/code/rop/example/ropOut69/tools.log 
```


