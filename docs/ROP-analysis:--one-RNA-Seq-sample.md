In this tutorial we will show how run ROP for one RNA-Seq sample. We use RNA-Seq sample of normal skin (SRR1146076)  downloaded from [here](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54456). Please note 
that proposed sample is not necessarily the most typical RNA-Seq sample and is provided for demonstration purposes. 

In this tutorial we provide the mapped and unmapped reads of the RNA-Seq sample. The size of the unmapped reads in .fastq format is  1.4G. The size of the unmapped reads in .bam format is 0.3G. The size of the mapped reads (.bam) is 1.2G. 

In general case, you will need to map the reads with any of available high-throughput aligners (e.g. [tophat2](https://ccb.jhu.edu/software/tophat/index.shtml), [STAR](https://github.com/alexdobin/STAR)) and save unmapped reads in .bam (binary file, requires less space) or .fastq (text format) format. The instructions how to map the reads and save the unmapped reads are provided [here](https://github.com/smangul1/rop/wiki/How-to-map-reads-and-save-unmapped-reads). 

Please make sure that the basic unix commands (wget, python, perl) are available on the cluster.  Please install ROP first. The instructions how to install ROP are provided [here](https://github.com/smangul1/rop/wiki/How-to-install-ROP%3F) 

The first operation consists in navigating to ROP directory and creating a subdirectory for storing the training data. 

```
cd rop
mkdir tutorial
cd tutorial
mkdir data
```

Now, download the unmapped reads from RNA-Seq (please see https://github.com/prasmussen/gdrive regarding the gdrive utility)

```
gdrive download 0B_NUyiE86yDwdFFwZWFTck9QQjg
tar -xvf skinExample.tar
```

Now, you are ready to analyze the RNA-Seq sample using ROP. We are running ROP using the default options. ROP is an intensive pipeline requiring substantial amount of computational resources. Thus we don't recommend to run ROP from login nodes. Please check the policy of you cluster, from where to run the ROP pipeline. For hoffman2 (UCLA cluster) read the policy [here] (http://ccn.ucla.edu/wiki/index.php/Hoffman2:Interactive_Sessions). 

ROP requires two mandatory command line arguments, i.e. (1) the unmapped reads and (2) the directory to save the results of ROP.

```
usage: python rop.py [-h] [--qsub] [--qsubArray] [--b] [--skipLowq] [--skipQC]
                     [--circRNA] [--immune] [--gzip] [--quiet] [--dev]
                     [--license]
                     unmappedReads dir
```

To run ROP for unmapped reads in .bam format, use --b option
```
python rop.py --b tutorial/data/unmapped_SR_1146076.bam /tutorial/ropOut/
```

To run ROP for unmapped reads in .fastq format 

```
python rop.py  tutorial/data/unmapped_SR_1146076.fastq /tutorial/ropOut/
```

You should expect the following output of the ROP pipeline on your screen: 

```
*********************************************
ROP is a computational protocol aimed to discover the source of all reads, originated from complex RNA molecules, recombinant antibodies and microbial communities. Written by Serghei Mangul (smangul@ucla.edu) and Harry Taegyun Yang (harry2416@gmail.com), University of California, Los Angeles (UCLA). (c) 2016. Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see:
https://sergheimangul.wordpress.com/rop/
https://github.com/smangul1/rop/wiki
*********************************************
Processing 7687003 unmapped reads
1. Quality Control...
--filtered 6966209 low quality reads
--filtered 12457 low complexity reads (e.g. ACACACAC...)
--filtered 55299 rRNA reads
In toto : 7033965 reads failed QC and are filtered out
2. Remaping to human references...
--identified 16278 lost human reads from unmapped reads 
3. Maping to repeat sequences...
-Identify 3862 lost repeat sequences from unmapped reads
***Note : Repeat sequences classification into classes (e.g. LINE) and families (e.g. Alu) will be available in next release
4. Non-co-linear RNA profiling
***Note : Trans-spicing and gene fusions  are currently not supported, but will be in the next release.
--identified 1924 reads from circRNA
5a. B lymphocytes profiling...
--identified 720 reads mapped to immunoglobulin heavy (IGH) locus
--identified 672 reads mapped to immunoglobulin kappa (IGK) locus 
--identified 403 reads mapped to immunoglobulin lambda (IGL) locus
5b. T lymphocytes profiling...
--identified 33 reads mapped to T cell receptor alpha (TCRA) locus
--identified 55 reads mapped to T cell receptor beta (TCRB) locus
--identified 10 reads mapped to T cell receptor delta (TCRD) locus
--identified 8 reads mapped to T cell receptor gamma locus (TCRG) locus
In toto : 1901 reads mapped to antibody repertoire loci
***Note : Combinatorial diversity of the antibody repertoire (recombinations of the of VJ gene segments)  will be available in the next release.
6.  Microbiome profiling...
--identified 1577 reads mapped bacterial genomes
--identified 34 reads mapped viral genomes
--identified 5683 reads mapped ameoba genomes
--identified 3157 reads mapped crypto genomes
--identified 20 reads mapped giardia genomes
--identified 133 reads mapped microsporidia genomes
--identified 359 reads mapped piroplasma genomes
--identified 2961 reads mapped plasmo genomes
--identified 8076 reads mapped toxo genomes
--identified 1 reads mapped trich genomes
--identified 133 reads mapped tritryp genomes
In toto : 22134 reads mapped to microbial genomes
Summary:   The ROP protocol is able to account for 7078140 reads
***Unaccounted reads (not explained by ROP) are saved to /u/scratch/b/brigitta/skin/out/unmapped_SR_1146076_unaccountedReads.fasta
```

Alternatively, you can access the ROP message from the `unmapped_SR_1146076.log` logfile. 
 
The `/ropOut/` directory now contains the output of ROP. To navigate to the ropOut directory use this command 

```
cd /tutorial/ropOut/
```

The directory contains individual directions for each types of the ROP analysis. The ropOut directory now contains the output of ROP. The structure of the output is explained here. For example it contains /antibodyProfile/ directory with reads spanning antigen receptor gene rearrangement in the variable domain being identified by IgBLAST. The structure of the ROP output is explained [here](https://github.com/smangul1/rop/wiki/ROP-output-details)


###Genomic profile of RNA-Seq

To get the genomic profile of the mapped reads use `gprofile.py`. It requires pysam to be installed. On the hoffman2 you may use `module` command to load pysam as follows:


```
module load python/2.7.3
```


First navigate to '/tutorial/data/' a subdirectory to store the training data. 

```
cd /tutorial/data/
```

Then download the mapped reads using the following command:


```
wget https://googledrive.com/host/0B_NUyiE86yDwaUxoVjhlSjN5SkE/skinExample.tar
tar -xvf skinExample.tar
```

This is how to run `gprofile.py`:

```
python gprofile.py /tutorial/data//mapped_SR_1146076.bam /tutorial/data/mapped_SR_1146076_genomicProfile.csv
```

The output of the module is number of reads assigned to each genomic category saved into the `/tutorial/data/mapped_SR_1146076_genomicProfile.csv`


```
sampleName,nTotalMapped,nJunction,nCDS,nUTR3,nUTR5,nUTR_,nIntron,nIntergenic,nDeep,nMT,nMultiMapped
mapped_SR_1146076,17904083,5425835,4535762,3589107,362853,963938,765359,195246,32345,1061075,972563
```

You can use `/tutorial/data/mapped_SR_1146076_genomicProfile.csv` to create pie chart. The pie chart corresponding to ``/tutorial/data/mapped_SR_1146076_genomicProfile.csv`` is presented bellow:

![](https://sergheimangul.files.wordpress.com/2016/05/gprofile1.png?w=1280)







![](https://sergheimangul.files.wordpress.com/2016/05/rprofile_class4.png)


![](https://sergheimangul.files.wordpress.com/2016/05/rprofile_family1.png)


More details about additional options and strategies of the ROP are available [here](https://github.com/smangul1/rop/wiki/Additional-options)