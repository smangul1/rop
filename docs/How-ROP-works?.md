ROP protocol consists of two (optional) modules to categorize the mapped reads:
* Genomic profile of RNA-Seq. We developed [gprofile](https://github.com/smangul1/gprofile), a  tool to categorize mapped reads into genomic categories (CDS, UTR, intons, etc)  (details)
* Profile of repeat elements. We developed [rprofile](https://github.com/smangul1/rprofile), a tool to profile repetitive elements (e.g. SINEs, LINEs, LTRs)

ROP protocol consists of six steps to characterize the unmapped reads:

1. Quality control. Exclude low-quality, low-complexity and rRNA reads ([FASTX](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html), [SEQCLEAN](https://sourceforge.net/projects/seqclean/), [Megablast](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/))
2. Identify lost human reads, which are missed due to the heuristics implemented for computational speed in conventional aligners. These include reads with mismatches and short gaps relative to the reference set, but can also include perfectly matched reads ([Megablast](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/))
1. Identify lost repeat sequences, by mapping unmapped reads onto the database of repeat sequences ([Megablast](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) )
1. Identify ‘non-co-linear’ RNAs reads from circRNAs, gene fusions, and trans-splicing events, which combine sequence from distant elements (ncSplice, [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) , [CIRI](https://github.com/Frenzchen/ncSplice))
1. Identify reads from recombinations of B and T cell receptors i.e. V(D)J recombinations ([IgBLAST](http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/igblast/release/1.4.0/))
1. Profile taxonomic composition of microbial communities using the microbial reads mapped onto the microbial genomes and marker genes ([Megablast](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), [MetaPhlAn](http://huttenhower.sph.harvard.edu/metaphlan))


##Genomic profile of RNA-Seq

We developed [gprofile](https://github.com/smangul1/gprofile) to categorize the mapped reads into genomic categories based on the compatibility of each read with the features defined by gene annotations. Genomic profile can be used to benchmark different sequencing platform and library preparation methods, as well as assess the efficiency of rRNA depletion and level of the sample degradation.  

Those are the categories of the genomic profile:  

* read mapped to multiple locations on the reference genome is categorized as a multi-mapped read
* read fully contained within the CDS, intron, UTR3, or UTR5 boundaries of a least one transcript is classified as a CDS, intronic, UTR3, or UTR5, respectively
* read simultaneously overlapping UTR3 and UTR5 regions is classified as a UTR read
* read spanning exon-exon boundary is defined as a junction read
* read mapped outside of gene boundaries and within a proximity of 1Kb is defined as a (proximal) inter-genic read
* read mapped outside of gene boundaries and beyond the proximity of 1Kb is defined as a deep inter-genic read
* read mapped to mitochondrial DNA is classified as a mitochondrial read
* reads from a pair mapped to different chromosomes are classified as a fusion reads

Given the assignment of reads into the genomic categories, ROP will estimate relative proportions of the categories  based on the number of reads from the category. 


##Profile of repeat elements 

We developed [rprofile](https://github.com/smangul1/rprofile) to categorize the mapped reads into repeat categories based on the compatibility of each read with repeat instances defined by [RepeatMasker](http://www.repeatmasker.org/) [[more details] (http://serghei.bioinformatics.ucla.edu/rop/repeats/)]

Reads are categorized into the following classes: 

* LINE (Long Interspersed Nuclear Elements)
* SINE (Short Interspersed Nuclear Elements)
* LTR (Long terminal repeat)
* RC
* SVA
* RNA
* Satellite
* Retroposon


Reads are categorized into the following families: 
* acro
* Alu
* centr
* CR1
* Deu
* DNA
* Dong-R4
* ERV
* ERV1
* ERVK
* ERVL
* ERVL-MaLR
* Gypsy
* hAT
* hAT-Blackjack
* hAT-Charlie
* hAT-Tip100
* Helitron
* L1
* L2
* LTR
* Merlin
* MIR
* MuDR
* Penelope
* PiggyBac
* RNA
* RTE
* RTE-BovB
* Satellite
* SINE
* SVA_A
* SVA_B
* SVA_C
* SVA_D
* SVA_E
* SVA_F
* TcMar
* TcMar-Mariner
* TcMar-Tc2
* TcMar-Tigger
* telo

Read are categorized into the individual repeat instances (e.g. L1P4c). 

Given the assignment of reads into the individual repeat categories, ROP will estimate relative proportions of the categories  based on the number of reads from the category. 


## Step 1. Quality control

Filter out low quality, low complexity (e.g. ACACACAC...), and rRNA reads

* low quality reads are identified by [FASTX](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html). Low quality reads are defined as reads that have quality lower than 30 in at least 75% of their base pairs
* Low complexity are identified by [SEQCLEAN](https://sourceforge.net/projects/seqclean/)
* rRNA reads are identified as reads mapped to the [rRNA repeat sequence](http://www.ncbi.nlm.nih.gov/nuccore/U13369) by [Megablast](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

## Step 2. Remap to human sequences

Identify lost human reads, which are missed due to the heuristics implemented for computational speed in conventional aligners. These include reads with mismatches and short gaps relative to the reference set, but can also include perfectly matched reads 

## Step 3. Map to repeat sequences
Identify lost repeat sequences, by mapping unmapped reads onto the [database](http://www.girinst.org/repbase/) of repeat sequences using [Megablast](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 

## Step 4. Non-co-linear (NCL) RNA profiling

Identify ‘non-co-linear’ RNAs reads from circular RNAs, gene fusions, and trans-splicing events, which combine sequence from distant elements ([ncSplice](https://github.com/Frenzchen/ncSplice))

We identify:

* read spliced distantly on the same chromosome supports trans-splicing event
* read spliced across different chromosomes supports gene fusion event
* reads spliced in a head-to-tail configuration supports circRNAs

##Step 5. B and T lymphocytes profiling

Mapped and unmapped reads are used to survey the human antibody repertoire. Reads entirely aligned to B cell receptors (BCR) and T cell receptors (TCR)  genes are extracted from the mapped reads (.bam). Reads with extensive somatic hyper mutations (SHM) and reads arising from V(D)J recombination are identified by  [IgBLAST](http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/igblast/release/1.4.0/) from the unmapped reads (.fastq or .bam). 

Thus reads mapped to the constant regions of the antigen receptors allows to estimate the relative proportions of major antibodies classes. For example, based on mapped reads ROP estimates the relative proportions of the main antibody isotypes (IgA, IgD, IgE, IgG, and IgM).


##Step 6. Microbiome profiling

The remaining unmapped reads are potentially microbial in origin.  We use microbial genomes ([Megablast](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)) and phylogenetic marker genes ([MetaPhlAn](http://huttenhower.sph.harvard.edu/metaphlan)) to identify microbial reads and assign them to corresponding taxa. Microbial reads could have been introduced by contamination or the natural microbiome of the samples, which includes viral, bacterial, or other microbial species. 



