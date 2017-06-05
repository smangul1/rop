# ROP: Read Origin Protocol

```
Don’t let your unmapped reads go to waste
```

ROP is a computational protocol aimed to discover the source of all reads, 
which originate from complex RNA molecules, recombinant B and T cell receptors and 
microbial communities. The ROP accounts for 98.8% of all reads across poly(A) 
and ribo-depletion protocols, compared to 83.8% by conventional reference-based 
protocols. ROP profiles repeats, circRNAs, gene fusions, trans-splicing events, 
recombined B/T-cell receptor sequences and microbial communities. The 'dumpster 
diving' profile of unmapped reads output by our method is not limited to 
RNA-Seq technology and may be applied to whole-exome and whole-genome 
sequencing.

In general case, you need to map the reads with any of available high-throughput aligners (e.g. STAR, tophat2) and save unmapped reads in .bam (binary format, requires less space) or .fastq (text format). The instructions how to map the reads and save the unmapped reads are provided in the ROP tutorial. The comphehensive benchmarking of RNA alligners is provided in the recent Nature Methods study : Baruzzo, Giacomo, et al. "Simulation-based comprehensive benchmarking of RNA-seq aligners." Nature methods 14.2 (2017): 135-139.

# ROP Tutorial

ROP Tutorial is available here : https://github.com/smangul1/rop/wiki

# ROP protocol consists of two (optional) modules to categorize the mapped reads:

- We developed gprofiler, a tool to categorize the mapped reads into genomic categories (CDS, UTR, intons, etc).
- We developed rprofiler, a tool to profile repetitive elements (e.g. SINEs, LINEs, LTRs).

# ROP protocol consists of six steps to characterize the unmapped reads:

- Quality control. Exclude low-quality, low-complexity and rRNA reads (in-house code, SEQLEAN, Megablast )
- Identify lost human reads, which are missed (not mapped) due to the heuristics implemented for computational speed in conventional aligners. These include reads with mismatches and short gaps relative to the reference set, but can also include perfectly matched reads (Bowtie2, previously Megablast)
- Identify lost repeat sequences, by mapping unmapped reads onto the database of repeat sequences (BWA, previosly Megablast )
- Identify ‘non-co-linear’ RNAs reads from circRNAs, gene fusions, and trans-splicing events, which combine sequences from distant elements (CIRI)
- Identify reads from recombinations of B and T cell receptors i.e. V(D)J recombinations (ImReP, previously IgBLAST)
- Profile microbial communities using the microbial reads mapped onto the microbial genomes and marker genes (Megablast, MetaPhlAn)


## Databases

- Ensembl GRCh37/hg19
- Human ribosomal DNA complete repeating unit 
- Repeat elements (RepBase20.07)
- V(D)J genes of  B and T cell receptors
- Bacterial genomes
- Viral genomes
- Genomes of eukaryotic pathogens
- Repeat annotations

## Third Party software 

ROP is distributed with several open source components that were developed by other groups. These components are (c) their respective developers and are redistributed with ROP to provide ease-of-use. Please see the list of tools (not exhaustive) for licensing details:

- SEQLEAN
- BLAST+ (Megablast option)
- Bowtie2
- CIRI
- MetaPhlAn
- Samtools
- Bamtools
- HTSeq



We switched to SourceForge to store the official releases of ROP: 
https://sourceforge.net/projects/rop2/. The code at github is in permanent 
development and should be used for development only. Consider contacting 
smangul@ucla.edu.

For more details see: https://sergheimangul.wordpress.com/rop/

