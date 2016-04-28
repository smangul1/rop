# ROP : Read origin protocol
Authors: Serghei Mangul (smangul@ucla.edu), Harry Yang (Taegyun) (harry2416@gmail.com)

##Description
ROP is a computational protocol for profiling the composition of unmapped reads, which failed to map to the human references. ROP profiles repeats, circRNAs, gene fusions, trans-splicing events, recombined B and T cell receptors and microbial communities. ROP leverages accurate alignment methods for both host and microbial sequences and able to account for 98.6%  of the reads compared to 83.8% by conventional reference-based protocols. The ‘dumpster diving’ profile of unmapped reads output by our method is not limited to RNA-Seq technology and might be applied further to whole-exome and whole-genome sequencing. 
 
ROP prococol consist of six steps to characterize the unmapped reads:

1. Quality control. Exclude low-quality, low-complexity and rRNA reads (reads mathing rRNA repeat unit). 
2. Identify lost human reads, which are missed due to the heuristics implemented for computational speed in conventional aligners. These include reads with mismatches and short gaps relative to the reference set, but can also include perfectly matched reads.  
3. Identify lost repeat reads, by mapping unmapped reads to the database of repeat sequences (Repbase 20.7).
4. Identify ‘non-co-linear’ RNAs reads from circRNAs, gene fusions, and trans-splicing events, which combine sequence from distant elements.
5. Identify reads from recomobinations of B and T cell receptors i.e. V(D)J recombinations.
6. Profile taxonomic composition of microbial communities using the microbial reads mapped onto the  microbial genomes and marker genes.

ROP prococol consists of two steps to characterize mapped reads:

1. We developed gProfile. This tools categorize human reads into genomic categories (CDS, UTR, intons, etc).   
2. We developed rProfile, a tool to profile repetitive elements (e.g. SINEs, LINEs, LTRs).

If you use this software, please cite :

##Pre-requisites

##Installation

For standard installation: 
```bash 
python install.py --standard
```

For immune (i.e. TCR, BCR) analysis only: 
```bash
python install.py --immune
``` 


##Training Dataset

Use RNA-Seq sample of normal skin (SRR1146076) as a training sample. RNA-Seq sample was downloaded from here:
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54456


##How to run ROP

##3rd party software

ROP is distributed with several open source components that were developed by other groups.
These components are (c) their respective developers and are redistributed with PhyloSift to provide ease-of-use.

Please see the following web sites for licensing details and source code for these other components:


```bash
test
```

##Extract VDJ combinations from IgBLAST for miltiple samples

```bash
for f in data/*csv; do python <pathROP>/rop/src/iprofile/extract_VDJ.py ${f} $PWD/VJ_e20/ IGK 1e-20;done
```



##Contact

For issues with this software, contact smangul@ucla.edu. Otherwise, submit an issue directly through Github.

##License

This project is released under the GPL-3 License. Please view the LICENSE file for more details.
