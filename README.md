# ROP 1.0 : Read origin protocol
Authors: Serghei Mangul (smangul@ucla.edu), Harry Yang (Taegyun) (harry2416@gmail.com)

##Description
ROP is a computational protocol for profiling the composition of unmapped reads, which failed to map to the human references. 
 
 
 
ROP prococol consist of six steps to characterize the unmapped reads:

1. Quality control. Exclude low-quality, low-complexity and rRNA reads (reads mathing rRNA repeat unit) 
2. Identify lost human reads, which are missed due to the heuristics implemented for computational speed in conventional aligners. These include reads with mismatches and short gaps relative to the reference set, but can also include perfectly matched reads.  
3. Identify lost repeat reads, by mapping unmapped reads to the database of repeat sequences (Repbase 20.7)
4. Identify ‘non-co-linear’ RNAs reads from circRNAs, gene fusions, and trans-splicing events, which combine sequence from distant elements.
5. Identify reads from recomobinations of B and T cell receptors i.e. V(D)J recombinations
6. Identify mircobial reads using microbial genomes and phylogenetic marker genes to identify microbial reads and assign them to corresponding taxa

If you use this software, please cite :

##Pre-requisites

##Installation

