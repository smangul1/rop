## antibodyProfile directory

### Note: This is relevant to versions of ROP using IgBLAST (versions <1.0.6). Starting with ROP v1.0.6, we we have switched from IgBLAST to ImReP to profile B and T cell receptor repertoires.  ImReP shows superior accuracy compared to existing tools  (see our manuscript “Profiling adaptive immune repertoires across multiple human tissues by RNA Sequencing” available at bioRxiv).

This directory contains the output of [Step 5. B and T lymphocytes profiling](https://github.com/smangul1/rop/wiki/What-is-ROP%3F). 

Starting from release ROP v1.0.6 we are using accurate in-house tool to detect reads spanning V(D)J recombination to quanity individual immune response. 

It contains a separate directory for B cell receptors (BCR) : `BCR` directory and a separate directory for T cell receptor (TCR) :  `TCR` directory. 

BCR directory contains:

* IGH a directory for immunoglobulin heavy locus profile 
* IGK a directory for immunoglobulin kappa locus
* IGL a directory for immunoglobulin lambda locus profile 

TCR directory contains:

* TCRA a directory for T cell receptor alpha locus
* TCRB a directory for T cell receptor beta locus
* TCRG a directory for T cell receptor gamma locus
* TCRD a directory for T cell receptor delta locus

Reads spanning antigen receptor gene rearrangement in the variable domain are identified using [IgBLAST](http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/igblast/release/1.4.0/).  IgBLAST reports
alignment of the reads to the variable (V) gene, the diversity (D) gene and the joining (J) gene, or the recombination of those. Reads alignment is saved into the `_igblast.csv` (modified tabular output format 6).

n | id| What does it mean? 
:-- | :-- | :--
0 | VDJ | gene segment 
1 | qseqid |read name  
2 | sseqid |reference genome   
3 | pident |percentage of identical matches  
4 | length |alignment length
5 | mismatch | number of mismatches
6|	 gapopen	| number of gap openings
 7|	 qstart	| start of alignment in query
 8|	 qend	| end of alignment in query
 9|	 sstart	| start of alignment in subject
 10|	 send	| end of alignment in subject
 11|	 evalue	| expect value
 12|	 bitscore	| bit score

An example of the unmapped read aligned to the immunoglobulin heavy variable 1-46 gene (IGHV1-46)

```
V       SRR1146076.56325        IGHV1-46*01     93.67   79      5       0       0       1       79      203     281     7e-27     109
```
