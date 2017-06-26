ROP uses the following reference databases:

- [Ensembl GRCh37](http://www.ensembl.org/Homo_sapiens/Info/Index)
- [Human ribosomal DNA complete repeating unit](ftp://ftp.ncbi.nih.gov/)
- [RepBase20.07](http://www.girinst.org/repbase/)
- [IMGT V(D)J genes of  B and T cell receptors](http://www.imgt.org/vquest/refseqh.html#V-D-J-C-sets)
- [NCBI bacterial reference genomes](https://ftp.ncbi.nlm.nih.gov/genomes/)
- [NCBI viral reference genomes](https://ftp.ncbi.nlm.nih.gov/genomes/Viruses)
- [EuPathDB  eukaryotic reference genomes](http://eupathdb.org/eupathdb/)
- [Repeat annotations](http://labshare.cshl.edu/shares/mhammelllab/www-data/TEToolkit/TE_GTF/)

# Repeats
We use GTF files generated from the RepeatMasker annotations by Jin, Ying, et al. and downloaded from http://labshare.cshl.edu/shares/mhammelllab/www-data/TEToolkit/TE_GTF/.  Repeat elements overlapping CDS regions are excluded from the analysis. We filtered out 6,873 repeat elements overlapping CDS regions.

Prepared repeat annotations (bed formatted file) are available at https://github.com/smangul1/rop/tree/master/source/rprofile/annotations/human/bedPreparedhg19_rmsk_TE_prepared_noCDS.bed.gz

The prepared repeat annotations contain 8 Classes and 43 Families.

# Third Party Software

ROP is distributed with several open source components that were developed by other groups. These components are (c) their respective developers and are redistributed with ROP to provide ease-of-use. Please see the list of tools (not exhaustive) for licensing details:

- [SEQLEAN](https://sourceforge.net/projects/seqclean/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [BWA](http://bio-bwa.sourceforge.net/)
- [CIRI](https://sourceforge.net/projects/ciri/)
- [Samtools](http://www.htslib.org/)

