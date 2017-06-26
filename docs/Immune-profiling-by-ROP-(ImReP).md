Some mapping tools produce partially-mapped reads (i.e. STAR). In case read is mapped to BCR or TCR genes and is partially mapped to V or J gene, such read may contain CDR3 sequence and be used to assemble full-length CDR3 sequences. Because this read is not among the unmapped reads, ImReP will miss it. 

Starting from release 1.0.7 we offer rop-imrep.sh and rop-imrepGRCh38.sh, which will extract reads mapped to BCR and TCR loci and will merge with the unmapped reads to perform immune profiling using ImReP.

Given the bam file with mapped and unmapped reads you can run ROP-ImRep using this command:

```
~/code2/rop/rop-imrep.sh <bam> ROP ""
```

```
********************************************************************************
rop-imrep.sh is now available under ROP release v1.0.7
rop-imrep.sh is a script to run ROP-ImReP for bam file with mixture of mapped and unmapped reads.


ImReP is written by Igor Mandric and Serghei Mangul.

Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see: https://sergheimangul.wordpress.com/rop/
ROP Tutorial: https://github.com/smangul1/rop/wiki 
********************************************************************************
[1] bam with mapped and unmapped reads, BAm file needs to be indexed!!!
[2] dir to save results of ROP
[3] addditional options for ROP. Should be in double quaotes. In case you don't need addditional options for ROP use ""
********************************************************************************
```

For example: 

```
~/code2/rop/rop-imrep.sh G21033.GTEX-T6MO-0426-SM-32QOI.1.bam ROP5 ""
```