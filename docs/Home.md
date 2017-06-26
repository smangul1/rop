#  ROP : Read Origin Protocol 

![](https://github.com/smangul1/rop/blob/master/docs/rop.png)


## Quick Start 
ROP discovers the source of unmapped reads, which originated from complex RNA molecules, recombinant B and T cell receptors and microbial communities.
 
Download ROP using 
```
github clone https://github.com/smangul1/rop.git
```
Install ROP from the base directory

```
cd rop
./install.sh
```

Download reference database 
```
python getDB.py ~/
```

Run ROP analysis by a single command

```
python rop.py example/unmappedExample.fastq example/ropOut/
```

Find ROP analysis in `example/ropOut/` directory. Learn more about ROP output  [here](https://github.com/smangul1/rop/wiki/ROP-output-details) 

# ROP tutorial

Use the sidebar on the right to navigate ROP tutorial. Get started with a toy example of 2000 reads distributed with ROP package. 

# What is ROP

ROP is a computational protocol aimed to discover the source of all unmapped, which originate from complex RNA molecules, recombinant B and T cell receptors and microbial communities. We have tested ROP on 1 trillion reads from 10641 RNA-Seq samples across at least 54 tissues and 2630 individuals.  The ROP accounts for 99.9% of all reads, compared to 82.9% by conventional mapping-based protocols. ROP is able to profile: 
- repeats
- hyper-edited RNAs
- circRNAs, gene fusions, trans-splicing events
- recombined B and T cell receptor repertoires
- microbial communities

The 'dumpster diving' profile of unmapped reads output by our method is not limited to RNA-Seq technology and may be applied to whole-exome and whole-genome sequencing.


# Preprint 


Mangul, Serghei, et al. "[Comprehensive analysis of RNA-sequencing to find the source of 1 trillion reads across diverse adult human tissues](http://biorxiv.org/content/early/2017/06/12/053041)" bioRxiv (2016): 053041.

# Presentations

Slides from my talk at ASHG are available [here](https://sergheimangul.files.wordpress.com/2016/10/ashg2016_public.pdf) 

# Contact Us

Please do not hesitate to contact us (smangul@ucla.edu) if you have any comments, suggestions, or clarification requests regarding the tutorial or if you would like to contribute to this resource.

# Licence

See the [LICENCE.txt](https://github.com/smangul1/rop/blob/master/LICENSE.txt) file



