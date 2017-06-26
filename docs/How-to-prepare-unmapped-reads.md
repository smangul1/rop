In general case, you need to map the reads with any of available high-throughput aligners (e.g. STAR, tophat2) and save unmapped reads in .bam (binary format, requires less space) or .fastq (text format). 



Given the bam file with both mapped and unmapped reads use the following command to extract the unmapped reads:

```
samtools view -f 0x4 -bh  all.bam | samtools bam2fq - >unmapped.fastq
```

In case you need to extract each read from the pair into to a separate file  use the following commands:

```
samtools view -uf64 TCGA-CZ-4862.bam |samtools bam2fq - | gzip >x_1.fq.gz
samtools view -uf128 TCGA-CZ-4862.bam |samtools bam2fq - |gzip >x_2.fq.gz
```


The benchmarking of RNA aligners is provided in the recent Nature Methods study: Baruzzo, Giacomo, et al. "Simulation-based comprehensive benchmarking of RNA-seq aligners." Nature methods 14.2 (2017): 135-139.
