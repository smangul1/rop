Starting from release v1.0.7 we provide te possibility to assign multi-mapped reads (i.e reads mapped to multiple locations of the genome) into genomic categories. Use `--multi` and `--perCategory` options for `gprofile.py` and run `gprofilePlus.py` afterwards. 

```
usage: gprofile.py [-h] [--perCategory] [--mouse] [--multi] bam out
gprofile.py: error: too few arguments
```

File `out` contains statistics based on 

Number of reads per category reported in the `out` file is based on uniquely mapped reads (i.e. reads mapped to a single position in the genome). Multi-mapped reads are reported under `nMultiMapped` category. 

To assign multi-mapped reads into genomic categories you need to run `gprofilePlus.py`. It will randomly assign multi-mapped reads into genomic categories considering expression level of the genes. 