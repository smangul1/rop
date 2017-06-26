This is relevant to release ROP v1.0.6. We are planning to incorporate this is the next ROP release

Prior to moving all assembled CDR3s into one place, we need to rename full_cdr3.txt to incorporate the sample name.

```
while read line; do mv ${line}/antibodyProfile/full_cdr3.txt full_cdr3_${line}.txt;done<samples.txt
```


Starting from ROP v1.0.7 you can run rop-imrep which will consider reads mapped to BCR and TCr loci together with unmapped reads. 

To run across multiple samples, please use this command:

```
while read line ; do echo "~/code2/rop/rop-imrepGRCh38.sh /u/home/s/serghei/project/Catie/bam_sorted/HM-baseline1C_sort.bam $PWD/${line}/ \"\" ">run_${line}.sh;done<bam_sorted/samples.txt
```