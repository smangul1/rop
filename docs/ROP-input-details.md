ROP accepts : `.fastq`, '.fastq', and `.bam`. Input can be in gzip format :  `.fastq.gz` and `.fasta.gz`

To run ROP for `.fastq.gz` input use `--gzip` option as follows:

```
python rop.py --fastqGz example/unmappedExample.fastq.gz test
``` 

TO run ROP for `.fasta.gz` skipping STEP1-2 use `--gzip` and `--skipPreliminary` options as follows:

```
python rop.py --gzip --skipPreliminary example/unmappedExampleFasta.fasta.gz testNew1/
```
