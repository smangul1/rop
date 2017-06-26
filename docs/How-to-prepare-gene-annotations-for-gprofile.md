This is internal document. Use it in case you are planing to prepare gene annotations for a new organism. Note human and mouse annotations are already prepared. 

To extract the transcript names from gtf:

```
awk -F "transcript_id" '{print $2}' genes.gtf | awk -F "transcript_name" '{print $1}' | sed 's/"//g' | sed 's/;//' >transcripts.txt
```

To extract gene names from gtf:

```
awk -F "gene_id" '{print $2}' genes.gtf | awk -F "gene_name" '{print $1}' | sed 's/"//g' | sed 's/;//' >genes.txt
```

Merge them into a single file:

```
paste genes.txt transcripts.txt | awk '{print $1","$2}' | sort | uniq  >genes_transcripts.txt
```

Please download UT3, UTR5, and CDS from [here](https://genome.ucsc.edu/cgi-bin/hgTables).

Prepare them in the correct format:

```
awk -F "_" '{print $1}' CDS_NCBIM37.bed | awk '{if($2<0) print}' >CDS_NCBIM37_v2.bed
awk -F "_" '{print $1}' UTR5_NCBIM37.bed | awk '{if($2<0) print}' >UTR5_NCBIM37_v2.bed
awk -F "_" '{print $1}' UTR3_NCBIM37.bed | awk '{if($2<0) print}' >UTR3_NCBIM37_v2.bed
```

Now let's make a new bed file with the information of the gene each element belong to:

```
python ../../../scripts/prepareBed.py UTR3_NCBIM37_v2.bed genes_transcripts.txt ../bedPrepared/geneCoordinates_NCBIM37.bed ../bedPrepared/UTR3_NCBIM37_prepared.bed m

python ../../../scripts/prepareBed.py UTR5_NCBIM37_v2.bed genes_transcripts.txt ../bedPrepared/geneCoordinates_NCBIM37.bed ../bedPrepared/UTR5_NCBIM37_prepared.bed m

python ../../../scripts/prepareBed.py CDS_NCBIM37_v2.bed genes_transcripts.txt ../bedPrepared/geneCoordinates_NCBIM37.bed ../bedPrepared/CDS_NCBIM37_prepared.bed m
```
