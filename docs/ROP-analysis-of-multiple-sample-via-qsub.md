In case you have many files you can submit multiple jobs using qsub.


Run ROP for all the samples at once :

- Run a single sample from the command line to make sure there are no unexpected errors. 

- Create a text file with all of the file names (no extensions) that will be run
```
ls *fasta | awk -F ".fasta" '{print $1}'>sample.txt 
```

- Create a run.sh file for each sample that will be run. 

```
#general example
while read line; do echo "<your command here>" > run_${line}.sh; done < <path to file with all sample names>

#implementation for ROP tool
while read line; do echo "python /u/home/s/serghei/code2/rop/rop.py --qsubArray --skipQC /u/home/b/brigitta/scratch/gtex/data/${line}.fasta $PWD/${line}" > run_${line}.sh; done<../sample.txt


```



Please follow the following steps:



```
ls run*sh | awk '{i+=1;print "qsub -cwd -V -N rop"i" -l h_data=16G,time=24:00:00 "$1}' > all.sh
```

Now you have `all.sh` with individual qsub command per sample. To run all the qsubs from `all.sh` do the following:

- chmod u+x all.sh
- nohup ./all.sh &