Here we describe how to run ROP for multiple samples using qsub (hoffman2). Consider using this approach if the number of samples is high. This is also alternative to getting an interactive session node for each sample. 

Prior to running qsub complete the following steps: 

* Run a single sample from the command line to make sure there are no unexpected errors.

* Create a text file with all of the file names (no extensions) that will be run


```
ls *fastq | awk -F ".fastq" '{print $1}'>sample.txt 
```

* Create a run.sh file for each sample that will be run.

```
#general example
while read line; do echo "<your command here>" > run_${line}.sh; done < <path to file with all sample names>

#implementation for ROP tool
while read line; do echo "python /u/home/s/serghei/code2/rop/rop.py /u/home/b/brigitta/scratch/gtex/data/${line}.fastq /u/home/b/brigitta/scratch/gtex/out/${line}" > run_${line}.sh; done<sample.txt

```

Where 
* /u/home/s/serghei/code2/ - directory where rop is installed
* /u/home/b/brigitta/scratch/gtex/data/ - dir where fastq files are 
* /u/home/b/brigitta/scratch/gtex/out/ - dir where you want to save the results



For each run.sh file created, generate a command to run each file and store it on a separate line in a file called all.sh

```
ls run*sh | awk '{i+=1;print "qsub -cwd -V -N name"i" -l h_data=16G,time=24:00:00 "$1}' > all.sh
```

where 
* name is the name of the jobs
* 16G is the amount of memory you would like to use for each job
* 24:00:00 is the time you would like to require for each job. 

To run all.sh you can use 

```
nohup all.sh &
```