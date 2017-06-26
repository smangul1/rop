Here we describe how to run ROP for multiple samples using job array (hoffman2). Consider using this approach if the number of samples is higher than the number of jobs you are allowed to submit on the cluster. This tutorial about ROP analysis of multiple samples was prepared by  William Van (wvanderwey@ucla.edu ) during [B.I.G. Summer Program](http://qcb.ucla.edu/big-summer/) at QCBio, UCLA. 

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
- For each *run.sh*  file created, generate a command to run each file and store it on a separate line in a file called *myFunc.sh*

```
ls *sh | awk '{i+=1;print "if [ $1 == "i" ];then ./"$1" ;fi"}' > myFunc.sh
```

- Copy the file *myFuncFastWrapper.sh* from the ROP tool (/rop/source) to the current directory 

```
cp /u/home/s/serghei/code2/rop/myFuncFastWrapper.sh ./
```
- Change file permissions as necessary
```
#change file permissions for myFunc.sh
chmod 755 myFunc.sh

#change file permissions for all .sh files in current directory
chmod 755 *sh
```
- Find the total number of samples that will be run (number of lines in your myFunc.sh file)

```
wc -l myFunc.sh
```
- Lastly, use qsub to submit all of the samples

```
#general example
qsub -cwd -V -N <your job name> -l h_data=<memory needed per job>,express,time=<time needed per job> -t 1-<number of samples>:1 myFuncFastWrapper.sh

#implementation 
qsub -cwd -V -N rop -l h_data=12G,express,time=04:00:00 -t 1-176:1 myFuncFastWrapper.sh
```

