Identification of hyper-edited reads is not included in the current release of ROP. We are planning to include it in next the release.

Here we provide detailed instructions (prepared by Kevin Hsieh and Linus Chen) on how to identify hyper-edited reads using HE-pipeline downloaded from here:http://levanonlab.ls.biu.ac.il/resources/zip.


When running hyper-editing pipeline (HE-pipeline), additional changes can be made to parallelize the scripts for use with UCLA's Hoffman2 cluster. Before proceeding, follow the instructions in the `README` included with the scripts to prepare the reference and provide  the necessary third-party tools. Ensure that output directory is set correctly in `config_file.sh` (it is acceptable to use a single output directory), and that the list of input files has been prepared correctly.
  
###1. Modifications to script files
In `config_file.sh`, set `source_file="./_file_lists/file_$SGE_TASK_ID"`. This allows multiple instances to run simultaneously.   
  
In `run_hyper_editing.sh`, remove the last line `rm -r $unmap_dir $Trans_run_dir $analyseMM_dir` and insert `rm -r $unmap_dir/${out_pre}* $Trans_run_dir/${out_pre}* $analyseMM_dir/${out_pre}*` immediately before `done<$source_dir` so that intermediate files are deleted at the correct time and different instances do not interfere with each other.
  
###2. Submitting a job array
In order to submit a job array, we need to create separate input file lists for each instance. For that, this script, which we call `run.sh`, can be used:

	# usage: ./run.sh job_name config_file.sh
	
	mkdir _file_lists
	i=1
	while read sample; do
		echo Preparing for job $i: `printf "$sample" | sed "s .*\t  "`
		echo "$sample" > "./_file_lists/file_$i"
		((i++))
	done <./file_list
	echo "Submitting job-array $1."
	qsub -cwd -V -N $1 -l h_data=2G,time=24:00:00 -pe shared 8 -p -10 -t 1-$(($i-1)) -o ~ -e ~ $2

This script creates a directory called `_file_lists`, which it populates with individual file lists containing only one file each using the original input file list, `file_list`. It then submits `config_file.sh` as a job array using the `-t` option, with 8 cores, 8 × 2 GB = 16 GB of memory, 24 hours of time, and a lower priority so as not to interfere with other jobs. When invoking this script, specify the `job_name` and the name of `config_file.sh`, even if you have not changed it. Output will be sent to your home directory; change appropriately if you prefer output to go elsewhere. Beware that Hoffman2 may run several samples at once and intermediate files can grow very large, so we suggest having free space in your output directory equivalent to at least 50× the typical input size of all simultaneously running jobs (e.g., 50 × 150 simultaneous jobs × 500 MB typical input files = 3.75 TB).

###3. Analysis  

When using an analysis script (in this case, `Analyse_UE_basic.sh`), we suggest making an `analysis` folder within the Hyper Editing Scripts directory and using this script, which we call `run_analysis.sh`, to run the analysis:

	# usage: ./run_analysis.sh bed_directory analysis_file
	
	list=`ls $1 | grep bed_files | sed "s \.UE\..*  " | sed "s \.ES\..*  " | uniq`
	for sample in $list; do
		./Analyse_UE_basic.sh $1/$sample $2
	done

When invoking this script, specify the `bed_directory` (e.g., `../output/UEdetect.SE_0.05_0.6_30_0.6_0.1_0.8_0.2`) and the output `analysis_file` (e.g., `analysis_output.tsv`).