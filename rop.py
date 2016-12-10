print """********************************************************************************
ROP is a computational protocol aimed to discover the source of all reads, 
originated from complex RNA molecules, recombinant antibodies and microbial 
communities. 

Written by Serghei Mangul (smangul@ucla.edu), Harry Taegyun Yang 
(harry2416@gmail.com), Kevin Hsieh (kevin.hsieh@ucla.edu), and Linus Chen 
(u6.30cl@gmail.com). University of California, Los Angeles (UCLA). (c) 2016. 

Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see: https://sergheimangul.wordpress.com/rop/
ROP Tutorial: https://github.com/smangul1/rop/wiki 
********************************************************************************"""

from rop_functions import *

################################################################################
# Prepare for analysis

n, readLength, unmapped_file = prepare_for_analysis(ARGS.unmappedReads)
	# unmapped_file shall be passed reductively from step to step except when 
	# running non-reductively, in which case it shall not be modified after 
	# step 2

os.chdir(ARGS.dir)
write2Log("Processing " + str(n) + " unmapped reads of length " +\
  str(readLength) + ".", LOGFNS["gLogfile"], ARGS.quiet)
nReads = {	"LowQ": 0,
			"LowC": 0,
			"rRNA": 0,
			"lost": 0,
			"repeat": 0,
			"NCL": 0,
			"immune": 0,
			"bacteria": 0,
			"virus": 0,
			"ep": 0 }

################################################################################
# 1. Quality Control

if not ARGS.skipPreliminary and not ARGS.skipQC:
	write2Log("1. Quality Control...", LOGFNS["gLogfile"], ARGS.quiet)  
	write2Log("1. Quality Control...", LOGFNS["cmdLogfile"], True)
	os.chdir(DIRS["QC"])
	
	# 1a. lowQ
	if ARGS.skipLowq or not readsPresent("1a", unmapped_file): 
		write2Log("--low quality filtering step is skipped", LOGFNS["gLogfile"],
		  ARGS.quiet)
	else:
		nReads["LowQ"] = step_1a(unmapped_file, n)
		clean(unmapped_file)
		unmapped_file = INTFNS["lowQFileFasta"]
		write2Log("--filtered " + str(nReads["LowQ"]) + " low quality reads",
		  LOGFNS["gLogfile"], ARGS.quiet)
	
	# 1b. lowC
	if readsPresent("1b", unmapped_file):
		nReads["LowC"] = step_1b(unmapped_file)
		clean(unmapped_file)
		unmapped_file = INTFNS["lowQCFile"]
		write2Log("--filtered " + str(nReads["LowC"]) + " low complexity reads " +\
		  "(e.g. ACACACAC...)", LOGFNS["gLogfile"], ARGS.quiet)
		
	# 1c. rRNA
	if readsPresent("1c", unmapped_file):
		nReads["rRNA"], n_rRNATotal = step_1c(unmapped_file, readLength)
		clean(INTFNS["rRNAFile"])
		clean(unmapped_file)
		unmapped_file = INTFNS["afterrRNAFasta"]
		write2Log("--filtered " + str(nReads["rRNA"]) + " rRNA reads", 
		  LOGFNS["gLogfile"], ARGS.quiet)
	else:
		n_rRNATotal = 0
	
	# write results
	write2Log("In total: " + str(nReads["LowQ"] + nReads["LowC"] +\
	  nReads["rRNA"]) + " reads failed QC and are filtered out.", 
	  LOGFNS["gLogfile"], ARGS.quiet)
	write2Log("Number of entries in " + INTFNS["rRNAFile"] + " is " +\
	  str(n_rRNATotal), LOGFNS["cmdLogfile"], True)
	write2File("done!", ARGS.dir + "/step1_QC.done")
else:
	write2Log("1. Quality Control is skipped.", LOGFNS["gLogfile"], ARGS.quiet)

################################################################################
# 2. Remapping to reference

if not readsPresent("2", unmapped_file):
	ARGS.skipPreliminary = True

if not ARGS.skipPreliminary:
	write2Log("2. Remapping to reference...", LOGFNS["cmdLogfile"], True)
	write2Log("2. Remapping to reference...", LOGFNS["gLogfile"], ARGS.quiet)
	os.chdir(DIRS["lostReads"])
	
	nReads["lost"], lostReads0_len, lostReads1_len, lostReads2_len =\
	  step_2(unmapped_file)
	write2Log("--identified " + str(nReads["lost"]) + " lost reads from " +\
	  "unmapped reads. Among those: " +\
	  str(lostReads0_len) + " reads with 0 mismatches, " +\
	  str(lostReads1_len) + " reads with 1 mismatch, and " +\
	  str(lostReads2_len) + " reads with 2 mismatches", LOGFNS["gLogfile"], 
	  ARGS.quiet)
	write2Log("***Note: Complete list of lost reads is available from sam " +\
	  "files: " + INTFNS["gBamFile"] + ", " + INTFNS["tBamFile"], 
	  LOGFNS["gLogfile"], ARGS.quiet)
	clean(unmapped_file)
	unmapped_file = INTFNS["afterlostReadsFasta"]
	write2File("done!", ARGS.dir + "/step2_lostReads.done")
	
	if ARGS.clean:
		write2Log("Clean mode selected - removing analysis sam files.", 
		  LOGFNS["gLogfile"], ARGS.quiet)
		os.remove(INTFNS["gBamFile"])
		os.remove(INTFNS["tBamFile"])
else:
	write2Log("2. Remapping to reference is skipped.", LOGFNS["gLogfile"], 
	  ARGS.quiet)
		
if ARGS.nonReductive or ARGS.qsub or ARGS.qsubArray:
	write2Log("*******************************", LOGFNS["gLogfile"], ARGS.quiet)
	write2Log("Non-reductive mode is selected: Low quality, low complexity, " +\
	  "and rRNA reads and lost reads are filtered out. Resulting high " +\
	  "quality non-matching reads are provided as input for each of steps 3-6.", 
	  LOGFNS["gLogfile"], ARGS.quiet)
	write2Log("*******************************", LOGFNS["gLogfile"], ARGS.quiet)
	

################################################################################
# 3. Map repeat sequences

if not readsPresent("3", unmapped_file):
	ARGS.repeat = False

if ARGS.repeat:
	write2Log("3. Mapping repeat sequences...", LOGFNS["cmdLogfile"], True)
	write2Log("3. Mapping repeat sequences...", LOGFNS["gLogfile"], ARGS.quiet)
	os.chdir(DIRS["lostRepeat"])
	
	cmd = CD + "/tools/blastn -task megablast -index_name " + CD + "/" +\
	  DB_FOLDER + "/repeats/repbase.fa -use_index true -query " +\
	  unmapped_file + " -db " + CD + "/" + DB_FOLDER + "/repeats/repbase.fa " +\
	  "-outfmt 6 -evalue 1e-05 >" + INTFNS["repeatFile"] +\
	  " 2>" + LOGFNS["logLostRepeat"]
	if ARGS.qsub or ARGS.qsubArray:
		write2Log(cmd, RUNFNS["runLostRepeatFile"], True)
		write2Log("echo \"done!\" >" + DIRS["lostRepeat"] + "/" + BASENAME +\
		  "_lostRepeat.done", RUNFNS["runLostRepeatFile"], True)
		if ARGS.qsub:
			qsub("3", RUNFNS["runLostRepeatFile"])
	else:
		nReads["repeat"] = step_3(unmapped_file, readLength, cmd)
		write2Log("--identified " + str(nReads["repeat"]) + " lost repeat " +\
		  "sequences from unmapped reads.", LOGFNS["gLogfile"], ARGS.quiet)
		write2Log("***Note: Repeat sequences classification into classes " +\
		  "(e.g. LINE) and families (e.g. Alu) will be available in next " +\
		  "release", LOGFNS["gLogfile"], ARGS.quiet)
		if not ARGS.nonReductive:
			clean(unmapped_file)
			unmapped_file = INTFNS["afterlostRepeatFasta"]
		write2File("done!", ARGS.dir + "/step3_lostRepeatSequences.done")
else:
	write2Log("3. Mapping to repeat sequences is skipped.", LOGFNS["gLogfile"], 
	  ARGS.quiet)


################################################################################
# 4. Non-co-linear RNA profiling

if not readsPresent("4", unmapped_file):
	ARGS.circRNA = False

if ARGS.circRNA:
	write2Log("4. Non-co-linear RNA profiling", LOGFNS["cmdLogfile"], True)
	write2Log("4. Non-co-linear RNA profiling", LOGFNS["gLogfile"], ARGS.quiet)
	os.chdir(DIRS["NCL"])
	write2Log("***Note: Trans-spicing and gene fusions are currently not " +\
	  "supported.", LOGFNS["gLogfile"], ARGS.quiet)
	
	cmd = CD + "/tools/bwa mem -T -S " + CD + "/" + DB_FOLDER +\
	  "/BWAIndex/genome.fa " + unmapped_file + " >" + INTFNS["NCL_CIRI_file"] +\
	  " 2>" + LOGFNS["logNCL"] + " \n"
	cmd = cmd + "perl " + CD + "/tools/CIRI_v1.2.pl -S -I " +\
	  INTFNS["NCL_CIRI_file"] +" -O " + INTFNS["after_NCL_CIRI_file_prefix"] +\
	  " -F " + CD + "/" + DB_FOLDER + "/BWAIndex/genome.fa 1>>" +\
	  LOGFNS["logNCL"] + " 2>>" + LOGFNS["logNCL"]
	if ARGS.qsub or ARGS.qsubArray:
		write2Log(cmd, RUNFNS["runNCL_CIRIfile"], True)
		write2Log("echo \"done!\" >" + DIRS["NCL"] + "/" + BASENAME +\
		  "_NCL_CIRI.done", RUNFNS["runNCL_CIRIfile"], True)
		if ARGS.qsub:
			qsub("4", RUNFNS["runNCL_CIRIfile"])
	else:
		nReads["NCL"] = step_4(unmapped_file, cmd)
		write2Log("--identified " + str(nReads["NCL"]) + " reads from circRNA.", 
		  LOGFNS["gLogfile"], ARGS.quiet)
		write2Log("***Note: circRNAs detected by CIRI are available here: " +\
		  INTFNS["after_NCL_CIRI_file_prefix"], LOGFNS["gLogfile"], ARGS.quiet)
		if not ARGS.nonReductive:
			clean(unmapped_file)
			unmapped_file = INTFNS["afterNCLFasta"]
		write2File("done!", ARGS.dir + "/step4_NCL.done")
else:
	write2Log("4. Non-co-linear RNA profiling is skipped.", LOGFNS["gLogfile"], 
	  ARGS.quiet)

	
################################################################################
# 5. Lymphocyte profiling

if not readsPresent("5", unmapped_file):
	ARGS.immune = False
	
if ARGS.immune and ARGS.organism == "human":
	write2Log("5. Lymphocyte profiling", LOGFNS["cmdLogfile"], True)
	write2Log("5. Lymphocyte profiling", LOGFNS["gLogfile"], ARGS.quiet)
	os.chdir(DIRS["antibody"])
	
	cmd = CD + "/tools/imrep/imrep_wrapper.sh " + unmapped_file
	if ARGS.qsub or ARGS.qsubArray:
		write2Log(cmd, RUNFNS["runAntibodyFile"], True)
		write2Log("echo \"done!\" >" + DIRS["antibody"] + "/" + BASENAME +\
		  "_antibodyProfile.done", RUNFNS["runAntibodyFile"], True)
		if ARGS.qsub:
			qsub("5", RUNFNS["runAntibodyFile"])
	else:
		write2Log(cmd, LOGFNS["cmdLogfile"], True)
		nReads["immune"] = step_5(unmapped_file, cmd)
		if not ARGS.nonReductive:
			clean(unmapped_file)
			unmapped_file = INTFNS["afterImmuneFasta"]
		write2Log("In total: " + str(nReads["immune"]) + " (unique) reads " +\
		  "mapped to antibody repertoire loci.", LOGFNS["gLogfile"], ARGS.quiet)
		write2File("done!", ARGS.dir + "/step5_antibodyProfile.done")
elif ARGS.immune and ARGS.organism == "mouse":
	write2Log("5. Lymphocyte profiling is skipped (not supported for mouse).", 
	  LOGFNS["gLogfile"], ARGS.quiet)
else:
	write2Log("5. Lymphocyte profiling is skipped.", LOGFNS["gLogfile"], 
	  ARGS.quiet)


################################################################################
# 6. Microbiome profiling

if not readsPresent("6", unmapped_file):
	ARGS.microbiome = False

if ARGS.microbiome:
	write2Log("6. Microbiome profiling...", LOGFNS["cmdLogfile"], True)
	write2Log("6. Microbiome profiling...", LOGFNS["gLogfile"], ARGS.quiet)
	
	# 6a. metaphlan
	os.chdir(DIRS["metaphlan"])
	cmd = "python " + CD + "/tools/metaphlan2.py " + unmapped_file + " " +\
	  INTFNS["metaphlan_intermediate_map"] + " --mpa_pkl " + CD + "/" +\
	  DB_FOLDER + "/metaphlan/mpa_v20_m200.pkl --bowtie2_exe " + CD +\
	  "/tools/bowtie2 --input_type multifasta --bowtie2db " + CD + "/" +\
	  DB_FOLDER + "/metaphlan/mpa_v20_m200 -t reads_map --nproc 8 " +\
	  "--bowtie2out " + INTFNS["metaphlan_intermediate_bowtie2out"] + " >>" +\
	  LOGFNS["logMetaphlan"] + " 2>>" + LOGFNS["logMetaphlan"]
	cmd = cmd + "\n" + "python " + CD + "/tools/metaphlan2.py --mpa_pkl " +\
	  CD + "/" + DB_FOLDER + "/metaphlan/mpa_v20_m200.pkl --bowtie2_exe " +\
	  CD + "/tools/bowtie2 --input_type bowtie2out " +\
	  INTFNS["metaphlan_intermediate_bowtie2out"] + " -t rel_ab >" +\
	  INTFNS["metaphlan_output"] + " 2>>" + LOGFNS["logMetaphlan"]
	write2Log(cmd, LOGFNS["cmdLogfile"], True)
	if ARGS.qsub or ARGS.qsubArray:
		write2Log(cmd, RUNFNS["run_metaphlan_file"], True)
		write2Log("echo \"done!\" >" + DIRS["metaphlan"] + "/" + BASENAME +\
		  "_metaphlan.done", RUNFNS["run_metaphlan_file"], True)
		if ARGS.qsub:
			qsub("6a_metaphlan", RUNFNS["run_metaphlan_file"])
	else:
		if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
		write2Log("***Microbiome profiling by Metaphlan2: taxonomic profile " +\
		  "of microbial communities detected by Metaphlan2 is available " +\
		  "here: " + DIRS["metaphlan"], LOGFNS["gLogfile"], ARGS.quiet)
	
	# 6b. bacteria	
	os.chdir(DIRS["bacteria"])
	write2Log("Megablast was used to map the reads onto the	bacterial " +\
	  "reference genomes.", LOGFNS["logBacteria"], True)
	write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides " +\
	  "aligned are considered microbial reads and are saved into " +\
	  LOGFNS["bacteriaFileFiltered"], LOGFNS["logBacteria"], True)
	write2Log("---------------", LOGFNS["logBacteria"], True)
	cmd = CD + "/tools/blastn -task megablast -index_name " + CD + "/" +\
	  DB_FOLDER + "/bacteria/bacteria -use_index true -query " +\
	  unmapped_file + " -db " + CD + "/" + DB_FOLDER + "/bacteria/bacteria " +\
	  "-outfmt 6 -evalue 1e-05 >" + LOGFNS["bacteriaFile"] + " 2>>" +\
	  LOGFNS["logBacteria"]
	write2Log(cmd, LOGFNS["cmdLogfile"], True)
	if ARGS.qsub or ARGS.qsubArray:
		write2Log(cmd, RUNFNS["runBacteriaFile"], True)
		write2Log("echo \"done!\" >" + DIRS["bacteria"] + "/" + BASENAME +\
		  "_bacteria.done", RUNFNS["runBacteriaFile"], True)
		if ARGS.qsub:
			qsub("6b_bacteria", RUNFNS["runBacteriaFile"])
	else:
		nReads["bacteria"] = step_6b(unmapped_file, readLength, cmd)
		write2Log("--identified " + str(nReads["bacteria"]) + " reads " +\
		"mapped to bacterial genomes", LOGFNS["gLogfile"], ARGS.quiet)
		if not ARGS.nonReductive:
			clean(unmapped_file)
			unmapped_file = INTFNS["afterBacteriaFasta"]
	
	# 6c. virus	
	os.chdir(DIRS["virus"])
	write2Log("Megablast was used to map the reads onto the	viral reference	" +\
	  "genomes.", LOGFNS["logVirus"], True)
	write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides " +\
	  "aligned are considered microbial reads and are saved into " +\
	  LOGFNS["virusFileFiltered"], LOGFNS["logVirus"], True)
	write2Log("----------------------", LOGFNS["logVirus"], True)
	cmd = CD + "/tools/blastn -task megablast -index_name " + CD + "/" +\
	  DB_FOLDER + "/virus/viruses -use_index true -query " + unmapped_file +\
	  " -db " + CD + "/" + DB_FOLDER + "/virus/viruses -outfmt 6 -evalue " +\
	  "1e-05 >" + LOGFNS["virusFile"] + " 2>>" + LOGFNS["logVirus"]
	write2Log(cmd, LOGFNS["cmdLogfile"], True)
	if ARGS.qsub or ARGS.qsubArray:
		write2Log(cmd, RUNFNS["runVirusFile"], True)
		write2Log("echo \"done!\" >" + DIRS["virus"] + "/" + BASENAME +\
		  "_virus.done", RUNFNS["runVirusFile"], True)
		if ARGS.qsub:
			qsub("6c_virus", RUNFNS["runVirusFile"])
	else:
		nReads["virus"] = step_6c(unmapped_file, readLength, cmd)
		write2Log("--identified " + str(nReads["virus"]) + " reads mapped " +\
		"to viral genomes", LOGFNS["gLogfile"], ARGS.quiet)
		if not ARGS.nonReductive:
			clean(unmapped_file)
			unmapped_file = INTFNS["afterVirusFasta"]

	# 6d. eukaryotic pathogens	
	os.chdir(DIRS["eupathdb"])
	dbList = ["ameoba", "crypto", "giardia", "microsporidia", "piroplasma",
	  "plasmo", "toxo", "trich", "tritryp"]
	for db in dbList:
		afterFasta = DIRS["eupathdb"] + BASENAME + "_after_" + db + ".fasta"
		eupathdbFile = DIRS["eupathdb"] + BASENAME + "_" + db +\
		  "_blastFormat6.csv"
		eupathdbFileFiltered = DIRS["eupathdb"] + BASENAME + "_" + db +\
		  "Filtered_blastFormat6.csv"
		runEupathdbFile = DIRS["eupathdb"] + "/run_" + BASENAME + "_" + db +\
		  ".sh"
		write2Log("Megablast was used to map the reads onto the	reference " +\
		  "genomes of eukaryotes (http://eupathdb.org/eupathdb/).",
		  LOGFNS["logEukaryotes"], True)
		write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides " +\
		  "aligned are considered microbial reads and are saved into " +\
		  afterFasta, LOGFNS["logEukaryotes"], True)
		write2Log("-------------", LOGFNS["logEukaryotes"], True)
		cmd = CD + "/tools/blastn -task megablast -index_name " + CD + "/" +\
		  DB_FOLDER + "/eupathdb/" + db + " -use_index true -query " +\
		  unmapped_file + " -db " + CD + "/" + DB_FOLDER + "/eupathdb/" + db +\
		  " -outfmt 6 -evalue 1e-05 >" + eupathdbFile + " 2>>" +\
		  LOGFNS["logEukaryotes"]
		write2Log(cmd, LOGFNS["cmdLogfile"], True)
		if ARGS.qsub or ARGS.qsubArray:
			write2Log(cmd, runEupathdbFile, True)
			write2Log("echo \"done!\"  >" + DIRS["eupathdb"] +\
			  "/" + DIRS["eupathdb"] + ".done", runEupathdbFile, True)
			if ARGS.qsub:
				qsub("6d_" + db, runEupathdbFile)
		else:
			nEupathdbReads = step_6d(unmapped_file, readLength, cmd, 
			  eupathdbFile, eupathdbFileFiltered, afterFasta)
			write2Log("--identified " + str(nEupathdbReads) + " reads " +\
			  "mapped to " + db + " genomes", LOGFNS["gLogfile"], ARGS.quiet)
			nReads["ep"] = nEupathdbReads
			if not ARGS.nonReductive:
				clean(unmapped_file)
				unmapped_file = afterFasta
	if not (ARGS.qsub or ARGS.qsubArray):
		os.rename(DIRS["eupathdb"] + BASENAME + "_after_tritryp.fasta", 
		  INTFNS["unaccountedReadsFasta"])
		write2Log("In total: " + str(nReads["bacteria"] + nReads["virus"] +\
		  nReads["ep"]) + " reads mapped to microbial genomes", 
		  LOGFNS["gLogfile"], ARGS.quiet)
		write2File("done!", ARGS.dir + "/step6_microbiomeProfile.done")
else:
	write2Log("6.  Microbiome profiling is skipped.", LOGFNS["gLogfile"], 
	  ARGS.quiet)

################################################################################
# Wrap-up

if not ARGS.qsubArray and not ARGS.qsub:
	write2Log("Summary: The ROP protocol is able to account for " +\
	  str(sum(nReads.values())) + " reads", 
	  LOGFNS["gLogfile"], ARGS.quiet)
	write2Log("***Unaccounted reads (not explained by ROP) are saved to " +\
	  INTFNS["unaccountedReadsFasta"], LOGFNS["gLogfile"], ARGS.quiet)
	write2Log("***Log file with all the commands used is available here: " +\
	  LOGFNS["cmdLogfile"], LOGFNS["gLogfile"], ARGS.quiet)
	tLog = ARGS.dir + "/" + "numberReads_" + BASENAME + ".log"
	write2Log("sample,totalUnmapped,nReads[\"LowQ\"],nReads[\"LowC\"]," +\
	  "nReads[\"rRNA\"],nReads[\"lost\"],nReads[\"repeat\"]," +\
	  "nReads[\"NCL\"],nReads[\"immune\"],nMicrobialReads", tLog, True)
	write2Log(BASENAME + "," + str(n) + "," + str(nReads["LowQ"]) + "," +\
	  str(nReads["LowC"]) + "," + str(nReads["rRNA"]) + "," +\
	  str(nReads["lost"]) + "," + str(nReads["repeat"]) + "," +\
	  str(nReads["NCL"]) + "," + str(nReads["immune"]) + "," +\
	  str(nReads["bacteria"] + nReads["virus"] + nReads["ep"]), tLog, True)

write2Log("""The list of the tools used by ROP is provided below.
************
**We have used FastQC (version 0.0.13, with the default parameters) downloaded from  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ to filter out low quality reads"
**We have used SEQLEAN (seqclean-x86_64, with the default parameters) downloaded from https://sourceforge.net/projects/seqclean/ to filter out low complexity reads
**We have used Megablast (BLAST+ version 2.2.30, with the following parameters: task = megablast, use_index = true; -outfmt 6 ;-evalue 1e-05; perc_identity = 100) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to filter out reads mapped to rRNA repeat sequence
**We have used Bowtie2 (version 2.0.5, with the following parameters: -k 1; -p 8; -f) downloaded from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml to identify lost reads mapped to reference transcriptome and genome
**We have used Megablast (BLAST+ version 2.2.30, with the following options: task=megablast, use_index=true, -outfmt 6 -evalue 1e-05, perc_identity	= 90) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to identify lost repeat reads mapped to database of repeat sequences (RepBase 20.07)
**We have used CIRI (version 1.2 with the following parameters: -S -I ) downloaded from https://sourceforge.net/projects/ciri/ to identify reads from circRNAs
**We have used ImReP (version 2.0) to identify immune reads spanning BCR/TCR receptor gene rearrangement in the variable domain (V(D)J recombinations)
**We have used Megablast (BLAST+ version 2.2.30, with the following parameters: task=megablast, use_index=true; -outfmt 6 ;-evalue 1e-05) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to identify microbial reads mapped onto the microbial genomes (bacteria, viruses, and eukaryotic pathogens)
**We have used Metaphlan2 (version 2.2.0, with the following parameters: --mpa_pkl; --bowtie2_exe; --input_type multifasta; --bowtie2db ; -t reads_map/rel_ab ; --nproc 8) downloaded from https://bitbucket.org/biobakery/metaphlan2 to obtain taxonomic profile of microbial communities
************
For more information about the paramemers and databases used by ROP, please see the preprint: Dumpster diving in RNA-sequencing to find the source of every last read http://biorxiv.org/content/early/2016/05/13/053041"
********************""", LOGFNS["toolsLogfile"], True)

write2Log("Important: ROP relies on several open source tools that were " +\
  "developed by other groups. These components are (c) their respective " +\
  "developers and are redistributed with ROP to provide ease-of-use. The " +\
  "list of the tools used by ROP and the parameters/reference databases are " +\
  "provided here: " + LOGFNS["toolsLogfile"], LOGFNS["gLogfile"], ARGS.quiet)
