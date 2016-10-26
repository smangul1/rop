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
			"bacteria": 0,
			"virus": 0,
			"ep": 0 }
nReadsImmuneTotal = 0

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
		write2Log("-- Identified " + str(nReads["repeat"]) + " lost repeat " +\
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

# temporary (until MiXCR switch)
INTFNS["afterNCLFasta"] = unmapped_file
branch_point_file = unmapped_file

if ARGS.immune:
	immuneReads = set()
	write2Log("5a. B lymphocytes profiling", LOGFNS["cmdLogfile"], "False")
	write2Log("5a. B lymphocytes profiling", LOGFNS["gLogfile"], ARGS.quiet)
	
	# IGH
	os.chdir(DIRS["igh"])
	cmd="ln -s %s//%s/antibody/internal_data/ ./" %(CD,DB_FOLDER)
	os.system(cmd)
	if ARGS.nonReductive or ARGS.qsub or ARGS.qsubArray:
		input_file = branch_point_file
	else:
		input_file = INTFNS["afterNCLFasta"]
	write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of immunoglobulin heavy locus (IGH)",LOGFNS["logIGH"],"False")
	write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",LOGFNS["logIGH"],"False")
	if ARGS.organism == "human":
		cmd = "%s/tools/igblastn -germline_db_V %s/%s/antibody/IGHV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGHJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05  2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(CD,CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file, LOGFNS["logIGH"],INTFNS["ighFile"])
	else:
		cmd = CD + "/tools/igblastn -organism mouse -auxiliary_data " + CD + "/" + DB_FOLDER + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/IGHV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGHJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05  2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file, LOGFNS["logIGH"],INTFNS["ighFile"])
	write2Log(cmd,LOGFNS["cmdLogfile"],"False")
	if ARGS.qsub or ARGS.qsubArray:
		f = open(RUNFNS["runIGHFile"],'w')
		f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(CD,DB_FOLDER))
		f.write(cmd+"\n")
		f.write("echo \"done!\"> %s/%s_igh.done \n" %(DIRS["igh"],BASENAME))
		f.close()
		if ARGS.qsub:
			if ARGS.maui:
				cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(RUNFNS["runIGHFile"])
			else:
				cmdQsub="qsub -cwd -V -N igh -l h_data=16G,time=24:00:00 %s" %(RUNFNS["runIGHFile"])
			os.system(cmdQsub)
			write2Log("Job for STEP5a(IGH) has been submitted via qsub",LOGFNS["gLogfile"], ARGS.quiet)
	else:
		os.chdir(DIRS["igh"])
		os.system(cmd)
		immuneReadsIGH=nReadsImmune(INTFNS["ighFile"])
		nReadsImmuneIGH=len(immuneReadsIGH)
		write2Log("--identified %s reads mapped to immunoglobulin heavy (IGH) locus" %(nReadsImmuneIGH) ,LOGFNS["gLogfile"],ARGS.quiet)

	#IGK	
	os.chdir(DIRS["igk"])
	cmd="ln -s %s//%s/antibody/internal_data/ ./" %(CD,DB_FOLDER)
	os.system(cmd)
	write2Log("IgBlast	was	used	to	identify	reads	spanning	VJ recombinations of immunoglobulin kappa locus (IGK)",LOGFNS["logIGK"],"False")
	write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",LOGFNS["logIGK"],"False")
	if ARGS.organism == "human":
		cmd = "%s/tools/igblastn -germline_db_V %s/%s/antibody/IGKV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGKJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05  2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(CD,CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file, LOGFNS["logIGK"],INTFNS["igkFile"])
	else:
		cmd = CD + "/tools/igblastn -organism mouse -auxiliary_data " + CD + "/" + DB_FOLDER + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/IGKV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGKJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05  2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file, LOGFNS["logIGK"],INTFNS["igkFile"])
	write2Log(cmd,LOGFNS["cmdLogfile"],"False")
	if ARGS.qsub or ARGS.qsubArray:
		f = open(RUNFNS["runIGKFile"],'w')
		f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(CD,DB_FOLDER))
		f.write(cmd+"\n")
		f.write("echo \"done!\"> %s/%s_igk.done \n" %(DIRS["igk"],BASENAME))
		f.close()
		if ARGS.qsub:
			if ARGS.maui:
				cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(RUNFNS["runIGKFile"])
			else:
				cmdQsub="qsub -cwd -V -N igk -l h_data=16G,time=24:00:00 %s" %(RUNFNS["runIGKFile"])
			os.system(cmdQsub)
			write2Log("Job for STEP5a(IGK) has been submitted via qsub",LOGFNS["gLogfile"], ARGS.quiet)
	else:
		os.chdir(DIRS["igk"])
		os.system(cmd)
		immuneReadsIGK=nReadsImmune(INTFNS["igkFile"])
		nReadsImmuneIGK=len(immuneReadsIGK)
		write2Log("--identified %s reads mapped to immunoglobulin kappa (IGK) locus " %(nReadsImmuneIGK) ,LOGFNS["gLogfile"],ARGS.quiet)
				
	# IGL	
	os.chdir(DIRS["igl"])
	cmd="ln -s %s//%s/antibody/internal_data/ ./" %(CD,DB_FOLDER)
	os.system(cmd)
	write2Log("IgBlast	was	used	to	identify	reads	spanning	VJ recombinations of immunoglobulin lambda locus (IGL)",LOGFNS["logIGL"],"False")
	write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",LOGFNS["logIGL"],"False")
	if ARGS.organism == "human":
		cmd = "%s/tools/igblastn -germline_db_V %s/%s/antibody/IGLV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGLJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(CD,CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file,LOGFNS["logIGL"],INTFNS["iglFile"])
	else:
		cmd = CD + "/tools/igblastn -organism mouse -auxiliary_data " + CD + "/" + DB_FOLDER + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/IGLV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGLJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file,LOGFNS["logIGL"],INTFNS["iglFile"])
	write2Log(cmd,LOGFNS["cmdLogfile"],"False")
	if ARGS.qsub or ARGS.qsubArray:
		f = open(RUNFNS["runIGLFile"],'w')
		f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(CD,DB_FOLDER))
		f.write(cmd+"\n")
		f.write("echo \"done!\">%s/%s_igl.done \n" %(DIRS["igl"],BASENAME))
		f.close()
		if ARGS.qsub:
			if ARGS.maui:
				cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(RUNFNS["runIGLFile"])
			else:
				cmdQsub="qsub -cwd -V -N igl -l h_data=16G,time=24:00:00 %s" %(RUNFNS["runIGLFile"])
			os.system(cmdQsub)
			write2Log("Job for STEP5a(IGL) has been submitted via qsub",LOGFNS["gLogfile"], ARGS.quiet)
	else:
		os.chdir(DIRS["igl"])
		os.system(cmd)
		immuneReadsIGL=nReadsImmune(INTFNS["iglFile"])
		nReadsImmuneIGL=len(immuneReadsIGL)
		write2Log("--identified %s reads mapped to immunoglobulin lambda (IGL) locus" %(nReadsImmuneIGL) ,LOGFNS["gLogfile"],ARGS.quiet)

	write2Log("5b. T lymphocytes profiling...",LOGFNS["cmdLogfile"],"False")
	write2Log("5b. T lymphocytes profiling...",LOGFNS["gLogfile"],ARGS.quiet)
	
	# TCRA	
	os.chdir(DIRS["tcra"])
	cmd="ln -s %s//%s/antibody/internal_data/ ./" %(CD,DB_FOLDER)
	os.system(cmd)
	write2Log("IgBlast	was	used	to	identify	reads	spanning	VJ recombinations of T cell receptor alpha locus (TCRA)",LOGFNS["logTCRA"],"False")
	write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",LOGFNS["logTCRA"],"False")
	if ARGS.organism == "human":
		cmd = "%s/tools/igblastn -ig_seqtype TCR -germline_db_V %s/%s/antibody/TRAV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRAJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(CD,CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file,LOGFNS["logTCRA"],INTFNS["tcraFile"])
	else:
		cmd = CD + "/tools/igblastn -organism mouse -ig_seqtype TCR -auxiliary_data " + CD + "/" + DB_FOLDER + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/TRAV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRAJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file,LOGFNS["logTCRA"],INTFNS["tcraFile"])
	write2Log(cmd,LOGFNS["cmdLogfile"],"False")
	if ARGS.qsub or ARGS.qsubArray:
		f = open(RUNFNS["runTCRAFile"],'w')
		f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(CD,DB_FOLDER))
		f.write(cmd+"\n")
		f.write("echo \"done!\">%s/%s_tcra.done \n"%(DIRS["tcra"],BASENAME))
		f.close()
		if ARGS.qsub:
			if ARGS.maui:
				cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(RUNFNS["runTCRAFile"])
			else:
				cmdQsub="qsub -cwd -V -N tcra -l h_data=16G,time=24:00:00 %s" %(RUNFNS["runTCRAFile"])
			os.system(cmdQsub)
			write2Log("Job for STEP5b(TCRA) has been submitted via qsub",LOGFNS["gLogfile"], ARGS.quiet)
	else:
		os.chdir(DIRS["tcra"])
		os.system(cmd)
		immuneReadsTCRA=nReadsImmune(INTFNS["tcraFile"])
		nReadsImmuneTCRA=len(immuneReadsTCRA)
		write2Log("--identified %s reads mapped to T cell receptor alpha (TCRA) locus" %(nReadsImmuneTCRA) ,LOGFNS["gLogfile"],ARGS.quiet) 

	# TCRB
	os.chdir(DIRS["tcrb"])
	cmd="ln -s %s//%s/antibody/internal_data/ ./" %(CD, DB_FOLDER)
	os.system(cmd)
	write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of T cell receptor beta locus (TCRB)",LOGFNS["logTCRB"],"False")
	write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",LOGFNS["logTCRB"],"False")
	if ARGS.organism == "human":
		cmd = "%s/tools/igblastn -ig_seqtype TCR -germline_db_V %s/%s/antibody/TRBV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRBJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(CD,CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file,LOGFNS["logTCRB"],INTFNS["tcrbFile"])
	else:
		cmd = CD + "/tools/igblastn -organism mouse -ig_seqtype TCR -auxiliary_data " + CD + "/" + DB_FOLDER + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/TRBV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRBJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file,LOGFNS["logTCRB"],INTFNS["tcrbFile"])
	write2Log(cmd,LOGFNS["cmdLogfile"],"False")
	if ARGS.qsub or ARGS.qsubArray:
		f = open(RUNFNS["runTCRBFile"],'w')
		f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(CD, DB_FOLDER))
		f.write(cmd+"\n")
		f.write("echo \"done!\">%s/%s_tcrb.done \n"%(DIRS["tcrb"],BASENAME))
		f.close()
		if ARGS.qsub:
			if ARGS.maui:
				cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(RUNFNS["runTCRBFile"])
			else:
				cmdQsub="qsub -cwd -V -N tcrb -l h_data=16G,time=24:00:00 %s" %(RUNFNS["runTCRBFile"])
			os.system(cmdQsub)
			write2Log("Job for STEP5b(TCRB) has been submitted via qsub",LOGFNS["gLogfile"], ARGS.quiet)
	else:
		os.chdir(DIRS["tcrb"])
		os.system(cmd)
		immuneReadsTCRB=nReadsImmune(INTFNS["tcrbFile"])
		nReadsImmuneTCRB=len(immuneReadsTCRB)
		write2Log("--identified %s reads mapped to T cell receptor beta (TCRB) locus" %(nReadsImmuneTCRB) ,LOGFNS["gLogfile"],ARGS.quiet)
		
	#TCRD	
	os.chdir(DIRS["tcrd"])
	cmd="ln -s %s//%s/antibody/internal_data/ ./" %(CD, DB_FOLDER)
	os.system(cmd)
	write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of T cell receptor delta locus (TCRD)",LOGFNS["logTCRD"],"False")
	write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",LOGFNS["logTCRD"],"False")
	if ARGS.organism == "human":
		cmd = "%s/tools/igblastn -ig_seqtype TCR -germline_db_V %s/%s/antibody/TRDV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRDJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(CD,CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file,LOGFNS["logTCRD"],INTFNS["tcrdFile"])
	else:
		cmd = CD + "/tools/igblastn -organism mouse -ig_seqtype TCR -auxiliary_data " + CD + "/" + DB_FOLDER + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/TRDV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRDJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file,LOGFNS["logTCRD"],INTFNS["tcrdFile"])
	write2Log(cmd,LOGFNS["cmdLogfile"],"False")
	if ARGS.qsub or ARGS.qsubArray:
		f = open(RUNFNS["runTCRDFile"],'w')
		f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(CD, DB_FOLDER))
		f.write(cmd+"\n")
		f.write("echo \"done!\">%s/%s_tcrd.done \n" %(DIRS["tcrd"], BASENAME))
		f.close()
		if ARGS.qsub:
			if ARGS.maui:
				cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(RUNFNS["runTCRDFile"])
			else:
				cmdQsub="qsub -cwd -V -N tcrd -l h_data=16G,time=24:00:00 %s" %(RUNFNS["runTCRDFile"])
			os.system(cmdQsub)
			write2Log("Job for STEP5b(TCRD) has been submitted via qsub",LOGFNS["gLogfile"], ARGS.quiet)
	else:
		os.chdir(DIRS["tcrd"])
		os.system(cmd)
		immuneReadsTCRD=nReadsImmune(INTFNS["tcrdFile"])
		nReadsImmuneTCRD=len(immuneReadsTCRD)
		write2Log("--identified %s reads mapped to T cell receptor delta (TCRD) locus" %(nReadsImmuneTCRD) ,LOGFNS["gLogfile"],ARGS.quiet)
				
	# TCRG	
	os.chdir(DIRS["tcrg"])
	cmd="ln -s %s//%s/antibody/internal_data/ ./" %(CD, DB_FOLDER)
	os.system(cmd)
	write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of T cell receptor gamma locus (TCRG)",LOGFNS["logTCRG"],"False")
	write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",LOGFNS["logTCRG"],"False")
	if ARGS.organism == "human":
		cmd = "%s/tools/igblastn -ig_seqtype TCR -germline_db_V %s/%s/antibody/TRGV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRGJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(CD,CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file,LOGFNS["logTCRG"],INTFNS["tcrgFile"])
	else:
		cmd = CD + "/tools/igblastn -organism mouse -ig_seqtype TCR -auxiliary_data " + CD + "/" + DB_FOLDER + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/TRGV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRGJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(CD,DB_FOLDER,CD,DB_FOLDER,CD,DB_FOLDER,input_file,LOGFNS["logTCRG"],INTFNS["tcrgFile"])
	write2Log(cmd,LOGFNS["cmdLogfile"],"False")
	if ARGS.qsub or ARGS.qsubArray:
		f = open(RUNFNS["runTCRGFile"],'w')
		f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(CD, DB_FOLDER))
		f.write(cmd+"\n")
		f.write("echo \"done!\">%s/%s_tcrg.done \n" %(DIRS["tcrg"],BASENAME))
		f.close()
		if ARGS.qsub:
			if ARGS.maui:
				cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(RUNFNS["runTCRGFile"])
			else:
				cmdQsub="qsub -cwd -V -N tcrg -l h_data=16G,time=24:00:00 %s" %(RUNFNS["runTCRGFile"])
			os.system(cmdQsub)
			write2Log("Job for STEP5b(TCRG) has been submitted via qsub",LOGFNS["gLogfile"], ARGS.quiet)
	else:
		os.chdir(DIRS["tcrg"])
		os.system(cmd)
		immuneReadsTCRG=nReadsImmune(INTFNS["tcrgFile"])
		nReadsImmuneTCRG=len(immuneReadsTCRG)
		write2Log("--identified %s reads mapped to T cell receptor gamma locus (TCRG) locus" %(nReadsImmuneTCRG) ,LOGFNS["gLogfile"],ARGS.quiet)
	nReadsImmuneTotal=0
	if not ARGS.qsub and not ARGS.qsubArray:
		nReadsImmuneTotal=nReadsImmuneIGH+nReadsImmuneIGL+nReadsImmuneIGK+nReadsImmuneTCRA+nReadsImmuneTCRB+nReadsImmuneTCRD+nReadsImmuneTCRG
		write2Log("In total: %s reads mapped to antibody repertoire loci" %(nReadsImmuneTotal) ,LOGFNS["gLogfile"],ARGS.quiet)
		write2Log("***Note: Combinatorial diversity of the antibody repertoire (recombinations of the of VJ gene segments)  will be available in the next release.",LOGFNS["gLogfile"],ARGS.quiet)
		write2File("done!",ARGS.dir+"/step5_antibodyProfile.done")
		immuneReads=set().union(immuneReadsTCRA,immuneReadsTCRB,immuneReadsTCRD,immuneReadsTCRG)
		excludeReadsFromFasta(input_file,immuneReads,INTFNS["afterImmuneFasta"])
		if not ARGS.dev:
			if ARGS.circRNA:		   
				os.remove(INTFNS["afterNCLFasta"])
else:
	print "5a. B lymphocytes profiling is skipped."
	print "5b. T lymphocytes profiling is skipped."

# temporary (until MiXCR switch)
if os.path.exists(INTFNS["afterImmuneFasta"]):
	unmapped_file = INTFNS["afterImmuneFasta"]
	
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
	  str(nReadsImmuneTotal + sum(nReads.values())) + " reads", 
	  LOGFNS["gLogfile"], ARGS.quiet)
	write2Log("***Unaccounted reads (not explained by ROP) are saved to " +\
	  INTFNS["unaccountedReadsFasta"], LOGFNS["gLogfile"], ARGS.quiet)
	write2Log("***Log file with all the commands used is available here: " +\
	  LOGFNS["cmdLogfile"], LOGFNS["gLogfile"], ARGS.quiet)
	tLog = ARGS.dir + "/" + "numberReads_" + BASENAME + ".log"
	write2Log("sample,totalUnmapped,nReads[\"LowQ\"],nReads[\"LowC\"]," +\
	  "nReads[\"rRNA\"],nReads[\"lost\"],nReads[\"repeat\"]," +\
	  "nReads[\"NCL\"],nReadsImmuneTotal,nMicrobialReads", tLog, True)
	write2Log(BASENAME + "," + str(n) + "," + str(nReads["LowQ"]) + "," +\
	  str(nReads["LowC"]) + "," + str(nReads["rRNA"]) + "," +\
	  str(nReads["lost"]) + "," + str(nReads["repeat"]) + "," +\
	  str(nReads["NCL"]) + "," + str(nReadsImmuneTotal) + "," +\
	  str(nReads["bacteria"] + nReads["virus"] + nReads["ep"]), tLog, True)

write2Log("""The list of the tools used by ROP and the paramemers is provided below.
************
**We have used FastQC (version 0.0.13, with the default parameters) downloaded from  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ to filter out low quality reads"
**We have used SEQLEAN (seqclean-x86_64, with the default parameters) downloaded from https://sourceforge.net/projects/seqclean/ to filter out low complexity reads
**We have used Megablast (BLAST+ version 2.2.30, with the following parameters: task = megablast, use_index = true; -outfmt 6 ;-evalue 1e-05; perc_identity = 100) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to filter out reads mapped to rRNA repeat sequence
**We have used Bowtie2 (version 2.0.5, with the following parameters: -k 1; -p 8; -f) downloaded from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml to identify lost reads mapped to reference transcriptome and genome
**We have used Megablast (BLAST+ version 2.2.30, with the following options: task=megablast, use_index=true, -outfmt 6 -evalue 1e-05, perc_identity	= 90) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to identify lost repeat reads mapped to database of repeat sequences (RepBase 20.07)
**We have used CIRI (version 1.2 with the following parameters: -S -I ) downloaded from https://sourceforge.net/projects/ciri/ to identify reads from circRNAs
**We have used IgBLAST (version v 1.4.0 with the following parameters: -germline_db_V;	germline_db_D; -germline_db_J; -outfmt 7 std qseq sseq; -evalue = 1e-05 ) downloaded from http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/igblast/release/1.4.0/ to identify immune reads spanningBCR/TCR receptor gene rearrangement in the variable domain (V(D)J recombinations)
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
