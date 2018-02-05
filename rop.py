print ("""********************************************************************************
ROP (v1.0.8) is a computational protocol aimed to discover the source of all reads,
originated from complex RNA molecules, recombinant B and T cell receptors and microbial
communities. 

Written by Serghei Mangul (smangul@ucla.edu), Harry Taegyun Yang 
(harry2416@gmail.com), Kevin Hsieh (kevin.hsieh@ucla.edu), and Linus Chen 
(u6.30cl@gmail.com). University of California, Los Angeles (UCLA). (c) 2016. 

Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see: https://github.com/smangul1/rop/wiki
********************************************************************************""")





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
    os.chdir(DIRS["QC"])
	
    # 1a. lowQ
    if ARGS.skipLowq or not readsPresent("1a", unmapped_file):
        write2Log("--low quality filtering step is skipped", LOGFNS["gLogfile"],ARGS.quiet)
    else:
        nReads["LowQ"] = step_1a(unmapped_file, n)
        clean(unmapped_file)
        unmapped_file = INTFNS["lowQFileFasta"]
        write2Log("--identified " + str(nReads["LowQ"]) + " low quality reads. Those reads are marked as lowQuality in the read name and are used in the donwstream analyis",LOGFNS["gLogfile"], ARGS.quiet)

		
    # 1c. rRNA
    if readsPresent("1c", unmapped_file):
        nReads["rRNA"], n_rRNATotal = step_1c(unmapped_file, readLength)
        clean(INTFNS["rRNAFile"])
        clean(unmapped_file)
        unmapped_file = INTFNS["afterrRNAFasta"]
        write2Log("--filtered " + str(nReads["rRNA"]) + " rRNA reads", LOGFNS["gLogfile"], ARGS.quiet)
    else:
        n_rRNATotal = 0
	
    # write results
    #write2Log("In total: " + str(nReads["rRNA"]) + " reads failed QC and are filtered out.", LOGFNS["gLogfile"], ARGS.quiet)
    #write2Log("In total: " + strn(Reads["LowQ"]) + " are marked as low quality reads",LOGFNS["gLogfile"], ARGS.quiet)



    #write2Log("Number of entries in " + INTFNS["rRNAFile"] + " is " +str(n_rRNATotal), LOGFNS["cmdLogfile"], True)
    #write2File("done!", ARGS.dir + "/step1_QC.done")
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
    print (ARGS.max,ARGS.pe)
    
    nReads["lost"], lostReads0_len, lostReads1_len, lostReads2_len =step_2(unmapped_file,ARGS.max,ARGS.pe,readLength)
    write2Log("--identified " + str(nReads["lost"]) + " lost reads from " + "unmapped reads. Among those: " +str(lostReads0_len) + " reads with 0 mismatches, " +str(lostReads1_len) + " reads with 1 mismatch, and " +str(lostReads2_len) + " reads with 2 mismatches", LOGFNS["gLogfile"],ARGS.quiet)
    write2Log("***Note: Complete list of lost reads is available from sam " +"files: " + INTFNS["gBamFile"] + ", " + INTFNS["tBamFile"],LOGFNS["gLogfile"], ARGS.quiet)
    clean(unmapped_file)
    unmapped_file = INTFNS["afterlostReadsFasta"]
    write2File("done!", ARGS.dir + "/step2_lostReads.done")
	
    if ARGS.clean:
	    write2Log("Clean mode selected - removing analysis sam files.", LOGFNS["gLogfile"], ARGS.quiet)
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
		nReads["repeat"] = step_3(unmapped_file, readLength, cmd,ARGS.max,ARGS.pe)
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


    command=read_commands()
    cmdNCL=command[2]
    
    print (RUNFNS["runNCL_CIRIfile"])
    
    cmd = cmdNCL +" "+ unmapped_file +" 2>" + LOGFNS["logNCL"] + " \n"
    cmd+=CD+"/tools/samtools bam2fq accepted_hits.bam >accepted_hits.fastq"
    
    
    
    
    print (cmd)

    if ARGS.qsub or ARGS.qsubArray:
	    write2Log(cmd, RUNFNS["runNCL_CIRIfile"], True)
	    write2Log("echo \"done!\" >" + DIRS["NCL"] + "/" + BASENAME +"_NCL_CIRI.done", RUNFNS["runNCL_CIRIfile"], True)
	    if ARGS.qsub:
		    qsub("4", RUNFNS["runNCL_CIRIfile"])
    else:
	    nReads["NCL"] = step_4(unmapped_file, cmd,ARGS.pe)
	    write2Log("--identified " + str(nReads["NCL"]) + " reads from NCLs.", LOGFNS["gLogfile"], ARGS.quiet)
	    write2Log("***Note: circRNAs detected by CIRI are available here: " +INTFNS["after_NCL_CIRI_file_prefix"], LOGFNS["gLogfile"], ARGS.quiet)
	    if not ARGS.nonReductive:
		    clean(unmapped_file)
		    unmapped_file = INTFNS["afterNCLFasta"]
	    write2File("done!", ARGS.dir + "/step4_NCL.done")
else:
    write2Log("4. Non-co-linear RNA profiling is skipped.", LOGFNS["gLogfile"], ARGS.quiet)

	
################################################################################
# 5.BCR/TCR

if not readsPresent("5", unmapped_file):
    ARGS.immune = False
	
if ARGS.immune and ARGS.organism == "human":
    write2Log("5. T and B cell repetoires profiling", LOGFNS["cmdLogfile"], True)
    write2Log("5. T and B cell repetoires profiling", LOGFNS["gLogfile"], ARGS.quiet)
    os.chdir(DIRS["antibody"])
    cmd = " python "+ CD + "/prerequisite_software/imrep/imrep.py  -f 0 --extendedOutput " + unmapped_file +" imrep.cdr3 >imrep.log 2>>imrep.log"
    
    print ("******")
    print (unmapped_file)
    
    if ARGS.qsub or ARGS.qsubArray:
        write2Log(cmd, RUNFNS["runAntibodyFile"], True)
        write2Log("echo \"done!\" >" + DIRS["antibody"] + "/" + BASENAME +"_antibodyProfile.done", RUNFNS["runAntibodyFile"], True)
        if ARGS.qsub:
            qsub("5", RUNFNS["runAntibodyFile"])
    else:
        write2Log(cmd, LOGFNS["cmdLogfile"], True)
        nReads["immune"] = step_5(unmapped_file, cmd,ARGS.pe)
        if not ARGS.nonReductive:
            clean(unmapped_file)
            unmapped_file = INTFNS["afterImmuneFasta"]
        write2Log("In total: " + str(nReads["immune"]) + " reads " +"mapped to recombined TCRs rand BCRs,", LOGFNS["gLogfile"], ARGS.quiet)
        write2File("done!", ARGS.dir + "/step5_antibodyProfile.done")
elif ARGS.immune and ARGS.organism == "mouse":
	write2Log("5. T and B cell repetoires profiling is skipped (not supported for mouse).",
	  LOGFNS["gLogfile"], ARGS.quiet)
else:
	write2Log("5. T and B cell repetoires profiling is skipped.", LOGFNS["gLogfile"],
	  ARGS.quiet)


################################################################################
# 6. Microbiome profiling

if not readsPresent("6", unmapped_file):
	ARGS.microbiome = False

if ARGS.microbiome:
    write2Log("6. Microbiome profiling...", LOGFNS["cmdLogfile"], True)
    write2Log("6. Microbiome profiling...", LOGFNS["gLogfile"], ARGS.quiet)
	
    
    # 6b. bacteria
    os.chdir(DIRS["bacteria"])
    write2Log("BWA was used to map the reads onto the	bacterial reference genomes.", LOGFNS["logBacteria"], True)
    write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides aligned are considered microbial reads and are saved into " + LOGFNS["bacteriaFileFiltered"], LOGFNS["logBacteria"], True)
    write2Log("---------------", LOGFNS["logBacteria"], True)
    
    
    
    command=read_commands()
    cmd_bacteria=command[3]
    
    cmd = cmd_bacteria + " " + unmapped_file + "|"+ CD + "/tools/samtools view -SF4 - > bacteria.sam 2>>" +LOGFNS["logBacteria"]

    print (cmd)


    write2Log(cmd, LOGFNS["cmdLogfile"], True)
    if ARGS.qsub or ARGS.qsubArray:
        write2Log(cmd, RUNFNS["runBacteriaFile"], True)
        write2Log("echo \"done!\" >" + DIRS["bacteria"] + "/" + BASENAME +"_bacteria.done", RUNFNS["runBacteriaFile"], True)
        if ARGS.qsub:
            qsub("6b_bacteria", RUNFNS["runBacteriaFile"])
    else:
        nReads["bacteria"] = step_6b(unmapped_file, readLength, cmd,ARGS.max,ARGS.pe)
        write2Log("--identified " + str(nReads["bacteria"]) + " reads " + "mapped to bacterial genomes", LOGFNS["gLogfile"], ARGS.quiet)
        if not ARGS.nonReductive:
            clean(unmapped_file)
            unmapped_file = INTFNS["afterBacteriaFasta"]



    
    # 6b. virus
    os.chdir(DIRS["virus"])
    write2Log("BWA was used to map the reads onto the	viral reference genomes.", LOGFNS["logVirus"], True)
    write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides aligned are considered microbial reads and are saved into " + LOGFNS["virusFileFiltered"], LOGFNS["logVirus"], True)
    write2Log("---------------", LOGFNS["logVirus"], True)

    command=read_commands()
    cmd_virus=command[4]

    cmd =   cmd_virus  + " "+ unmapped_file + "|"+ CD + "/tools/samtools view -SF4 - > virus.sam 2>>" +LOGFNS["logVirus"]

    
    write2Log(cmd, LOGFNS["cmdLogfile"], True)
    if ARGS.qsub or ARGS.qsubArray:
        write2Log(cmd, RUNFNS["runBacteriaFile"], True)
        write2Log("echo \"done!\" >" + DIRS["bacteria"] + "/" + BASENAME +"_bacteria.done", RUNFNS["runBacteriaFile"], True)
        if ARGS.qsub:
            qsub("6b_bacteria", RUNFNS["runVirusFile"])
    else:
        nReads["virus"] = step_6c(unmapped_file, readLength, cmd,ARGS.max,ARGS.pe)
        write2Log("--identified " + str(nReads["virus"]) + " reads " + "mapped to viral genomes", LOGFNS["gLogfile"], ARGS.quiet)
        if not ARGS.nonReductive:
            clean(unmapped_file)
            unmapped_file = INTFNS["afterVirusFasta"]





	# 6d. eukaryotic pathogens	
    os.chdir(DIRS["eupathdb"])
    dbList = ["amoebadb", "cryptodb", "giardiadb", "microsporidiadb", "piroplasmadb","plasmodb", "toxodb", "trichdb", "tritrypdb"]




    for db in dbList:
        print (db)
        afterFasta = DIRS["eupathdb"] + BASENAME + "_after_" + db + ".fasta"
        eupathdbFile = DIRS["eupathdb"] + "/" + db +".sam"
        
        eupathdbFileFiltered = DIRS["eupathdb"] + BASENAME + "_" + db +"Filtered.csv"
        runEupathdbFile = DIRS["eupathdb"] + "/run_" + BASENAME + "_" + db +".sh"
        

        write2Log("BWA was used to map the reads onto the	reference " +"genomes of single-cell eukaryotes (http://eupathdb.org/eupathdb/).",LOGFNS["logEukaryotes"], True)
        write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides " +"aligned are considered microbial reads and are saved into " +afterFasta, LOGFNS["logEukaryotes"], True)
        write2Log("-------------", LOGFNS["logEukaryotes"], True)
        cmd = CD + "/prerequisite_software/bwa/bwa mem  " + CD + "/" + DB_FOLDER +"/single_cell_eukaryotes/EuPathDB_Merged_" +db + "_ConcatContigs.fasta "+ unmapped_file + "|"+ CD + "/tools/samtools view -SF4 - > "+db+".sam 2>>" +LOGFNS["logEukaryotes"]
        print (cmd)

        write2Log(cmd, LOGFNS["cmdLogfile"], True)
        if ARGS.qsub or ARGS.qsubArray:
            write2Log(cmd, runEupathdbFile, True)
            write2Log("echo \"done!\"  >" + DIRS["eupathdb"] +"/" + DIRS["eupathdb"] + ".done", runEupathdbFile, True)
            if ARGS.qsub:
                qsub("6d_" + db, runEupathdbFile)
        else:
            
            
            
            nEupathdbReads = step_6d(unmapped_file, readLength, cmd, eupathdbFile, eupathdbFileFiltered, afterFasta,db,ARGS.max,ARGS.pe)
            
            
            print (db,nEupathdbReads)
            
            write2Log("--identified " + str(nEupathdbReads) + " reads " +"mapped to " + db + " genomes", LOGFNS["gLogfile"], ARGS.quiet)
            nReads["ep"] += nEupathdbReads
            if not ARGS.nonReductive:
                clean(unmapped_file)
                unmapped_file = afterFasta
        


    if not (ARGS.qsub or ARGS.qsubArray):
        os.rename(DIRS["eupathdb"] + BASENAME + "_after_tritrypdb.fasta",INTFNS["unaccountedReadsFasta"])
        write2Log("In total: " + str(nReads["bacteria"] + nReads["virus"] +nReads["ep"]) + " reads mapped to microbial genomes",LOGFNS["gLogfile"], ARGS.quiet)
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


set_lost_human=set()
set_lost_repeat=set()
set_NCL=set()
set_immune=set()
set_bacteria=set()
set_virus=set()
set_eupathdb=set()




if os.path.isfile(DIRS["lostReads"]+'human_reads_SE.txt'):
    set_lost_human=set(line.strip() for line in open(DIRS["lostReads"]+'human_reads_SE.txt'))
if os.path.isfile(DIRS["lostReads"]+'repeat_reads_SE.txt'):
    set_lost_repeat=set(line.strip() for line in open(DIRS["lostRepeat"]+'repeat_reads_SE.txt'))
if os.path.isfile(DIRS["lostReads"]+'NLC_reads_SE.txt'):
    set_NCL=set(line.strip() for line in open(DIRS["NCL"]+'NLC_reads_SE.txt'))
if os.path.isfile(DIRS["lostReads"]+'immune_reads_SE.txt'):
    set_immune=set(line.strip() for line in open(DIRS["antibody"]+'immune_reads_SE.txt'))
if os.path.isfile(DIRS["lostReads"]+'bacteria_reads_SE.txt'):
    set_bacteria=set(line.strip() for line in open(DIRS["bacteria"]+'bacteria_reads_SE.txt'))
    set_virus=set(line.strip() for line in open(DIRS["virus"]+'viral_reads_SE.txt'))
    set_eupathdb=set(line.strip() for line in open(DIRS["eupathdb"]+'eukaryotic_reads_SE.txt'))



print (len(set_lost_human),len(set_lost_repeat),len(set_NCL),len(set_immune),len(set_bacteria),len(set_virus),len(set_eupathdb))
print (nReads["lost"]-len(set_lost_human))
print (nReads["repeat"]-len(set_lost_repeat))
print (nReads["NCL"]-len(set_NCL))
print (nReads["immune"]-len(set_immune))
print (nReads["bacteria"]-len(set_bacteria))
print (nReads["virus"]-len(set_virus))
print (nReads["ep"]-len(set_eupathdb))


write2Log("""The list of the tools used by ROP is provided below.
************
**We have used in-house script to filter out low quality reads"
**We have used Megablast (BLAST+ version 2.2.30, with the following parameters: task = megablast, use_index = true; -outfmt 6 ;-evalue 1e-05; perc_identity = 100) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to filter out reads mapped to rRNA repeat sequence
**We have used Bowtie2 (version 2.0.5, with the following parameters: -k 1; -p 8; -f) downloaded from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml to identify lost reads mapped to reference transcriptome and genome
**We have used Megablast (BLAST+ version 2.2.30, with the following options: task=megablast, use_index=true, -outfmt 6 -evalue 1e-05, perc_identity	= 90) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to identify lost repeat reads mapped to database of repeat sequences (RepBase 20.07)
**We have used CIRI (version 1.2 with the following parameters: -S -I ) downloaded from https://sourceforge.net/projects/ciri/ to identify reads from circRNAs
**We have used ImReP (version 0.3) to identify immune reads spanning BCR/TCR receptor gene rearrangement in the variable domain (V(D)J recombinations)
**We have used Megablast (BLAST+ version 2.2.30, with the following parameters: task=megablast, use_index=true; -outfmt 6 ;-evalue 1e-05) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to identify microbial reads mapped onto the microbial genomes (bacteria, viruses, and eukaryotic pathogens)
************
For more information about the paramemers and databases used by ROP, please see the preprint: Dumpster diving in RNA-sequencing to find the source of every last read http://biorxiv.org/content/early/2016/05/13/053041"
********************""", LOGFNS["toolsLogfile"], True)

write2Log("Important: ROP relies on several open source tools that were " +\
  "developed by other groups. These components are (c) their respective " +\
  "developers and are redistributed with ROP to provide ease-of-use. The " +\
  "list of the tools used by ROP and the parameters/reference databases are " +\
  "provided here: " + LOGFNS["toolsLogfile"], LOGFNS["gLogfile"], ARGS.quiet)
