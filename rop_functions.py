"""********************************************************************************
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

from rop_globals import *

################################################################################
### UTILITY
################################################################################

class InputError(Exception):
    def __init__(self, message=""):
        print "ERROR: Invalid input. " + message

class SubprocessError(Exception):
    def __init__(self, message=""):
        print "ERROR: Subprocess crashed. " + message

# deletes a file except under specified conditions
def clean(unmapped_file):
	if unmapped_file != ARGS.unmappedReads and not ARGS.dev:
		os.remove(unmapped_file)

def qsub(step_no, run_file):
	if ARGS.maui:
		cmdQsub = "qsub -d `pwd` -l walltime=10:00:00 -l " +\
		  "nodes=1:m16G:ppn=12 " + run_file
	else:
		cmdQsub = "qsub -cwd -V -N rop_" + step_no + " -l " +\
		  "h_data=16G,time=10:00:00 " + run_file
	if subprocess.Popen([cmdQsub], shell=True).wait(): raise SubprocessError()
	write2Log("Job for step " + str(step_no) + " has been submitted via qsub.",
	  LOGFNS["gLogfile"], ARGS.quiet)
	  
def readsPresent(step_no, unmapped_file):
	if not os.path.getsize(unmapped_file) > 1:
		write2Log("WARNING: No unmapped reads! Skipping step " + step_no + ".", 
		  LOGFNS["gLogfile"], ARGS.quiet)
		return False
	else:
		return True

	  
################################################################################
### I/O 
################################################################################
# Write to log or file

def write2Log(message, logFile_name, quiet):
	with open(logFile_name, "a") as logFile:
		logFile.write(message + '\n')
	if not quiet:
		print message

def write2File(content, file_name):
	with open(file_name, "w") as file:
		file.write(content)

################################################################################
# Exclude reads

def excludeReadsFromFasta(inFasta_name, reads, outFasta_name):
	fasta_seqs = SeqIO.parse(open(inFasta_name), "fasta")
	with open(outFasta_name, "w") as outFasta:
		for seq in fasta_seqs:
			if seq.name not in reads:
				SeqIO.write([seq], outFasta, "fasta")

def excludeReadsFromFastaGzip(inFasta_name, reads, outFasta_name):
	fasta_seqs = SeqIO.parse(open(inFasta_name), "fasta")
	with gzip.open(outFasta_name, "w") as outFasta:
		for seq in fasta_seqs:
			if seq.name not in reads:
				SeqIO.write([seq], outFasta, "fasta")

################################################################################
# .bam

def bam2fasta(cd, bam_name, fasta_name):
	message = "Convert bam to fasta"
	write2Log(message, LOGFNS["gLogfile"], ARGS.quiet)
	cmd = cd + "/tools/bamtools convert -in " + bam_name + " -format fasta >" +\
	  fasta_name
	write2Log(cmd, LOGFNS["cmdLogfile"], "False")
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()

def bam2fastq(cd, bam_name, fastq_name):
	message = "Convert bam to fastq"
	write2Log(message, LOGFNS["gLogfile"], ARGS.quiet)
	cmd = cd + "/tools/bamtools convert -in " + bam_name + " -format fastq >" +\
	  fastq_name
	write2Log(cmd, LOGFNS["cmdLogfile"], "False")
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()

################################################################################
# .gz

def write_gzip_into_readable(gz_name, outFile_name): 
	with gzip.open(gz_name, "r") as gz:
		with open(outFile_name, "w") as outFile:
			for line in gz:
				outFile.write(line)

def write_gzip(inFasta_name, outFasta_name):
	fasta_seqs = SeqIO.parse(open(inFasta_name), "fasta")
	with gzip.open(outFasta_name, "w") as outFasta:
		for seq in fasta_seqs:
			SeqIO.write([seq], outFasta, "fasta")


################################################################################
### RESULT COUNTING
################################################################################

def nReadsImmune(inFile_name):
	readsImmune = set()
	with open(inFile_name, "r") as inFile:
		for line in csv.reader(inFile, delimiter='\t'):
			read = line[1]
			eValue = float(line[12])
			if eValue < 1e-05:
				readsImmune.add(read)
	return readsImmune

def nMicrobialReads(inFile_name, readLength, outFile_name):
	readsMicrobiome = set()
	with open(inFile_name, "r") as inFile:
		with open(outFile_name, "w") as outFile:
			for line in csv.reader(inFile, delimiter='\t'):
				read = line[0]
				identity = float(line[2])
				alignmentLength = float(line[3])
				eValue = float(line[10])
				if (eValue < 1e-05 
				  and alignmentLength >= 0.8*readLength 
				  and identity >= 0.9*readLength):
					readsMicrobiome.add(read)
					outFile.write("\t".join(line) + "\n")
	return readsMicrobiome

def nReadsMetaphlan(inFile_name):
	readsMetaphlan = set()
	with open(inFile_name, "r") as inFile:
		for line in csv.reader(inFile, delimiter='\t'):
			read = line[0]
			readsMetaphlan.add(read)
	return readsMetaphlan

def nCircularReads(inFile_name):
	reads=set()
	with open(inFile_name, "r") as inFile:
		inFile.next()
		for line in csv.reader(inFile, delimiter='\t'):
			for kv in line[10].split(","):
				if kv:
					reads.add(kv)
	return reads

	
################################################################################
### PROGRAM FUNCTIONS
################################################################################

# input: name of unmapped reads file (bam/fastx/gz)
# output: number of reads, read length, name of fastx file
def prepare_for_analysis(unmapped_file):	
	# conversions and checks
	to_fasta = ARGS.skipPreliminary or ARGS.skipQC or ARGS.skipLowq
	if ARGS.b:  # bam input
		if to_fasta:  # convert to fasta
			bam2fasta(CD, unmapped_file, INTFNS["lowQFileFasta"])
			unmapped_file = INTFNS["lowQFileFasta"]
		else:  # convert to fastq
			bam2fastq(CD, unmapped_file, INTFNS["unmappedFastq"])
			unmapped_file = INTFNS["unmappedFastq"]
	elif ARGS.gzip:  # gz input
		if to_fasta:  # convert to fasta
			INTFNS["lowQFileFasta"] = ARGS.dir + "/" + BASENAME + ".fasta"
			write_gzip_into_readable(unmapped_file, INTFNS["lowQFileFasta"])
			unmapped_file = INTFNS["lowQFileFasta"]
		else:  # convert to fastq
			INTFNS["unmappedFastq"] = ARGS.dir + "/" + BASENAME + ".fastq"
			write_gzip_into_readable(unmapped_file, INTFNS["unmappedFastq"])
			unmapped_file = INTFNS["unmappedFastq"]
	elif to_fasta:  # fasta input
		filename, file_extension = os.path.splitext(unmapped_file)
		if not (file_extension == ".fa" or file_extension == ".fasta"):
			raise InputError(".fasta file required.")
	else:  # fastq input
		filename, file_extension = os.path.splitext(unmapped_file)
		if not (file_extension == ".fq" or file_extension == ".fastq"):
			raise InputError(".fastq file required.")
	
	# detect number of reads and read length
	n = 0
	readLength = 0
	if to_fasta:
		with open(unmapped_file, "rU") as fastafile:
			for record in SeqIO.parse(fastafile, "fasta"):
				readLength=len(record)  
				  # assumes the same length, will not work for Ion Torrent or Pac Bio
				n += 1
	else:  # process the fastq
		with open(unmapped_file, "rU") as fastqfile:
			for record in SeqIO.parse(fastqfile, "fastq"):
				readLength=len(record) 
				  # assumes the same length, will not work for Ion Torrent or Pac Bio
				n += 1
	return n, readLength, unmapped_file
	
def step_1a(unmapped_file, n):
	valid = "ACTG"
	nAfterLowQReads = 0
	with open(INTFNS["lowQFileFasta"], "w") as lowQFileFasta:
		for record in SeqIO.parse(unmapped_file, "fastq"):
			  # assumes the same length, will not work for Ion Torrent or Pac Bio
			j = record.letter_annotations["phred_quality"]
			prc = len([i for i in j if i >= 20])/float(len(j))
			if prc > 0.75 and all(i in valid for i in record.seq):
				lowQFileFasta.write(str(">" + record.name) + "\n")
				lowQFileFasta.write(str(record.seq) + "\n")
				nAfterLowQReads += 1
	return n - nAfterLowQReads

def step_1b(unmapped_file):
	cmd = "export PATH=$PATH:" + CD + "/tools/seqclean-x86_64/bin"
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
	cmd = CD + "/tools/seqclean-x86_64/seqclean " + unmapped_file +\
	  " -l 50 -M -o " + INTFNS["lowQCFile"] + " 2>>" + LOGFNS["logQC"]
	write2Log(cmd, LOGFNS["cmdLogfile"], True)
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
	cmd = "rm -rf " + DIRS["QC"] + "/cleaning_1/ >temp 2>temp; " +\
	  "rm -f " + DIRS["QC"] + "/*.cln >temp 2>temp; " +\
	  "rm -f " + DIRS["QC"] + "/*.cidx >temp 2>temp; " +\
	  "rm -f " + DIRS["QC"] + "/*.sort >temp 2>temp" 
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
	cmd = "grep trashed " + LOGFNS["logQC"] + " | awk -F \":\" '{print $2}'"
	proc = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE)
	(nLowCReadsTemp, err) = proc.communicate()
	if proc.wait(): raise SubprocessError()
	return int(nLowCReadsTemp.strip())

def step_1c(unmapped_file, readLength):
	cmd = CD + "/tools/blastn -task megablast -index_name " + CD + "/" +\
	  DB_FOLDER + "/rRNA/rRNA -use_index true -query " + unmapped_file +\
	  " -db " + CD + "/" + DB_FOLDER + "/rRNA/rRNA -outfmt 6 -evalue 1e-05 >" +\
	  INTFNS["rRNAFile"]
	write2Log(cmd, LOGFNS["cmdLogfile"], True)
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
	n_rRNATotal = 0
	rRNAReads = set()
	with open(INTFNS["rRNAFile"], "r") as f:
		reader = csv.reader(f, delimiter='\t')
		for line in reader:
			n_rRNATotal += 1
			element = line[0]
			identity = float(line[2])
			alignmentLength = float(line[3])
			eValue = float(line[10])
			if (eValue < 1e-05 and alignmentLength == readLength and 
			  identity >= 0.94 * readLength):
				rRNAReads.add(element)
	excludeReadsFromFasta(unmapped_file, rRNAReads, INTFNS["afterrRNAFasta"])
	return len(rRNAReads), n_rRNATotal

def step_2(unmapped_file):
	cmdGenome = CD + "/tools/bowtie2 -k 1 -p 8 -f -x " + CD + "/" + DB_FOLDER +\
	  "/bowtie2Index/genome -U " + unmapped_file + " 2>>" +\
	  LOGFNS["log_bowtieWG"] + " | " + CD + "/tools/samtools view -SF4 - >" +\
	  INTFNS["gBamFile"]
	cmdTranscriptome = CD + "/tools/bowtie2 -k 1 -f -p 8 -x " + CD + "/" +\
	  DB_FOLDER + "/bowtie2Index/hg19KnownGene.exon_polya200 -U " +\
	  unmapped_file + " 2>>" + LOGFNS["log_bowtieTR"] + " | " + CD +\
	  "/tools/samtools view -SF4 -  >" + INTFNS["tBamFile"]
	write2Log(cmdGenome, LOGFNS["cmdLogfile"], True)
	write2Log(cmdTranscriptome, LOGFNS["cmdLogfile"], True)
	proc1 = subprocess.Popen([cmdTranscriptome], shell=True)
	proc2 = subprocess.Popen([cmdGenome], shell=True)
	if proc1.wait() or proc2.wait(): raise SubprocessError()
	lostReads = set()
	lostReads0 = set()
	lostReads1 = set()
	lostReads2 = set()
	with open(INTFNS["tBamFile"], "r") as f:
		for line in csv.reader(f, delimiter='\t'):
			mismatch = int(line[16].split(':')[2])
			if mismatch < 3:
				lostReads.add(line[0])
			if mismatch == 0:
				lostReads0.add(line[0])
			elif mismatch == 1:
				lostReads1.add(line[0])
			elif mismatch == 2:
				lostReads2.add(line[0])
	with open(INTFNS["gBamFile"], "r") as f:
		for line in csv.reader(f,delimiter='\t'):
				mismatch=int(line[16].split(':')[2])
				if mismatch < 3:
					lostReads.add(line[0])
				if mismatch == 0:
					lostReads0.add(line[0])
				elif mismatch == 1:
					lostReads1.add(line[0])
				elif mismatch == 2:
					lostReads2.add(line[0])
	write2Log("Unmapped reads mapped to genome and/or transcriptome (using " +\
	  "bowtie2) are categorized as lost reads and are excluded from the " +\
	  "further analysis. This includes: " +\
	  str(len(lostReads0)) + " reads with 0 mismatches, " +\
	  str(len(lostReads1)) + " reads with 1 mismatch, and" +\
	  str(len(lostReads2)) + " reads with 2 mismatches.", LOGFNS["logHuman"], 
	    True)
	write2Log("Complete list of lost reads is available from sam files: " +\
	  INTFNS["gBamFile"] + ", " + INTFNS["tBamFile"], LOGFNS["logHuman"], True)
	excludeReadsFromFasta(unmapped_file, lostReads, 
	  INTFNS["afterlostReadsFasta"])
	return len(lostReads), len(lostReads0), len(lostReads1), len(lostReads2)
	
def step_3(unmapped_file, readLength, cmd):
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
	write2Log(cmd, LOGFNS["cmdLogfile"], True)
	lostRepeatReads = set()
	with open(INTFNS["repeatFile"], "r") as f:
		reader = csv.reader(f, delimiter='\t')
		for line in reader:
			element = line[0]
			identity = float(line[2])
			alignmentLength = float(line[3])
			eValue = float(line[10])
			if (eValue < 1e-05 and alignmentLength >= 0.8*readLength and 
			  identity >= 0.9*readLength):
				lostRepeatReads.add(element)
	write2Log("Lost repeat reads are mapped to the repeat sequences " +\
	  "(using megablast)", LOGFNS["logLostRepeat"], True)
	write2Log("Complete list of lost repeat reads is available from tsv " +\
	  "file: " + INTFNS["repeatFile"], LOGFNS["logLostRepeat"], True)
	excludeReadsFromFasta(unmapped_file, lostRepeatReads, 
	  INTFNS["afterlostRepeatFasta"])
	return len(lostRepeatReads)

def step_4(unmapped_file, cmd):
	NCL_reads = set()
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
	NCL_reads = nCircularReads(INTFNS["after_NCL_CIRI_file_prefix"])
	excludeReadsFromFasta(unmapped_file, NCL_reads, INTFNS["afterNCLFasta"])
	return len(NCL_reads)
	
def step_6b(unmapped_file, readLength, cmd):
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
	bacteriaReads = nMicrobialReads(LOGFNS["bacteriaFile"], readLength, 
	  LOGFNS["bacteriaFileFiltered"])
	excludeReadsFromFasta(unmapped_file, bacteriaReads, 
	  INTFNS["afterBacteriaFasta"])
	return len(bacteriaReads)

def step_6c(unmapped_file, readLength, cmd):
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
	virusReads = nMicrobialReads(LOGFNS["virusFile"], readLength,
	  LOGFNS["virusFileFiltered"])
	excludeReadsFromFasta(unmapped_file, virusReads, INTFNS["afterVirusFasta"])
	return len(virusReads)

def step_6d(unmapped_file, readLength, cmd, eupathdbFile, eupathdbFileFiltered,\
  afterFasta):
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
	eupathdbReads = set()
	eupathdbReads = nMicrobialReads(eupathdbFile, readLength, 
	  eupathdbFileFiltered)
	excludeReadsFromFasta(unmapped_file, eupathdbReads, afterFasta)
	return len(eupathdbReads)
