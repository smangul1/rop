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
		print ("ERROR: Invalid input. " + message)

class SubprocessError(Exception):
	def __init__(self, message=""):
		print ("ERROR: Subprocess crashed. " + message)

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
		print (message)

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




def set2file(set,fileName):
    file=open(fileName,"a")
    
    
    set=sorted(set)
    for s in set:
        file.write(str(s))
        file.write("\n")

def pe2se(lostReads,fileName):
    file=open(fileName,"w")
    lostReads_se=set()
    for s in lostReads:
        lostReads_se.add(s.split("---")[0])

    for s in lostReads_se:
        file.write(s)
        file.write("\n")
        
    file.close()

################################################################################
# .bam

def bam2fasta(cd, bam_name, fasta_name):
	
	
	import pysam
	samfile = pysam.AlignmentFile(bam_name, "rb")
	message = "Convert bam to fasta"
	write2Log(message, LOGFNS["gLogfile"], ARGS.quiet)
	cmd = cd + "/tools/bamtools convert -in " + bam_name + " -format fasta >" +fasta_name
	write2Log(cmd, LOGFNS["cmdLogfile"], "False")
	if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()

def bam2fastq(cd, bam_name, fastq_name):
    message = "Extract unmapped reads and convert to fastq"
    write2Log(message, LOGFNS["gLogfile"], ARGS.quiet)
    
    cmd=cd + "/tools/samtools view -f 0x4 -bh " + bam_name + " | " + cd+ "/tools/samtools bam2fq ->"+fastq_name
    print cmd
    write2Log(cmd, LOGFNS["cmdLogfile"], "False")
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()



#Customize the tools used by ROP

def read_commands():
    f=open(CD+'/rop.commands.txt')
    command=[]
    reader=csv.reader(f)
    for line in reader:
        command.append(line[0])

    return command




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



def extract_from_read(line,readLength):

    alignment=line[12].split(':')[2]  #MD:Z:53T22
    ed=int(line[11].split(':')[2])  #NM:i:0
    countA=alignment.count('A')
    countC=alignment.count('C')
    countT=alignment.count('T')
    countG=alignment.count('G')
    alignment2=alignment.replace('A',',').replace('C',',').replace('T',',').replace('G',',').replace('^',',').split(',')
    alignment3=[]
            
    for a in alignment2:
        if a!="":
            alignment3.append(int(a))
        
    alignmentLength=sum(alignment3)+countA+countC+countG+countT

    clipped=readLength-alignmentLength
    ed=readLength-alignmentLength+int(line[11].split(':')[2])  #NM:i:0

    return (ed,alignmentLength,clipped)

def nMicrobialReads(inFile_name, readLength, outFile_name,flag):
    
    readsMicrobiome = set()
    with open(inFile_name, "r") as inFile:
        with open(outFile_name, "w") as outFile:
            for line in csv.reader(inFile, delimiter='\t'):
                read = line[0]
                
                
                
                alignment=line[12].split(':')[2]  #MD:Z:53T22

                ed=int(line[11].split(':')[2])  #NM:i:0
                
                
                
                
                countA=alignment.count('A')
                countC=alignment.count('C')
                countT=alignment.count('T')
                countG=alignment.count('G')
                
                alignment2=alignment.replace('A',',').replace('C',',').replace('T',',').replace('G',',').replace('^',',').split(',')
                
                
                alignment3=[]
                
                for a in alignment2:
                    if a!="":
                        alignment3.append(int(a))
            
                alignmentLength=sum(alignment3)+countA+countC+countG+countT
                
                
                identity=readLength-ed
                
                
                
                if flag!=True:
                    
                    if  alignmentLength >= 0.8*readLength and identity >= 0.9*readLength:
                        readsMicrobiome.add(read)
                        outFile.write("\t".join(line) + "\n")
                else:
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
    fasta_sequences = SeqIO.parse(open(inFile_name),'fastq')
    for seq in fasta_sequences:
        reads.add(seq.id)
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
    valid = "ACTGN"
    nLowQReads = 0
    with open(INTFNS["lowQFileFasta"], "w") as lowQFileFasta:
        for record in SeqIO.parse(unmapped_file, "fastq"):
            # assumes the same length, will not work for Ion Torrent or Pac Bio
            j = record.letter_annotations["phred_quality"]
            prc = len([i for i in j if i >= 20])/float(len(j))
            if prc > 0.75 and all(i in valid for i in record.seq):
                lowQFileFasta.write(str(">" + record.name.replace("/","---")) + "\n")
                lowQFileFasta.write(str(record.seq) + "\n")
                
            else:
                lowQFileFasta.write(str(">lowQuality_" + record.name.replace("/","---")) + "\n")
                lowQFileFasta.write(str(record.seq) + "\n")
                nLowQReads += 1
				
    return nLowQReads



def step_1c(unmapped_file, readLength):
    cmd = CD + "/tools/blastn -task megablast -index_name " + CD + "/" +DB_FOLDER + "/rRNA/rRNA -use_index true -query " + unmapped_file +" -db " + CD + "/" + DB_FOLDER + "/rRNA/rRNA -outfmt 6 -evalue 1e-05 >" +INTFNS["rRNAFile"] + " 2>log_megablast_rRNA.log"
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



def step_2(unmapped_file,flag,flag_PE,readLength):

    if flag:
        ed_human=10
    else:
        ed_human=6





    command=read_commands()
    bwa_WG=command[0]
    bwa_TR=command[1]


    cmdGenome = bwa_WG + " "+unmapped_file + " 2>>" + LOGFNS["log_bowtieWG"] + " | " + CD + "/tools/samtools view -SF4 -bh - >" +INTFNS["gBamFile"]
    cmdTranscriptome = bwa_TR + " "+unmapped_file + " 2>>" + LOGFNS["log_bowtieTR"] + " | " + CD + "/tools/samtools view -SF4 -bh -  >" + INTFNS["tBamFile"]





    write2Log(cmdGenome, LOGFNS["cmdLogfile"], True)
    write2Log(cmdTranscriptome, LOGFNS["cmdLogfile"], True)
    proc1 = subprocess.Popen([cmdTranscriptome], shell=True)
    proc2 = subprocess.Popen([cmdGenome], shell=True)
    if proc1.wait() or proc2.wait(): raise SubprocessError()

    cmd=CD+"/tools/samtools view " + INTFNS["gBamFile"] + ">" + INTFNS["gSamFile"]
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    
    cmd=CD+"/tools/samtools view " + INTFNS["tBamFile"] + ">" + INTFNS["tSamFile"]
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()


    lostReads = set()
    lostReads0 = set()
    lostReads1 = set()
    lostReads2 = set()
    with open(INTFNS["tSamFile"], "r") as f:
        for line in csv.reader(f, delimiter='\t'):
            flag_OK=True # flag: 0 consider read; 1 don't consider
            cigar=line[5]
            
            (ed,alignmentLength,clipped)=extract_from_read(line,readLength)
            
         


            if  ed<=ed_human:


                lostReads.add(line[0])
                if ed == 0:
                    lostReads0.add(line[0])
                elif ed == 1:
                    lostReads1.add(line[0])
                elif ed == 2:
                    lostReads2.add(line[0])
                        
                        
    with open(INTFNS["gSamFile"], "r") as f:
        for line in csv.reader(f,delimiter='\t'):
            
            (ed,alignmentLength,clipped)=extract_from_read(line,readLength)
            identity=readLength-ed
            
            
            
            
            
            if  ed<=ed_human:

                lostReads.add(line[0])
                if ed == 0:
                    lostReads0.add(line[0])
                elif ed == 1:
                    lostReads1.add(line[0])
                elif ed == 2:
                    lostReads2.add(line[0])
    
    
    write2Log("Unmapped reads mapped to genome and/or transcriptome (using " + "BWA, edit distance<6) are categorized as lost reads and are excluded from the further analysis. This includes : " + str(len(lostReads0)) + " reads with 0 mismatches, " + str(len(lostReads1)) + " reads with 1 mismatch, and" + str(len(lostReads2)) + " reads with 2 mismatches.", LOGFNS["logHuman"],True)
    write2Log("Complete list of lost reads is available from sam files: " + INTFNS["gBamFile"] + ", " + INTFNS["tBamFile"], LOGFNS["logHuman"], True)
    excludeReadsFromFasta(unmapped_file, lostReads,INTFNS["afterlostReadsFasta"])

    if flag_PE:
        pe2se(lostReads,"lost_human_reads_SE.txt")
    set2file(lostReads,"lost_human_reads.txt")

    os.remove(INTFNS["gSamFile"])
    os.remove(INTFNS["tSamFile"])





    return len(lostReads), len(lostReads0), len(lostReads1), len(lostReads2)

#-----------------------------------------------------------------------------------------------------------------------
def step_3(unmapped_file, readLength, cmd,flag,flag_PE):
    write2Log(cmd, LOGFNS["cmdLogfile"], True)
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    lostRepeatReads = set()
    with open(INTFNS["repeatFile"], "r") as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            element = line[0]
            identity = float(line[2])
            alignmentLength = float(line[3])
            eValue = float(line[10])
            


            if (eValue < 1e-05 and alignmentLength >= 0.8*readLength and identity >= 0.9*readLength):
                lostRepeatReads.add(element)
                
                
    write2Log("Lost repeat reads are mapped to the repeat sequences " +"(using megablast)", LOGFNS["logLostRepeat"], True)
    write2Log("Complete list of lost repeat reads is available from tsv " +"file: " + INTFNS["repeatFile"], LOGFNS["logLostRepeat"], True)
    excludeReadsFromFasta(unmapped_file, lostRepeatReads, INTFNS["afterlostRepeatFasta"])
      
      
    if flag_PE:
        pe2se(lostRepeatReads,"lost_repeat_reads_SE.txt")
    set2file(lostRepeatReads,"lost_repeat_reads.txt")
      
    return len(lostRepeatReads)




#-----------------------------------------------------------------------------------------------------------------------
def step_4(unmapped_file, cmd,flag_PE):
    NCL_reads = set()
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    NCL_reads = nCircularReads("accepted_hits.fastq")
    excludeReadsFromFasta(unmapped_file, NCL_reads, INTFNS["afterNCLFasta"])
    
    
    if flag_PE:
        pe2se(NCL_reads,"NLC_reads_SE.txt")
    set2file(NCL_reads,"NCL_reads.txt")
    
    return len(NCL_reads)

def step_5(unmapped_file, cmd,flag_PE):
    proc = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE)
    (output, err) = proc.communicate()
    if proc.wait(): raise SubprocessError()
    write2Log(output.strip(), LOGFNS["gLogfile"], ARGS.quiet)
    
    immuneReads=set()
    
    
    
    imrep_filename="full_cdr3_"+os.path.basename(unmapped_file).replace(".fasta","")+".txt"
    
    
    file=open(imrep_filename)
    reader=csv.reader(file, delimiter='\t')
    next(reader, None)
    for line in reader:
        immuneReads.add(line[0])
    
    
    
    excludeReadsFromFasta(unmapped_file, immuneReads,INTFNS["afterImmuneFasta"])

    if flag_PE:
        pe2se(immuneReads,"immune_reads_SE.txt")
    set2file(immuneReads,"immune_reads.txt")


    return len(immuneReads)


def step_6b(unmapped_file, readLength, cmd,flag,flag_PE):
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    bacteriaReads =nMicrobialReads("bacteria.sam", readLength,LOGFNS["bacteriaFileFiltered"],flag)
    excludeReadsFromFasta(unmapped_file, bacteriaReads, INTFNS["afterBacteriaFasta"])
    if flag_PE:
        pe2se(bacteriaReads,"bacteria_reads_SE.txt")
    set2file(bacteriaReads,"immune_reads.txt")

    return len(bacteriaReads)


#virus
def step_6c(unmapped_file, readLength, cmd,flag,flag_PE):
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    
    
    cmd=CD+"/tools/samtools view " + INTFNS["bam_viral"] + ">" + INTFNS["sam_viral"]
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    cmd=CD+"/tools/samtools view " + INTFNS["bam_viral_vipr"] + ">" + INTFNS["sam_viral_vipr"]
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    
    virusReads_NCBI = nMicrobialReads(INTFNS["sam_viral"], readLength,LOGFNS["virusFileFiltered"],flag)
    virusReads_VIPR = nMicrobialReads(INTFNS["sam_viral_vipr"], readLength,LOGFNS["virusFileFiltered"],flag)
    virusReads=virusReads_NCBI | virusReads_VIPR
    excludeReadsFromFasta(unmapped_file, virusReads, INTFNS["afterVirusFasta"])
    if flag_PE:
        pe2se(virusReads,"viral_reads_SE.txt")
    set2file(virusReads,"viral_reads.txt")

    os.remove(INTFNS["sam_viral"])
    os.remove(INTFNS["sam_viral_vipr"])

    return len(virusReads)

def step_6d(unmapped_file, readLength, cmd,flag,flag_PE):
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    
    cmd=CD+"/tools/samtools view " + INTFNS["bam_fungi"] + ">" + INTFNS["sam_fungi"]
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    
    fungiReads = nMicrobialReads(INTFNS["sam_fungi"], readLength,LOGFNS["fungiFileFiltered"],flag)
    excludeReadsFromFasta(unmapped_file, fungiReads, INTFNS["afterFungiFasta"])
    if flag_PE:
        pe2se(fungiReads,"fungi_reads_SE.txt")
    set2file(fungiReads,"fungi_reads.txt")

    os.remove(INTFNS["sam_fungi"])

    return len(fungiReads)


def step_6_protozoa(unmapped_file, readLength, cmd,flag,flag_PE):
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    
    cmd=CD+"/tools/samtools view " + INTFNS["bam_protozoa"] + ">" + INTFNS["sam_protozoa"]
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    
    protozoaReads = nMicrobialReads(INTFNS["sam_protozoa"], readLength,LOGFNS["protozoaFileFiltered"],flag)
    excludeReadsFromFasta(unmapped_file, protozoaReads, INTFNS["afterProtozoaFasta"])
    if flag_PE:
        pe2se(protozoaReads,"protozoa_reads_SE.txt")
    set2file(protozoaReads,"protozoa_reads.txt")

    os.remove(INTFNS["sam_protozoa"])
    return len(protozoaReads)





