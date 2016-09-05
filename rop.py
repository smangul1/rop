"""
ROP is a computational protocol aimed to discover the source of all reads, 
originated from complex RNA molecules, recombinant antibodies and microbial 
communities.

Written by Serghei Mangul (smangul@ucla.edu), Harry Taegyun Yang
(harry2416@gmail.com), Kevin Hsieh(kevin.hsieh@ucla.edu), Linus Chen (u6.30cl@gmail.com), University of California, Los Angeles (UCLA). (c) 2016.

Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see: https://sergheimangul.wordpress.com/rop/
ROP Tutorial: https://github.com/smangul1/rop/wiki 
"""

import sys, os  # system imports
import argparse, subprocess  # utilities
import csv, gzip  # file I/O

cd = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cd + "/tools/biopython/biopython-1.66/")
from Bio import SeqIO  # module needed for sequence input
sys.path.append(cd + "/tools/pysam-master/")
#import pysam (future use)


################################################################################
### I/O FUNCTIONS 
################################################################################
# Exclude reads from FASTA

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
# BAM to FASTA/FASTQ

def bam2fasta(cd, bam_name, fasta_name):
    message = "Convert bam to fasta"
    write2Log(message, logfns["gLogfile"], args.quiet)
    cmd = cd + "/tools/bamtools convert -in " + bam_name + " -format fasta >" + fasta_name
    write2Log(cmd, logfns["cmdLogfile"], "False")
    
    if subprocess.Popen([cmd], shell=True).wait():
        sys.exit(2)

def bam2fastq(cd, bam_name, fastq_name):
    message = "Convert bam to fastq"
    write2Log(message, logfns["gLogfile"], args.quiet)
    cmd = cd + "/tools/bamtools convert -in " + bam_name + " -format fastq >" + fastq_name
    write2Log(cmd, logfns["cmdLogfile"], "False")
    
    if subprocess.Popen([cmd], shell=True).wait():
        sys.exit(2)

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
# GZ functions

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
# nReads functions

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
### ARGUMENTS 
################################################################################

ap = argparse.ArgumentParser("python rop.py")

necessary_arguments = ap.add_argument_group("Necessary Inputs")
necessary_arguments.add_argument("unmappedReads", 
    help="unmapped reads in .fastq format")
necessary_arguments.add_argument("dir", 
    help="directory to save results of the analysis")

job_option_arguments = ap.add_argument_group('Job Options')
job_option_arguments.add_argument("--qsub", 
    help="submit qsub jobs on Hoffman2 (UCLA) cluster. If planning to use on your cluster, contact smangul@ucla.edu", 
    action="store_true")
job_option_arguments.add_argument("--qsubArray", 
    help="prepare qsub scripts to be run later using job array. Working on Hoffman2 (UCLA) cluster. If planning to use on your cluster contact smangul@ucla.edu", 
    action="store_true")
job_option_arguments.add_argument("--maui", 
    help="use this option together with --qsub to submit jobs via Maui scheduler. Maui is a job scheduler developped by Adaptive Computing. More details are here: https://wiki.calculquebec.ca/w/Maui/en", 
    action="store_true")
job_option_arguments.add_argument("--organism", 
    help="run ROP for a specified organism (mouse) instead of human", 
    action="store", default="human")

input_option_arguments = ap.add_argument_group('Input Options')
input_option_arguments.add_argument("--b", "-b", 
    help="unmapped reads in .bam format instead of .fastq", 
    action="store_true")
input_option_arguments.add_argument("--gzip", "-z", 
    help="unmapped reads in .gz format", 
    action="store_true")
input_option_arguments.add_argument("--skipLowq", 
    help="skip filtering low quality reads. The input reads need to be in .bam format", 
    action="store_true")
input_option_arguments.add_argument("--skipQC", 
    help="skip entire QC step (filtering low-quality, low-complexity and rRNA reads). The input reads need to be in .fasta format", 
    action="store_true")
input_option_arguments.add_argument("--skipPreliminary", "-s", 
    help="skip the preliminary steps including (1) QC and (2) remapping to references (lost reads). The input reads need to be in .fasta format", 
    action="store_true")

run_only_options = ap.add_argument_group('Run Options - select analysis (multiple selections possible)')
run_only_options.add_argument("--repeat", 
    help="Run lost repeat profiling ONLY", 
    action="store_true")
run_only_options.add_argument("--immune", 
    help="Run BCR/TCR profiling ONLY", 
    action="store_true")
run_only_options.add_argument("--metaphlan", 
    help="Run Metaphlan2 ONLY (to obtain taxonomic profile of microbial communities)", 
    action="store_true")
run_only_options.add_argument("--circRNA", 
    help="Run circular RNA profiling ONLY", 
    action="store_true")
run_only_options.add_argument("--microbiome", 
    help="Run microbime profiling ONLY", 
    action="store_true")

misc_option_arguments = ap.add_argument_group('Miscellenous Options')
misc_option_arguments.add_argument("--outGz", 
    help="beta mode. Intermediate .fasta files are stored as .fasta.gz", 
    action="store_true")
misc_option_arguments.add_argument("--rezip", 
    help="rezip fasta files after analysis", 
    action="store_true")
misc_option_arguments.add_argument("--clean", 
    help="clean all intermediate files for maximum space efficiency - use with caution", 
    action="store_true")
misc_option_arguments.add_argument("--quiet", 
    help="suppress progress report and warnings", 
    action="store_true")
misc_option_arguments.add_argument("--dev", 
    help="keep intermediate files", 
    action="store_true")
misc_option_arguments.add_argument("--nonReductive", 
    help="non-reductive analysis. Dev mode - please use with caution", 
    action="store_true")
misc_option_arguments.add_argument("--f", 
    help="overwrite the analysis directory (provided as dir option) - please use with caution  ", 
    action="store_true")

args = ap.parse_args()

# relative path to absolute path
args.unmappedReads = os.path.abspath(args.unmappedReads)
args.dir = os.path.abspath(args.dir)

# "ONLY" option configuration: if none are selected, make everything true
if (not args.repeat and not args.immune and not args.circRNA 
  and not args.microbiome and not args.metaphlan):
    args.repeat = True
    args.immune = True
    args.circRNA = True
    args.microbiome = True
    args.metaphlan = True
else:
    args.nonReductive = True  # it is gonna be non-reductive for now 

    
################################################################################
### MAIN CODE
################################################################################

print """********************************************************************************
ROP is a computational protocol aimed to discover the source of all reads, 
originated from complex RNA molecules, recombinant antibodies and microbial 
communities. 

Written by Serghei Mangul (smangul@ucla.edu), Harry Taegyun Yang
(harry2416@gmail.com), Kevin Hsieh(kevin.hsieh@ucla.edu), Linus Chen (u6.30cl@gmail.com), University of California, Los Angeles (UCLA). (c) 2016.

Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see: https://sergheimangul.wordpress.com/rop/
ROP Tutorial: https://github.com/smangul1/rop/wiki 
********************************************************************************"""

# check if args.dir exists
if os.path.exists(args.dir) and not args.f:
    print """ERROR ::: The directory """ + args.dir + """ exists. Please choose a different directory to save results of the analysis. Alternatively choose --f option to
overwrite the results into """ + args.dir
    sys.exit(2)

if os.path.exists(args.dir) and args.f:
    cmd = "rm -fr " + args.dir + " &>/dev/null"
    os.system(cmd)

# database folder path
db_folder = "/db_" + args.organism

# basename (name of unmapped reads file without extension)
basename = os.path.splitext(os.path.basename(args.unmappedReads))[0]

# check if the input format is gzip and option --gzip provided
filename, file_extension = os.path.splitext(args.unmappedReads)
if file_extension == ".gz" and not args.gzip:
    print "ERROR: --gzip option is not selected with gzip format input. Please try with --gzip option"
    sys.exit(1)

################################################################################
# Prepare variables

# analysis directories
dirs = dict()
dirs["QC"] = args.dir + "/QC/"
dirs["lostReads"] = args.dir + "/lostReads/"
dirs["lostRepeat"] = args.dir + "/lostRepeatSequences/"
dirs["antibody"] = args.dir + "/antibodyProfile/"
dirs["bcr"] = dirs["antibody"] + "/BCR/"
dirs["tcr"] = dirs["antibody"] + "/TCR/"
dirs["NCL"] = args.dir + "/NCL/"
dirs["igh"] = dirs["bcr"] + "/IGH/"
dirs["igk"] = dirs["bcr"] + "/IGK/"
dirs["igl"] = dirs["bcr"] + "/IGL/"
dirs["tcra"] = dirs["tcr"] + "/TCRA/"
dirs["tcrb"] = dirs["tcr"] + "/TCRB/"
dirs["tcrd"] = dirs["tcr"] + "/TCRD/"
dirs["tcrg"] = dirs["tcr"] + "/TCRG/"
dirs["microbiome"] = args.dir + "/microbiomeProfile/"
dirs["metaphlan"] =  dirs["microbiome"] + "/metaphlan/"
dirs["bacteria"] = args.dir + "/microbiomeProfile/bacteriaProfile/"
dirs["virus"] = args.dir + "/microbiomeProfile/viralProfile/"
dirs["eupathdb"] = args.dir + "/microbiomeProfile/eukaryoticPathogenProfile/"
for dir in dirs:
    if not os.path.exists(dirs[dir]):
        os.makedirs(dirs[dir])

# intermediate file names
intfns = dict()
intfns["unmappedFastq"] = args.dir + "/unmapped_" + basename + ".fastq"
intfns["lowQFile"] = dirs["QC"] + basename + "_lowQ.fastq"
intfns["lowQFileFasta"] = dirs["QC"] + basename + "_lowQ.fa"
intfns["lowQCFile"] = dirs["QC"] + basename + "_lowQC.fa"
intfns["rRNAFile"] = dirs["QC"] + basename + "_rRNA_blastFormat6.csv"
intfns["afterrRNAFasta"] = dirs["QC"] + basename + "_after_rRNA.fasta"
intfns["afterlostReadsFasta"] = dirs["lostReads"] + basename + "_after_rRNA_lostReads.fasta"
intfns["afterlostReadsFastaGzip"] = dirs["lostReads"] + basename + "_after_rRNA_lostReads.fasta.gz"
intfns["afterNCLFasta"] = dirs["NCL"] + basename + "_after_NCL.fasta"
intfns["afterImmuneFasta"] = dirs["bcr"] + basename + "_afterImmune.fasta"
intfns["afterBacteriaFasta"] = dirs["bacteria"] + basename + "_afterBacteria.fasta"
intfns["afterVirusFasta"] = dirs["virus"] + basename + "_afterVirus.fasta"
intfns["unaccountedReadsFasta"] = args.dir + "/" + basename + "_unaccountedReads.fasta"
intfns["metaphlan_intermediate_map"] = dirs["metaphlan"] + basename + "_metaphlan.map"
intfns["metaphlan_intermediate_bowtie2out"] =  dirs["metaphlan"] + basename + "_bowtie2out.txt"
intfns["metaphlan_output"] =  dirs["metaphlan"]  +  basename  +  "_metaphlan_output.tsv"
intfns["gBamFile"] = dirs["lostReads"] + basename + "_genome.sam"
intfns["tBamFile"] = dirs["lostReads"] + basename + "_transcriptome.sam"
intfns["repeatFile"] = dirs["lostRepeat"] + basename + "_lostRepeats_blastFormat6.tsv"
intfns["afterlostRepeatFasta"] = dirs["lostRepeat"] + basename + "_after_lostRepeat.fasta"
intfns["NCL_CIRI_file"] = dirs["NCL"] + basename +  "_NCL_CIRI_after_bwa.sam"
intfns["after_NCL_CIRI_file_prefix"] = basename + "_circRNA.tsv"
intfns["ighFile"] = dirs["igh"] + basename + "_IGH_igblast.tsv"
intfns["igkFile"] = dirs["igk"] + basename + "_IGK_igblast.tsv"
intfns["iglFile"] = dirs["igl"] + basename + "_IGL_igblast.tsv"
intfns["tcraFile"] = dirs["tcra"] + basename + "_TCRA_igblast.tsv"
intfns["tcrbFile"] = dirs["tcrb"] + basename + "_TCRB_igblast.tsv"
intfns["tcrdFile"] = dirs["tcrd"] + basename + "_TCRD_igblast.tsv"
intfns["tcrgFile"] = dirs["tcrg"] + basename + "_TCRG_igblast.tsv"

# log file names
logfns = dict()
logfns["logQC"] = dirs["QC"] + basename + "_QC.log"
logfns["logrRNA"] = dirs["QC"] + basename + "_rRNA.log"
logfns["logHuman"] = dirs["lostReads"] + basename + "_lostHuman.log"
logfns["log_bowtieWG"] = dirs["lostReads"] + basename + "_bowtieWG.log"
logfns["log_bowtieTR"] = dirs["lostReads"] + basename + "_bowtieTR.log"
logfns["logLostRepeat"] = dirs["lostRepeat"] + basename + "_lostRepeat.log"
logfns["logNCL"] = dirs["NCL"] + basename + "_NCL.log"
logfns["bacteriaFile"] = dirs["bacteria"] + basename + "_bacteria_blastFormat6.csv"
logfns["virusFile"] = dirs["virus"] + basename + "_virus_blastFormat6.csv"
logfns["bacteriaFileFiltered"] = dirs["bacteria"] + basename + "_bacteriaFiltered_blastFormat6.csv"
logfns["virusFileFiltered"] = dirs["virus"] + basename + "_virusFiltered_blastFormat6.csv"
logfns["gLogfile"] = args.dir + "/" + basename + ".log"
logfns["cmdLogfile"] = args.dir + "/" + "dev.log"
logfns["toolsLogfile"] = args.dir + "/"+"tools.log"

logfns["logIGH"] = dirs["igh"] + basename + "_igh.log"
logfns["logIGL"] = dirs["igl"] + basename + "_igl.log"
logfns["logIGK"] = dirs["igk"] + basename + "_igk.log"
logfns["logTCRA"] = dirs["tcra"] + basename + "_tcra.log"
logfns["logTCRB"] = dirs["tcrb"] + basename + "_tcrb.log"
logfns["logTCRG"] = dirs["tcrg"] + basename + "_tcrg.log"
logfns["logTCRD"] = dirs["tcrd"] + basename + "_tcrd.log"
logfns["logMetaphlan"] = dirs["metaphlan"] + basename + "_metaphlan.log"
logfns["logBacteria"] = dirs["bacteria"] + basename + "_bacteria.log"
logfns["logVirus"] = dirs["virus"] + basename + "_virus.log"
logfns["logEukaryotes"] = dirs["eupathdb"] + basename + "_eukaryotes.log"

# run file names
runfns = dict()
runfns["runLostReadsFile"] = dirs["lostReads"] + "/runLostReads_" + basename + ".sh"
runfns["runLostRepeatFile"] = dirs["lostRepeat"] + "/runLostRepeat_" + basename + ".sh"
runfns["runNCL_CIRIfile"] = dirs["NCL"] + "/run_NCL_CIRI" + basename + ".sh" 
runfns["runIGHFile"] = dirs["igh"] + "/runIGH_" + basename + ".sh"
runfns["runIGKFile"] = dirs["igk"] + "/runIGK_" + basename + ".sh"
runfns["runIGLFile"] = dirs["igl"] + "/runIGL_" + basename + ".sh"
runfns["runTCRAFile"] = dirs["tcra"] + "/runTCRA_" + basename + ".sh"
runfns["runTCRBFile"] = dirs["tcrb"] + "/runTCRB_" + basename + ".sh"
runfns["runTCRDFile"] = dirs["tcrd"] + "/runTCRD_" + basename + ".sh"
runfns["runTCRGFile"] = dirs["tcrg"] + "/runTCRG_" + basename + ".sh"
runfns["runBacteriaFile"] = dirs["bacteria"] +"/runBacteria_" + basename + ".sh"
runfns["runVirusFile"] = dirs["virus"] + "/runVirus_" + basename + ".sh"
runfns["run_metaphlan_file"] = dirs["metaphlan"] + "/run_metaphlan_" + basename + ".sh"

################################################################################
# Prepare for analysis

# pysam info (for future use)
# pys = pysam.AlignmentFile(args.unmappedReads, "rb")
# for read in pys.fetch(until_eof=True):
#     ...

os.chdir(args.dir)
readLength=0
n=0

if args.b:
    if args.skipLowq:
        bam2fasta(cd,args.unmappedReads,intfns["lowQFileFasta"])
        unmapped_file = intfns["lowQFileFasta"]
    else:
        bam2fastq(cd,args.unmappedReads,intfns["unmappedFastq"])
        unmapped_file = intfns["unmappedFastq"]
elif args.gzip and not args.skipPreliminary and not args.skipQC and not args.skipLowq:
    write2Log("--gzip option is selected. The input is in gzip format. It will be decompressed in the analysis dir.", logfns["gLogfile"], args.quiet)
    intfns["unmappedFastq"]=args.dir+"/"+basename+".fastq"
    write_gzip_into_readable(args.unmappedReads, intfns["unmappedFastq"])
    unmapped_file = intfns["unmappedFastq"]
else:
    unmapped_file = args.unmappedReads

# Get Num Reads in unmapped fastq/fastas
if not args.skipPreliminary and not args.skipQC and not args.skipLowq:
    fastqfile = open(unmapped_file, "rU")
    for record in SeqIO.parse(fastqfile,"fastq"):
            readLength=len(record) #assumes the same length, will not work for Ion Torrent or Pac Bio
            n+=1
    fastqfile.close()
elif not args.gzip:
    fastafile = open(unmapped_file, "rU")
    for record in SeqIO.parse(fastafile,"fasta"):
        readLength=len(record) #assumes the same length, will not work for Ion Torrent or Pac Bio
        n+=1
    fastafile.close()
elif args.gzip:
    fastafile = gzip.open(unmapped_file, "rU")
    for record in SeqIO.parse(fastafile,"fasta"):
        readLength=len(record) #assumes the same length, will not work for Ion Torrent or Pac Bio
        n+=1
    fastafile.close()

message="Processing %s unmapped reads of length %s" %(n,readLength)


write2Log(message,logfns["gLogfile"],args.quiet)


if n==0:
    write2Log("ERROR ::: The input file is empty.",logfns["gLogfile"],0)
    sys.exit(1)

################################################################################
# 1. QC

nLowQReads=0
nLowCReads=0
n_rRNAReads=0
nlostHumanReads=0

if args.skipPreliminary:
    # intfns["afterrRNAFasta"]=args.unmappedReads
    write2Log("1. Quality Control is skipped",logfns["gLogfile"],args.quiet)
    write2Log("2. Remapping to references is skipped",logfns["gLogfile"],args.quiet)
    filename, file_extension = os.path.splitext(args.unmappedReads)
    print file_extension
    if file_extension!=".fa" and file_extension!=".fasta" and not args.gzip :
        write2Log("ERROR ::: --skipPreliminary option is selected. Reads needs to be in FASTA format",logfns["gLogfile"],args.quiet)
        sys.exit(1)
    elif file_extension == ".gz":
        temp=args.unmappedReads.split(".gz")[0]
        filename2, file_extension2 = os.path.splitext(temp)
        if file_extension2!=".fa" and file_extension2!=".fasta":
            write2Log("ERROR ::: --skipPreliminary option is selected. Reads needs to be in FASTA format",logfns["gLogfile"],args.quiet)
            sys.exit(1)
        if not args.gzip:
            write2Log("ERROR ::: --gzip option is not selected with gzip format input. Please try with --gzip option", logfns["gLogfile"], args.quiet)
    elif file_extension == ".gz" and args.gzip: 
        write2Log("SkipPreliminary option is selected. The input is in gzip format. It will be decompressed in the analysis dir.", logfns["gLogfile"], args.quiet)
elif args.skipQC:
    filename, file_extension = os.path.splitext(args.unmappedReads)
    if file_extension!=".fa" and file_extension!=".fasta":
        write2Log("ERROR ::: --skipQC option is selected. Reads needs to be in FASTA format",logfns["gLogfile"],args.quiet)
        sys.exit(1)
    # intfns["afterrRNAFasta"] = args.unmappedReads
    write2Log("1. Quality Control is skipped",logfns["gLogfile"],args.quiet)
else:
    if not args.b and not args.gzip:
        intfns["unmappedFastq"]=args.unmappedReads
    valid = 'ACTG'
    if not args.skipLowq:
        
        try:
        
            #lowQ
            write2Log("1. Quality Control...",logfns["gLogfile"],args.quiet)  
            fastafile=open(intfns["lowQFileFasta"], 'w')
            readLength=0
            nLowQReads=0
            nAfterLowQReads=0
            for record in SeqIO.parse(intfns["unmappedFastq"], "fastq"):
                readLength=len(record) #assumes the same length, will not work for Ion Torrent or Pac Bio
                j=record.letter_annotations["phred_quality"]
                prc=len([i for i in j if i>=20])/float(len(j))
                if prc>0.75 and all(i in valid for i in record.seq):
                    fastafile.write(str(">" + record.name) + "\n")
                    fastafile.write(str(record.seq) + "\n")
                    nAfterLowQReads+=1
            if args.b:
                os.remove(intfns["unmappedFastq"])
            fastafile.close()
            nLowQReads=n-nAfterLowQReads
            write2Log("--filtered %s low quality reads" % (nLowQReads) ,logfns["gLogfile"],args.quiet)
        except lowQaulityReadsError:
            write2Log("ERROR ::: ROP was not able to process low quality reads. Consider reporting the bug here : https://github.com/smangul1/rop/issues ",logfns["gLogfile"],0)
            sys.exit(1)
    else:
        write2Log("skipLowq option is selected - low quality filtering step is skipped" ,logfns["gLogfile"], args.quiet)
    os.chdir(dirs["QC"])
    try:
        #lowC
        cmd="export PATH=$PATH:%s/tools/seqclean-x86_64/bin" %(cd)
        os.system(cmd)
        cmd=cd+"/tools/seqclean-x86_64/seqclean %s -l 50 -M -o %s 2>>%s" %(intfns["lowQFileFasta"], intfns["lowQCFile"],logfns["logQC"])
        write2Log(cmd,logfns["cmdLogfile"],"False")
        os.system(cmd)
        cmd = "rm -rf %s/cleaning_1/ >temp 2 >temp; rm -f %s/*.cln >temp 2 >temp; rm -f %s/*.cidx>temp 2 >temp; rm -f %s/*.sort>temp 2 >temp" % (dirs["QC"],dirs["QC"],dirs["QC"],dirs["QC"])
        os.system(cmd)
        proc = subprocess.Popen(["grep trashed %s | awk -F \":\" '{print $2}'" %(logfns["logQC"]) ], stdout=subprocess.PIPE, shell=True)
        (nLowCReadsTemp, err) = proc.communicate()
        nLowCReads=int(nLowCReadsTemp.rstrip().strip())
        write2Log("--filtered %s low complexity reads (e.g. ACACACAC...)" %(nLowCReads) ,logfns["gLogfile"],args.quiet)
    except lowCReadsError:
        write2Log("ERROR ::: ROP was not able to process low complexity reads. Consider reporting the bug here : https://github.com/smangul1/rop/issues ",logfns["gLogfile"],0)
        sys.exit(1)


    #rRNA
    cmd="%s/tools/blastn -task megablast -index_name %s/%s/rRNA/rRNA -use_index true -query %s -db %s/%s/rRNA/rRNA  -outfmt 6 -evalue 1e-05  >%s" %(cd,cd,db_folder,intfns["lowQCFile"],cd,db_folder,intfns["rRNAFile"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    os.system(cmd)
    n_rRNATotal=0
    rRNAReads = set()
    with open(intfns["rRNAFile"],'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
            n_rRNATotal+=1
            element=line[0]
            identity=float(line[2])
            alignmentLength=float(line[3])
            eValue=float(line[10])
            if eValue<1e-05 and alignmentLength==readLength and identity>=0.94*readLength:
                rRNAReads.add(element)
    excludeReadsFromFasta(intfns["lowQCFile"],rRNAReads,intfns["afterrRNAFasta"])
    n_rRNAReads=len(rRNAReads)
    write2Log("--filtered %s rRNA reads" %(n_rRNAReads) ,logfns["gLogfile"],args.quiet)
    write2Log("In toto : %s reads failed QC and are filtered out" %(nLowQReads+nLowCReads+n_rRNAReads) ,logfns["gLogfile"],args.quiet)
    write2File("done!",args.dir+"/step1_QC.done")
    message="Number of entries in %s is %s" %(intfns["rRNAFile"],n_rRNATotal)
    write2Log(message,logfns["cmdLogfile"],"False")
    if not args.dev:
        os.remove(intfns["lowQCFile"])
        os.remove(intfns["lowQFileFasta"])
        os.remove(intfns["rRNAFile"])


################################################################################
# 2. Remapping to reference

if not args.skipPreliminary:
    write2Log("2. Remapping to references...",logfns["cmdLogfile"],"False")
    write2Log("2. Remapping to references...",logfns["gLogfile"],args.quiet)
    
    # If input is afterQC fasta.gz
    if args.skipQC and args.outGz:
        write_gzip_into_readable(args.unmappedReads, intfns["afterrRNAFasta"])
    elif args.skipQC and not args.outGz:
        intfns["afterrRNAFasta"] = args.unmappedReads
    cmdGenome="%s/tools/bowtie2 -k 1 -p 8 -f -x %s/%s/bowtie2Index/genome -U %s 2>>%s | %s/tools/samtools view -SF4 -   >%s" %(cd,cd,db_folder,intfns["afterrRNAFasta"],logfns["log_bowtieWG"],cd,intfns["gBamFile"])

    #transcriptome
    cmdTranscriptome="%s/tools/bowtie2  -k 1 -f -p 8 -x %s/%s/bowtie2Index/hg19KnownGene.exon_polya200 -U %s 2>>%s | %s/tools/samtools view -SF4 -  >  %s " %(cd,cd,db_folder,intfns["afterrRNAFasta"],logfns["log_bowtieTR"], cd,intfns["tBamFile"])
    write2Log(cmdGenome,logfns["cmdLogfile"],"False")
    write2Log(cmdTranscriptome,logfns["cmdLogfile"],"False")
    os.system(cmdTranscriptome)
    os.system(cmdGenome)
    lostHumanReads = set()
    lostHumanReads0= set()
    lostHumanReads1= set()
    lostHumanReads2= set()
    with open(intfns["tBamFile"],'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
            mismatch=int(line[16].split(':')[2])
            if mismatch<3:
                lostHumanReads.add(line[0])
            if mismatch==0:
                lostHumanReads0.add(line[0])
            elif mismatch==1:
                lostHumanReads1.add(line[0])
            elif mismatch==2:
                lostHumanReads2.add(line[0])
    with open(intfns["gBamFile"],'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
                mismatch=int(line[16].split(':')[2])
                if mismatch<3:
                    lostHumanReads.add(line[0])
                if mismatch==0:
                    lostHumanReads0.add(line[0])
                elif mismatch==1:
                    lostHumanReads1.add(line[0])
                elif mismatch==2:
                    lostHumanReads2.add(line[0])
    write2Log("Unmapped reads mapped to genome and/or transcriptome (using bowtie2) are categorized as lost reads and are excluded from the further analysis. This includes :  %s reads with 0 mismatches; %s reads with 1 mismatch; %s reads with 2 mismatches" %(len(lostHumanReads0),len(lostHumanReads1),len(lostHumanReads2)),logfns["logHuman"],"False")
    write2Log("Complete list of lost reads is available from sam files: %s,%s" %(intfns["gBamFile"],intfns["tBamFile"]),logfns["logHuman"],"False")
    if args.outGz:
        excludeReadsFromFastaGzip(intfns["afterrRNAFasta"],lostHumanReads,intfns["afterlostReadsFastaGzip"])
    else:
        excludeReadsFromFasta(intfns["afterrRNAFasta"],lostHumanReads,intfns["afterlostReadsFasta"])
    nlostHumanReads=len(lostHumanReads)
    write2Log("--identified %s lost reads from unmapped reads. Among those: %s reads with 0 mismatches; %s reads with 1 mismatch; %s reads with 2 mismatches" %(len(lostHumanReads),len(lostHumanReads0),len(lostHumanReads1),len(lostHumanReads2)), logfns["gLogfile"], args.quiet)
    write2Log("***Note: Complete list of lost reads is available from sam files: %s,%s" %(intfns["gBamFile"],intfns["tBamFile"]), logfns["gLogfile"], args.quiet)
    write2File("done!",args.dir+"/step2_lostHumanReads.done")
    if not args.dev:
        if not args.skipQC and not args.skipPreliminary:
            os.remove(intfns["afterrRNAFasta"])
    if args.clean:
        write2Log("Clean mode selected - removing analysis sam files", logfns["gLogfile"], args.quiet)
        os.remove(intfns["gBamFile"])
        os.remove(intfns["tBamFile"])

# TODO
else:
    if args.outGz:
        write_gzip_into_readable(args.unmappedReads, intfns["afterlostReadsFasta"])
    else:
        intfns["afterlostReadsFasta"] = args.unmappedReads

if args.nonReductive or args.qsub or args.qsubArray:
    branch_point_file = intfns["afterlostReadsFasta"]
    write2Log("*********************************",logfns["gLogfile"], args.quiet)
    write2Log("Non-substractive mode is selected : Low quality, low complexity, rRNA reads and lost reads are filtered out. Resulting high quality non-matching reads are provided as input  for STEP3-STEP6",logfns["gLogfile"], args.quiet)
    write2Log("*********************************",logfns["gLogfile"], args.quiet)
    

################################################################################
# 3. Mapping to repeat sequences

if args.repeat:
    write2Log("3. Mapping to repeat sequences...",logfns["cmdLogfile"],"False")
    write2Log("3. Mapping to repeat sequences...",logfns["gLogfile"],args.quiet)

    #TO DO : make all fasta ->gzip
    #gzip -dc %s | , query -
    # CHANGED 
    # cmd="%s/tools/blastn -task megablast -index_name %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa -use_index true -query %s -db %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa  -outfmt 6 -evalue 1e-05  > %s" %(cd, cd, intfns["afterlostReadsFasta"], cd, intfns["repeatFile"])
    
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = intfns["afterlostReadsFasta"]
    #cmd="%s/tools/blastn -task megablast -index_name %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa -use_index true -query %s -db %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa  -outfmt 6 -evalue 1e-05  > %s 2>%s" %(cd, cd, input_file, cd, intfns["repeatFile"],logLostRepeat)
	cmd="%s/tools/blastn -task megablast -index_name %s/%s/repeats/repbase.fa -use_index true -query %s -db %s/%s/repeats/repbase.fa  -outfmt 6 -evalue 1e-05  > %s 2>%s" %(cd, cd, db_folder, input_file, cd, db_folder, intfns["repeatFile"],logfns["logLostRepeat"])
    if args.qsub or args.qsubArray:
        f = open(logfns["runLostRepeatFile"],'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_lostRepeat.done \n" %(dirs["lostRepeat"],basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runLostRepeatFile"])
            else:
                cmdQsub="qsub -cwd -V -N lostRepeat -l h_data=16G,time=10:00:00 %s" %(logfns["runLostRepeatFile"])
            os.system(cmdQsub)
            write2Log("Job for STEP3 has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.system(cmd)
        write2Log(cmd,logfns["cmdLogfile"],"False")
    if not args.qsub and not args.qsubArray:
        lostRepeatReads = set()
        with open(intfns["repeatFile"],'r') as f:
            reader=csv.reader(f,delimiter='\t')
            for line in reader:
                element=line[0]
                identity=float(line[2])
                alignmentLength=float(line[3])
                eValue=float(line[10])
                if eValue<1e-05 and alignmentLength>=0.8*readLength and identity>=0.9*readLength:
                    lostRepeatReads.add(element)
        nRepeatReads=len(lostRepeatReads)
        write2Log("-- Identified %s lost repeat sequences from unmapped reads" %(nRepeatReads) ,logfns["gLogfile"],args.quiet)
        write2Log("***Note : Repeat sequences classification into classes (e.g. LINE) and families (e.g. Alu) will be available in next release" ,logfns["gLogfile"],args.quiet)
        excludeReadsFromFasta(intfns["afterlostReadsFasta"],lostRepeatReads,intfns["afterlostRepeatFasta"])
        write2File("done!",args.dir+"/step3_lostRepeatSequences.done")
        if not args.dev:
            if not args.skipPreliminary:
                os.remove(intfns["afterlostReadsFasta"])
        write2Log("Lost repeat reads are mapped to the repeat sequences (using megablast)",logfns["logLostRepeat"],"False")
        write2Log("Complete list of lost repeat reads is available from tsv file: %s" %(intfns["repeatFile"]),logfns["logLostRepeat"],"False")
else:
    print "3. Maping to repeat sequences is skipped."


################################################################################
# 4. Non-co-linear RNA profiling

write2Log("4. Non-co-linear RNA profiling",logfns["cmdLogfile"],"False")
write2Log("4. Non-co-linear RNA profiling",logfns["gLogfile"],args.quiet)
write2Log("***Note : Trans-spicing and gene fusions  are currently not supported, but will be in the next release.",logfns["gLogfile"],args.quiet)

os.chdir(dirs["NCL"])
NCL_reads=set()
nNCLReads=0

if args.circRNA:
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = intfns["afterlostRepeatFasta"]
    cmd="%s/tools/bwa mem -T -S %s/%s/BWAIndex/genome.fa %s > %s 2>%s \n" %(cd,cd,db_folder,input_file,intfns["NCL_CIRI_file"],logfns["logNCL"])
    cmd = cmd + "perl %s/tools/CIRI_v1.2.pl -S -I %s -O %s -F %s/%s/BWAIndex/genome.fa 1>>%s 2>>%s" %(cd,intfns["NCL_CIRI_file"],intfns["after_NCL_CIRI_file_prefix"],cd,db_folder,logfns["logNCL"],logfns["logNCL"])
    if args.qsub or args.qsubArray:
        f = open(logfns["runNCL_CIRIfile"],'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_NCL_CIRI.done \n" %(dirs["NCL"],basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runNCL_CIRIfile"])
            else:
                cmdQsub="qsub -cwd -V -N NCL_CIRI -l h_data=8G,time=10:00:00 %s" %(logfns["runNCL_CIRIfile"])
            os.system(cmdQsub)
            write2Log("Job for STEP4 has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.chdir(dirs["NCL"])
        os.system(cmd)
        NCL_reads=nCircularReads(intfns["after_NCL_CIRI_file_prefix"])
        nReadsNCL=len(NCL_reads)
        write2Log("--identified %s reads from circRNA" %(nReadsNCL) ,logfns["gLogfile"],args.quiet)
        write2Log("***Note: circRNAs detected by CIRI are available here: %s" %(intfns["after_NCL_CIRI_file_prefix"]) ,logfns["gLogfile"],args.quiet)
        excludeReadsFromFasta(input_file,NCL_reads,intfns["afterNCLFasta"])
        nNCLReads=len(NCL_reads)
        write2File("done!",args.dir+"/step4_NCL.done")
        if not args.dev:
            os.remove(input_file)
else:
    print "4. Non-co-linear RNA profiling is skipped."

    
################################################################################
# 5. T and B lymphocytes profiling

if args.immune:
    immuneReads = set()
    write2Log("5a. B lymphocytes profiling", logfns["cmdLogfile"], "False")
    write2Log("5a. B lymphocytes profiling", logfns["gLogfile"], args.quiet)
    
    """
    # All B_lymphocytes (work in progress)
    
    B_lymphocytes = ["igh", "igk", "igl"]
    
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = intfns["afterNCLFasta"]

    for current in B_lymphocytes:
        os.chdir(dirs[current])
        cmd="ln -s " + cd + "/" + db_folder + "/antibody/internal_data ./"
        if subprocess.Popen([cmd], shell=True).wait():
            sys.exit(2)
        write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of immunoglobulin heavy locus (" + current.upper() + ")",logfns["log" + current.upper()],"False")
        write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logfns["log" + current.upper()],"False")
        
        aux_string = "" if args.organism == "human" else " -auxiliary_data " + cd + "/" + \
            db_folder + "/antibody/optional_file/" + args.organism + "_gl.aux"
            
        cmd = cd + "/tools/igblastn -organism " + args.organism + aux_string + \
          " -germline_db_V " + cd + "/" + db_folder + "/antibody/" + current.upper() + "V.fa" + \
          " -germline_db_D " + cd + "/" + db_folder + "/antibody/" + current.upper() + "D.fa" + \
          " -germline_db_J " + cd + "/" + db_folder + "/antibody/" + current.upper() + "J.fa" + \
          " -query " + input_file + " -outfmt '7 std qseq sseq' -evalue 1e-05" + \
          " 2>>" + logfns["log" + current.upper()] + \
          " | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }'" + \
          " >" + intfns[current + "File"]
        
        write2Log(cmd, logfns["cmdLogfile"],"False")
        if args.qsub or args.qsubArray:
            f = open(logfns["run" + current.upper() + "File"],'w')
            f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(cd,db_folder))
            f.write(cmd+"\n")
            f.write("echo \"done!\"> %s/%s_" + current + ".done \n" %(dirs[current],basename))
            f.close()
            if args.qsub:
                if args.maui:
                    cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["run" + current.upper() + "File"])
                else:
                    cmdQsub="qsub -cwd -V -N " + current + " -l h_data=16G,time=24:00:00 %s" %(logfns["run" + current.upper() + "File"])
                os.system(cmdQsub)
                write2Log("Job for STEP5a(" + current.upper() + ") has been submitted via qsub",logfns["gLogfile"], args.quiet)
        else:
            os.chdir(dirs[current])
            os.system(cmd)
            immuneReadsCurrent=nReadsImmune(intfns[current + "File"])
            nReadsImmuneCurrent=len(immuneReadsCurrent)
            write2Log("--identified %s reads mapped to "%(nReadsImmuneCurrent) + current.upper() + " locus"  ,logfns["gLogfile"],args.quiet)
    """
    
    # IGH
    os.chdir(dirs["igh"])
    cmd="ln -s %s//%s/antibody/internal_data/ ./" %(cd,db_folder)
    os.system(cmd)
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = intfns["afterNCLFasta"]
    write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of immunoglobulin heavy locus (IGH)",logfns["logIGH"],"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logfns["logIGH"],"False")
    if args.organism == "human":
        cmd = "%s/tools/igblastn -germline_db_V %s/%s/antibody/IGHV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGHJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05  2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(cd,cd,db_folder,cd,db_folder,cd,db_folder,input_file, logfns["logIGH"],intfns["ighFile"])
    else:
        cmd = cd + "/tools/igblastn -organism mouse -auxiliary_data " + cd + "/" + db_folder + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/IGHV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGHJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05  2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(cd,db_folder,cd,db_folder,cd,db_folder,input_file, logfns["logIGH"],intfns["ighFile"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    if args.qsub or args.qsubArray:
        f = open(logfns["runIGHFile"],'w')
        f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(cd,db_folder))
        f.write(cmd+"\n")
        f.write("echo \"done!\"> %s/%s_igh.done \n" %(dirs["igh"],basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runIGHFile"])
            else:
                cmdQsub="qsub -cwd -V -N igh -l h_data=16G,time=24:00:00 %s" %(logfns["runIGHFile"])
            os.system(cmdQsub)
            write2Log("Job for STEP5a(IGK) has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.chdir(dirs["igh"])
        os.system(cmd)
        immuneReadsIGH=nReadsImmune(intfns["ighFile"])
        nReadsImmuneIGH=len(immuneReadsIGH)
        write2Log("--identified %s reads mapped to immunoglobulin heavy (IGH) locus" %(nReadsImmuneIGH) ,logfns["gLogfile"],args.quiet)

    #IGK    
    os.chdir(dirs["igk"])
    cmd="ln -s %s//%s/antibody/internal_data/ ./" %(cd,db_folder)
    os.system(cmd)
    write2Log("IgBlast	was	used	to	identify	reads	spanning	VJ recombinations of immunoglobulin kappa locus (IGK)",logfns["logIGK"],"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logfns["logIGK"],"False")
    if args.organism == "human":
        cmd = "%s/tools/igblastn -germline_db_V %s/%s/antibody/IGKV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGKJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05  2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(cd,cd,db_folder,cd,db_folder,cd,db_folder,input_file, logfns["logIGK"],intfns["igkFile"])
    else:
        cmd = cd + "/tools/igblastn -organism mouse -auxiliary_data " + cd + "/" + db_folder + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/IGKV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGKJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05  2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(cd,db_folder,cd,db_folder,cd,db_folder,input_file, logfns["logIGK"],intfns["igkFile"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    if args.qsub or args.qsubArray:
        f = open(runintfns["igkFile"],'w')
        f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(cd,db_folder))
        f.write(cmd+"\n")
        f.write("echo \"done!\"> %s/%s_igk.done \n" %(dirs["igk"],basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runIGKFile"])
            else:
                cmdQsub="qsub -cwd -V -N igk -l h_data=16G,time=24:00:00 %s" %(logfns["runIGKFile"])
            os.system(cmdQsub)
            write2Log("Job for STEP5a(IGK) has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.chdir(dirs["igk"])
        os.system(cmd)
        immuneReadsIGK=nReadsImmune(intfns["igkFile"])
        nReadsImmuneIGK=len(immuneReadsIGK)
        write2Log("--identified %s reads mapped to immunoglobulin kappa (IGK) locus " %(nReadsImmuneIGK) ,logfns["gLogfile"],args.quiet)
                
    # IGL    
    os.chdir(dirs["igl"])
    cmd="ln -s %s//%s/antibody/internal_data/ ./" %(cd,db_folder)
    os.system(cmd)
    write2Log("IgBlast	was	used	to	identify	reads	spanning	VJ recombinations of immunoglobulin lambda locus (IGL)",logfns["logIGL"],"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logfns["logIGL"],"False")
    if args.organism == "human":
        cmd = "%s/tools/igblastn -germline_db_V %s/%s/antibody/IGLV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGLJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(cd,cd,db_folder,cd,db_folder,cd,db_folder,input_file,logfns["logIGL"],intfns["iglFile"])
    else:
        cmd = cd + "/tools/igblastn -organism mouse -auxiliary_data " + cd + "/" + db_folder + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/IGLV.fa -germline_db_D %s/%s/antibody/IGHD.fa  -germline_db_J %s/%s/antibody/IGLJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(cd,db_folder,cd,db_folder,cd,db_folder,input_file,logfns["logIGL"],intfns["iglFile"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    if args.qsub or args.qsubArray:
        f = open(logfns["runIGLFile"],'w')
        f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(cd,db_folder))
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_igl.done \n" %(dirs["igl"],basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runIGLFile"])
            else:
                cmdQsub="qsub -cwd -V -N igl -l h_data=16G,time=24:00:00 %s" %(logfns["runIGLFile"])
            os.system(cmdQsub)
            write2Log("Job for STEP5a(IGL) has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.chdir(dirs["igl"])
        os.system(cmd)
        immuneReadsIGL=nReadsImmune(intfns["iglFile"])
        nReadsImmuneIGL=len(immuneReadsIGL)
        write2Log("--identified %s reads mapped to immunoglobulin lambda (IGL) locus" %(nReadsImmuneIGL) ,logfns["gLogfile"],args.quiet)

    write2Log("5b. T lymphocytes profiling...",logfns["cmdLogfile"],"False")
    write2Log("5b. T lymphocytes profiling...",logfns["gLogfile"],args.quiet)
    
    # TCRA    
    os.chdir(dirs["tcra"])
    cmd="ln -s %s//%s/antibody/internal_data/ ./" %(cd,db_folder)
    os.system(cmd)
    write2Log("IgBlast	was	used	to	identify	reads	spanning	VJ recombinations of T cell receptor alpha locus (TCRA)",logfns["logTCRA"],"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logfns["logTCRA"],"False")
    if args.organism == "human":
        cmd = "%s/tools/igblastn -ig_seqtype TCR -germline_db_V %s/%s/antibody/TRAV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRAJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(cd,cd,db_folder,cd,db_folder,cd,db_folder,input_file,logfns["logTCRA"],intfns["tcraFile"])
    else:
        cmd = cd + "/tools/igblastn -organism mouse -ig_seqtype TCR -auxiliary_data " + cd + "/" + db_folder + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/TRAV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRAJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(cd,db_folder,cd,db_folder,cd,db_folder,input_file,logfns["logTCRA"],intfns["tcraFile"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    if args.qsub or args.qsubArray:
        f = open(logfns["runTCRAFile"],'w')
        f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(cd,db_folder))
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_tcra.done \n"%(dirs["tcra"],basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runTCRAFile"])
            else:
                cmdQsub="qsub -cwd -V -N tcra -l h_data=16G,time=24:00:00 %s" %(logfns["runTCRAFile"])
            os.system(cmdQsub)
            write2Log("Job for STEP5b(TCRA) has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.chdir(dirs["tcra"])
        os.system(cmd)
        immuneReadsTCRA=nReadsImmune(intfns["tcraFile"])
        nReadsImmuneTCRA=len(immuneReadsTCRA)
        write2Log("--identified %s reads mapped to T cell receptor alpha (TCRA) locus" %(nReadsImmuneTCRA) ,logfns["gLogfile"],args.quiet) 

    # TCRB
    os.chdir(dirs["tcrb"])
    cmd="ln -s %s//%s/antibody/internal_data/ ./" %(cd, db_folder)
    os.system(cmd)
    write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of T cell receptor beta locus (TCRB)",logfns["logTCRB"],"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logfns["logTCRB"],"False")
    if args.organism == "human":
        cmd = "%s/tools/igblastn -ig_seqtype TCR -germline_db_V %s/%s/antibody/TRBV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRBJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(cd,cd,db_folder,cd,db_folder,cd,db_folder,input_file,logfns["logTCRB"],intfns["tcrbFile"])
    else:
        cmd = cd + "/tools/igblastn -organism mouse -ig_seqtype TCR -auxiliary_data " + cd + "/" + db_folder + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/TRBV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRBJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(cd,db_folder,cd,db_folder,cd,db_folder,input_file,logfns["logTCRB"],intfns["tcrbFile"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    if args.qsub or args.qsubArray:
        f = open(logfns["runTCRBFile"],'w')
        f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(cd, db_folder))
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_tcrb.done \n"%(dirs["tcrb"],basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runTCRBFile"])
            else:
                cmdQsub="qsub -cwd -V -N tcrb -l h_data=16G,time=24:00:00 %s" %(logfns["runTCRBFile"])
            os.system(cmdQsub)
            write2Log("Job for STEP5b(TCRB) has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.chdir(dirs["tcrb"])
        os.system(cmd)
        immuneReadsTCRB=nReadsImmune(intfns["tcrbFile"])
        nReadsImmuneTCRB=len(immuneReadsTCRB)
        write2Log("--identified %s reads mapped to T cell receptor beta (TCRB) locus" %(nReadsImmuneTCRB) ,logfns["gLogfile"],args.quiet)
        
    #TCRD    
    os.chdir(dirs["tcrd"])
    cmd="ln -s %s//%s/antibody/internal_data/ ./" %(cd, db_folder)
    os.system(cmd)
    write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of T cell receptor delta locus (TCRD)",logfns["logTCRD"],"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logfns["logTCRD"],"False")
    if args.organism == "human":
        cmd = "%s/tools/igblastn -ig_seqtype TCR -germline_db_V %s/%s/antibody/TRDV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRDJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(cd,cd,db_folder,cd,db_folder,cd,db_folder,input_file,logfns["logTCRD"],intfns["tcrdFile"])
    else:
        cmd = cd + "/tools/igblastn -organism mouse -ig_seqtype TCR -auxiliary_data " + cd + "/" + db_folder + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/TRDV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRDJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(cd,db_folder,cd,db_folder,cd,db_folder,input_file,logfns["logTCRD"],intfns["tcrdFile"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    if args.qsub or args.qsubArray:
        f = open(logfns["runTCRDFile"],'w')
        f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(cd, db_folder))
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_tcrd.done \n" %(dirs["tcrd"], basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runTCRDFile"])
            else:
                cmdQsub="qsub -cwd -V -N tcrd -l h_data=16G,time=24:00:00 %s" %(logfns["runTCRDFile"])
            os.system(cmdQsub)
            write2Log("Job for STEP5b(TCRD) has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.chdir(dirs["tcrd"])
        os.system(cmd)
        immuneReadsTCRD=nReadsImmune(intfns["tcrdFile"])
        nReadsImmuneTCRD=len(immuneReadsTCRD)
        write2Log("--identified %s reads mapped to T cell receptor delta (TCRD) locus" %(nReadsImmuneTCRD) ,logfns["gLogfile"],args.quiet)
                
    # TCRG    
    os.chdir(dirs["tcrg"])
    cmd="ln -s %s//%s/antibody/internal_data/ ./" %(cd, db_folder)
    os.system(cmd)
    write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of T cell receptor gamma locus (TCRG)",logfns["logTCRG"],"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logfns["logTCRG"],"False")
    if args.organism == "human":
        cmd = "%s/tools/igblastn -ig_seqtype TCR -germline_db_V %s/%s/antibody/TRGV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRGJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(cd,cd,db_folder,cd,db_folder,cd,db_folder,input_file,logfns["logTCRG"],intfns["tcrgFile"])
    else:
        cmd = cd + "/tools/igblastn -organism mouse -ig_seqtype TCR -auxiliary_data " + cd + "/" + db_folder + "/antibody/optional_file/mouse_gl.aux -germline_db_V %s/%s/antibody/TRGV.fa -germline_db_D %s/%s/antibody/TRBD.fa  -germline_db_J %s/%s/antibody/TRGJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(cd,db_folder,cd,db_folder,cd,db_folder,input_file,logfns["logTCRG"],intfns["tcrgFile"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    if args.qsub or args.qsubArray:
        f = open(logfns["runTCRGFile"],'w')
        f.write("ln -s %s//%s/antibody/internal_data/ ./ \n" %(cd, db_folder))
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_tcrg.done \n" %(dirs["tcrg"],basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runTCRGFile"])
            else:
                cmdQsub="qsub -cwd -V -N tcrg -l h_data=16G,time=24:00:00 %s" %(logfns["runTCRGFile"])
            os.system(cmdQsub)
            write2Log("Job for STEP5b(TCRG) has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.chdir(dirs["tcrg"])
        os.system(cmd)
        immuneReadsTCRG=nReadsImmune(intfns["tcrgFile"])
        nReadsImmuneTCRG=len(immuneReadsTCRG)
        write2Log("--identified %s reads mapped to T cell receptor gamma locus (TCRG) locus" %(nReadsImmuneTCRG) ,logfns["gLogfile"],args.quiet)
    nReadsImmuneTotal=0
    if not args.qsub and not args.qsubArray:
        nReadsImmuneTotal=nReadsImmuneIGH+nReadsImmuneIGL+nReadsImmuneIGK+nReadsImmuneTCRA+nReadsImmuneTCRB+nReadsImmuneTCRD+nReadsImmuneTCRG
        write2Log("In toto : %s reads mapped to antibody repertoire loci" %(nReadsImmuneTotal) ,logfns["gLogfile"],args.quiet)
        write2Log("***Note : Combinatorial diversity of the antibody repertoire (recombinations of the of VJ gene segments)  will be available in the next release.",logfns["gLogfile"],args.quiet)
        write2File("done!",args.dir+"/step5_antibodyProfile.done")
        immuneReads=set().union(immuneReadsTCRA,immuneReadsTCRB,immuneReadsTCRD,immuneReadsTCRG)
        excludeReadsFromFasta(input_file,immuneReads,intfns["afterImmuneFasta"])
        if not args.dev:
            if args.circRNA:           
                os.remove(intfns["afterNCLFasta"])
else:
    print "5a. B lymphocytes profiling is skipped."
    print "5b. T lymphocytes profiling is skipped."
    
################################################################################
# 5. Metaphlan

if args.microbiome:
    write2Log("***Extra step.  Metaphlan profiling...",logfns["cmdLogfile"],"False")
    write2Log("***Extra step.  Metaphlan profiling...",logfns["gLogfile"],args.quiet)
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = intfns["afterImmuneFasta"]
    cmd = "python %s/tools/metaphlan2.py %s %s --mpa_pkl %s/%s/metaphlan/mpa_v20_m200.pkl --bowtie2_exe %s/tools/bowtie2 --input_type multifasta --bowtie2db %s/%s/metaphlan/mpa_v20_m200 -t reads_map --nproc 8 --bowtie2out %s>>%s 2>>%s" % (cd, input_file, intfns["metaphlan_intermediate_map"], cd, db_folder, cd, cd, db_folder, intfns["metaphlan_intermediate_bowtie2out"],logfns["logMetaphlan"],logfns["logMetaphlan"])
    cmd = cmd + "\n" + "python %s/tools/metaphlan2.py --mpa_pkl %s/%s/metaphlan/mpa_v20_m200.pkl --bowtie2_exe %s/tools/bowtie2 --input_type bowtie2out %s -t rel_ab > %s 2>>%s" %(cd, cd, db_folder, cd, intfns["metaphlan_intermediate_bowtie2out"], intfns["metaphlan_output"],logfns["logMetaphlan"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    if args.qsub or args.qsubArray:
        f= open(logfns["run_metaphlan_file"], 'w')
        f.write(cmd + "\n")
        f.write("echo \"done!\" > %s/%s_metaphlan.done \n" % (dirs["metaphlan"], basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["run_metaphlan_file"])
            else:
                cmdQsub="qsub -cwd -V -N metaphlan -l h_data=16G,time=24:00:00 %s" %(logfns["run_metaphlan_file"])
            os.system(cmdQsub)
            write2Log("Job for Metaphlan has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.chdir(dirs["metaphlan"])
        os.system(cmd)
        write2Log("***Microbiome profiling by Metaphlan2: taxonomic profile of microbial communities detected by Metaphlan2 is available here: %s" %(dirs["metaphlan"]) ,logfns["gLogfile"],args.quiet)
else:
    print "Extra step.  Metaphlan profiling is skipped."


################################################################################
# 6. Microbiome profiling

if args.microbiome:
    write2Log("6.  Microbiome profiling...",logfns["cmdLogfile"],"False")
    write2Log("6.  Microbiome profiling...",logfns["gLogfile"],args.quiet)
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = intfns["afterImmuneFasta"]

    # bacteria    
    os.chdir(dirs["bacteria"])
    readLength=0
    if not args.outGz:
        fastafile = open(input_file, "rU")
    else:
        fastafile = gzip.open(input_file, "rU")
    for record in SeqIO.parse(fastafile,"fasta"):
        readLength=len(record) #assumes the same length, will not work for Ion Torrent or Pac Bio
        break;
    write2Log("Megablast was used to map the reads onto the	bacterial reference	genomes ",logfns["logBacteria"],"False")
    write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides aligned are considered microbial reads and are saved into %s " %(logfns["bacteriaFileFiltered"]), logfns["logBacteria"],"False")
    write2Log("---------------",logfns["logBacteria"],"False")
    cmd="%s/tools/blastn -task megablast -index_name %s/%s/bacteria/bacteria -use_index true -query %s -db %s/%s/bacteria/bacteria  -outfmt 6 -evalue 1e-05  >%s 2>>%s" %(cd,cd,db_folder,input_file,cd,db_folder,logfns["bacteriaFile"],logfns["logBacteria"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    if args.qsub or args.qsubArray:
        f = open(logfns["runBacteriaFile"],'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_bacteria.done \n"%(dirs["bacteria"], basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runBacteriaFile"])
            else:
                cmdQsub="qsub -cwd -V -N bacteria -l h_data=16G,time=24:00:00 %s" %(logfns["runBacteriaFile"])
            os.system(cmdQsub)
            write2Log("Job for STEP6 (bacteriaProfile) has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.chdir(dirs["bacteria"])
        os.system(cmd)
        bacteriaReads=nMicrobialReads(logfns["bacteriaFile"],readLength,logfns["bacteriaFileFiltered"])
        nReadsBacteria=len(bacteriaReads)
        write2Log("--identified %s reads mapped bacterial genomes" %(nReadsBacteria) ,logfns["gLogfile"],args.quiet)
        excludeReadsFromFasta(input_file,bacteriaReads,intfns["afterBacteriaFasta"])
    
    # virus    
    os.chdir(dirs["virus"])
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = intfns["afterBacteriaFasta"]
        if not args.dev:
            os.remove(intfns["afterImmuneFasta"])
    write2Log("Megablast was used to map the reads onto the	viral reference	genomes ",logfns["logVirus"],"False")
    write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides aligned are considered microbial reads and are saved into %s" %(logfns["virusFileFiltered"]),logfns["logVirus"],"False")
    write2Log("----------------------",logfns["logVirus"],"False")
    cmd="%s/tools/blastn -task megablast -index_name %s/%s/virus/viruses -use_index true -query %s -db %s/%s/virus/viruses  -outfmt 6 -evalue 1e-05  >%s 2>>%s" %(cd,cd,db_folder,input_file,cd,db_folder,logfns["virusFile"],logfns["logVirus"])
    write2Log(cmd,logfns["cmdLogfile"],"False")
    if args.qsub or args.qsubArray:
        f = open(logfns["runVirusFile"],'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_virus.done \n"%(dirs["bacteria"], basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(logfns["runBacteriaFile"])
            else:
                cmdQsub="qsub -cwd -V -N virus -l h_data=16G,time=24:00:00 %s" %(logfns["runBacteriaFile"])
            os.system(cmdQsub)
            write2Log("Job for STEP6 (viralProfile) has been submitted via qsub",logfns["gLogfile"], args.quiet)
    else:
        os.system(cmd)
        virusReads=nMicrobialReads(logfns["virusFile"],readLength,logfns["virusFileFiltered"])
        nReadsVirus=len(virusReads)
        write2Log("--identified %s reads mapped viral genomes" %(nReadsVirus) ,logfns["gLogfile"],args.quiet)
        excludeReadsFromFasta(input_file,virusReads,intfns["afterVirusFasta"])

    # eukaryotic pathogens    
    os.chdir(dirs["eupathdb"])
    dbList=["ameoba",
            "crypto",
            "giardia",
            "microsporidia",
            "piroplasma",
            "plasmo",
            "toxo",
            "trich",
            "tritryp"]
    if args.nonReductive or args.qsub or args.qsubArray:
        inFasta = branch_point_file
    else:
        inFasta = intfns["afterVirusFasta"]
        if not args.dev:
            os.remove(intfns["afterBacteriaFasta"])
    nReadsEP=0
    for db in dbList:
        afterFasta=dirs["eupathdb"]+"%s_after_%s.fasta" %(basename,db)
        eupathdbFile=dirs["eupathdb"]+basename+"_"+db+"_blastFormat6.csv"
        eupathdbFileFiltered=dirs["eupathdb"]+basename+"_"+db+"Filtered_blastFormat6.csv"
        runEupathdbFile=dirs["eupathdb"]+"/run_"+basename+"_"+db+".sh"
        write2Log("Megablast was used to map the reads onto the	reference genomes of Eukaryotes (http://eupathdb.org/eupathdb/)",logfns["logEukaryotes"],"False")
        write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides aligned are considered microbial reads and are saved into %s" %(afterFasta), logfns["logEukaryotes"],"False")
        write2Log("-------------",logfns["logEukaryotes"],"False")
        cmd="%s/tools/blastn -task megablast -index_name %s/%s/eupathdb/%s -use_index true -query %s -db %s/%s/eupathdb/%s  -outfmt 6 -evalue 1e-05  >%s 2>>%s" %(cd,cd,db_folder,db,inFasta,cd,db_folder,db,eupathdbFile,logfns["logEukaryotes"])
        write2Log(cmd,logfns["cmdLogfile"],"False")
        if args.qsub or args.qsubArray:
            f = open(runEupathdbFile,'w')
            f.write(cmd+"\n")
            f.write("echo \"done!\">%s/%s.done" %(dirs["eupathdb"], db)+ "\n")
            f.close()
            if args.qsub:
                if args.maui:
                    cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runEupathdbFile)
                else:
                    cmdQsub="qsub -cwd -V -N %s -l h_data=16G,time=24:00:00 %s" %(db,runEupathdbFile)
                os.system(cmdQsub)
                write2Log("Job for STEP6 (%s) has been submitted via qsub" %(db),logfns["gLogfile"], args.quiet)
        else:
            os.chdir(dirs["eupathdb"])
            os.system(cmd)
            eupathdbReads=set()
            eupathdbReads=nMicrobialReads(eupathdbFile,readLength,eupathdbFileFiltered)
            nEupathdbReads=len(eupathdbReads)
            write2Log("--identified %s reads mapped %s genomes" %(nEupathdbReads,db) ,logfns["gLogfile"],args.quiet)
            excludeReadsFromFasta(inFasta,eupathdbReads,afterFasta)
            inFasta=afterFasta
            nReadsEP+=nEupathdbReads
    if not args.qsubArray and not args.qsub:
        os.rename(dirs["eupathdb"]+"%s_after_tritryp.fasta" %(basename), intfns["unaccountedReadsFasta"])
        if not args.dev:
            os.remove(intfns["afterVirusFasta"])
        for db in ["ameoba",
               "crypto",
               "giardia",
               "microsporidia",
               "piroplasma",
               "plasmo",
               "toxo",
               "trich"]:
            os.remove(dirs["eupathdb"]+"%s_after_%s.fasta" %(basename,db))
    if not args.qsub and  not args.qsubArray:
        write2File("done!",args.dir+"/step6_microbiomeProfile.done")
    if not args.qsub and  not args.qsubArray and not args.nonReductive:
        write2Log("In toto : %s reads mapped to microbial genomes" %(nReadsBacteria+nReadsVirus+nReadsEP) ,logfns["gLogfile"],args.quiet)
        nTotalReads=nLowQReads+nLowCReads+n_rRNAReads+nlostHumanReads+nRepeatReads+nNCLReads+nReadsImmuneTotal+nReadsBacteria+nReadsVirus+nReadsEP
        write2Log("Summary: The ROP protocol is able to account for %s reads" %(nTotalReads) ,logfns["gLogfile"],args.quiet)
        write2Log("***Unaccounted reads (not explained by ROP) are saved to %s" %(intfns["unaccountedReadsFasta"]) ,logfns["gLogfile"],args.quiet)
        write2Log("***Log file with all the commands used is available here: %s" %(logfns["cmdLogfile"]) ,logfns["gLogfile"],args.quiet)
        tLog=args.dir+"/"+"numberReads_"+basename+".log"
        tLogfile=open(tLog,'w')
        tLogfile.write("sample,totalUnmapped,nLowQReads,nLowCReads,n_rRNAReads,nlostHumanReads,nRepeatReads,nNCLReads,nReadsImmuneTotal,nMicrobialReads")
        tLogfile.write("\n")
        message=basename+","+str(n)+","+str(nLowQReads)+","+str(nLowCReads)+","+str(n_rRNAReads)+","+str(nlostHumanReads)+","+str(nRepeatReads)+","+str(nNCLReads)+","+str(nReadsImmuneTotal)+","+str(nReadsBacteria+nReadsVirus+nReadsEP)
        tLogfile.write(message)
        tLogfile.write("\n")
        tLogfile.close()
else:
    print "6.  Microbiome profiling is skipped."


################################################################################
# Wrap-up

write2Log("The list of the tools used by ROP and the paramemers are provided below" ,logfns["toolsLogfile"],"False")
write2Log("************" ,logfns["toolsLogfile"],"False")
write2Log("**We have used FastQC (version 0.0.13, with the default parameters) downloaded from  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ to filter out low quality reads" ,logfns["toolsLogfile"],"False")
write2Log("**We have used SEQLEAN (seqclean-x86_64, with the default parameters) downloaded from https://sourceforge.net/projects/seqclean/ to filter out low complexity reads , " ,logfns["toolsLogfile"],"False")
write2Log("**We have used Megablast (BLAST+ version 2.2.30, with the following parameters: task	=	megablast,	use_index	=	 true; -outfmt 6 ;-evalue 1e-05; perc_identity	=	100) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to filter out reads mapped to rRNA	repeat sequence	(HSU13369 Human ribosomal DNA	complete repeating	unit) " ,logfns["toolsLogfile"],"False")
write2Log("**We have used Bowtie2(version 2.0.5, with the following parameters: -k 1; -p 8; -f ) downloaded from  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml to identify lost human reads mapped to reference transcriptome and genome (Ensembl GRCh37GRCh37/hg19)" ,logfns["toolsLogfile"],"False")
write2Log("**We have used Megablast (BLAST+ version 2.2.30, with the following options : task=megablast, use_index=true, -outfmt 6 -evalue 1e-05, perc_identity	=	90) downloaded from  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/  to identify lost repeat reads mapped to database of 	repeat	sequences (RepBase20.07) " ,logfns["toolsLogfile"],"False")
write2Log("**We have used CIRI (version 1.2 with the following parameters : -S -I ) downloaded from https://sourceforge.net/projects/ciri/ to identify reads from circRNAs" ,logfns["toolsLogfile"],"False")
write2Log("**We have used IgBLAST (version v 1.4.0 with the following parameters: -germline_db_V;	germline_db_D; -germline_db_J;	-outfmt 7 std qseq sseq; -evalue = 1e-05 )  downloaded from http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/igblast/release/1.4.0/  to identify immune reads spanningBCR/TCR	 receptor gene	 rearrangement	in	 the	variable	domain (V(D)J	recombinations)",logfns["toolsLogfile"],"False")
write2Log("**We have used Megablast (BLAST+ version 2.2.30, with the following parameters: task=megablast, use_index=true; -outfmt 6 ;-evalue 1e-05) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to identify microbial reads mapped onto the microbial genomes (Bacteria, Viruses, and Eukaryotic Pathogens)" ,logfns["toolsLogfile"],"False")
write2Log("**We have used Metaphlan2 (version 2.2.0, with the following parameters: --mpa_pkl ; --bowtie2_exe;  --input_type multifasta; --bowtie2db ; -t reads_map/rel_ab ;  --nproc 8) downloaded from https://bitbucket.org/biobakery/metaphlan2 to obtain taxonomic profile of microbial communities" ,logfns["toolsLogfile"],"False")
write2Log("************" ,logfns["toolsLogfile"],"False")
write2Log("For more information about the paramemers and databases used by ROP, please see the preprint : Dumpster diving in RNA-sequencing to find the source of every last read http://biorxiv.org/content/early/2016/05/13/053041" ,logfns["toolsLogfile"],"False")
write2Log("********************",logfns["gLogfile"],args.quiet)
write2Log("Important: ROP relies on  several open source tools that were developed by other groups. These components are (c) their respective developers and are redistributed with ROP to provide ease-of-use. The list of the tools used by ROP and the parameters/reference databases are provided here: %s " %(logfns["toolsLogfile"]) ,logfns["gLogfile"],args.quiet)

if args.rezip:
    write_gzip(branch_point_file, intfns["afterlostReadsFastaGzip"])
    os.remove(branch_point_file)
