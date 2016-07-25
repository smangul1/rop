import sys
import csv
import os
import argparse
import subprocess
import gzip

#codeDir
codeDir=os.path.dirname(os.path.realpath(__file__))

sys.path.append('%s/tools/biopython/biopython-1.66/' %(codeDir))
import Bio
from Bio import SeqIO # module needed for sequence input






#import pysam




####################################################################
### I/O Functions 
####################################################################

def excludeReadsFromFasta(inFasta,reads,outFasta):

    fasta_sequences = SeqIO.parse(open(inFasta),'fasta')
    with open(outFasta, "w") as f:
        for seq in fasta_sequences:
            name = seq.name
            if name not in reads:
                SeqIO.write([seq], f, "fasta")

def excludeReadsFromFastaGzip(inFasta,reads,outFasta):
    
    fasta_sequences = SeqIO.parse(open(inFasta),'fasta')
    with gzip.open(outFasta, "w") as f:
        for seq in fasta_sequences:
            name = seq.name
            if name not in reads:
                SeqIO.write([seq], f, "fasta")

####################################################################

def bam2fasta(codeDir,inFile,outFile):
    message="Convert bam to fasta"
    cmdConvertBam2Fastq="%s/tools/bamtools convert -in %s -format fasta >%s" %(codeDir,inFile,outFile)
    write2Log(cmdConvertBam2Fastq,cmdLogfile,"False")
    os.system(cmdConvertBam2Fastq)

####################################################################

def bam2fastq(codeDir,inFile,outFile):
    message="Convert bam to fastq"
    write2Log(message,gLogfile,args.quiet)
    cmdConvertBam2Fastq="%s/tools/bamtools convert -in %s -format fastq >%s" %(codeDir,inFile,outFile)
    write2Log(cmdConvertBam2Fastq,cmdLogfile,"False")
    os.system(cmdConvertBam2Fastq)
    if not args.dev:
        os.remove(inFile)

#######################################################################

def write2Log(message,gLog,option):
    logFile=open(gLog,'a')
    if not option:
        print message
    logFile.write(message)
    logFile.write("\n")
    logFile.close()

#######################################################################

def write2File(message,logFile):
    gLogfile=open(logFile,'w')
    gLogfile.write(message)
    gLogfile.close()

#######################################################################

def write_gzip_into_readable(gz_input, output): 
    temp_file = open(output, 'w')
    with gzip.open(gz_input, 'r') as f:
        for line in f:
            temp_file.write(line)
    temp_file.close()

#######################################################################

def nReadsImmune(inFile):
    readsImmune=set()
    with open(inFile,'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
                read=line[1]
                eValue=float(line[12])
                if eValue<1e-05:
                    readsImmune.add(read)
    return readsImmune

#######################################################################

def nMicrobialReads(inFile,readLength,outFile):
    readsMicrobiome=set()
    out=open(outFile,'w')
    with open(inFile,'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
            read=line[0]
            identity=float(line[2])
            alignmentLength=float(line[3])
            eValue=float(line[10])
            if eValue<1e-05 and alignmentLength>=0.8*readLength and identity>=0.9*readLength:
                readsMicrobiome.add(read)
                out.write('\t'.join(line))
                out.write("\n")
    out.close()
    return readsMicrobiome

#######################################################################

def nReadsMetaphlan(inFile):
    readsMetaphlan=set()
    with open(inFile,'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
                read = line[0]
                readsMetaphlan.add(read)
    return readsMetaphlan

#######################################################################
#1:25169311|25169341     1       25169311        25169341        2       1_2_0   6       0.400   n/a     /n/a    SRR1146076.13939638,SRR1146076.25457964,
def nCirrcularReads(inFile):
    reads=set()
    with open(inFile,'r') as f:
        
        next(f)
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
            for kv in line[10].split(","):
                if kv:
                    reads.add(kv)

    return reads

#######################################################################
def write_gzip(inFasta,outFasta):
    
    fasta_sequences = SeqIO.parse(open(inFasta),'fasta')
    with gzip.open(outFasta, "w") as f:
        for seq in fasta_sequences:            
            SeqIO.write([seq], f, "fasta")

####################################################################

print "*********************************************"
print "ROP (version 1.0.4) is a computational protocol aimed to discover the source of all reads, originated from complex RNA molecules, recombinant antibodies and microbial communities. Written by Serghei Mangul (smangul@ucla.edu) and Harry Taegyun Yang (harry2416@gmail.com), University of California, Los Angeles (UCLA). (c) 2016. Released under the terms of the General Public License version 3.0 (GPLv3)"
print ""
print "For more details see:"
print "https://sergheimangul.wordpress.com/rop/"
print "*********************************************"

#######################################################################
### Arguments 
#######################################################################


ap = argparse.ArgumentParser('python rop.py')

necessary_arguments = ap.add_argument_group('Necessary Inputs')
necessary_arguments.add_argument('unmappedReads', help='unmapped Reads in the fastq format')
necessary_arguments.add_argument('dir', help='directory to save results of the analysis')

job_option_arguments = ap.add_argument_group('Job Options')
job_option_arguments.add_argument("--qsub", help="submit qsub jobs on hoffman2 (UCLA) cluster. If planning to use on your cluster contact smangul@ucla.edu", action="store_true")
job_option_arguments.add_argument("--qsubArray", help="prepare qsub scripts to be run later using job array. Working on hoffman2 (UCLA) cluster. If planning to use on your cluster contact smangul@ucla.edu", action="store_true")
job_option_arguments.add_argument("--maui", help = "use this option together with --qsub to submit jobs via  Maui scheduler. Maui is a job scheduler developped by Adaptive Computing. More details are here : https://wiki.calculquebec.ca/w/Maui/en", action = "store_true")


input_option_arguments = ap.add_argument_group('Input Options')
input_option_arguments.add_argument("--b", '-b', help="unmapped reads in bam format", action="store_true")
input_option_arguments.add_argument("--gzip", '-z', help="unmapped reads in .gz format.", action="store_true")
input_option_arguments.add_argument("--skipLowq", help="skip step filtering low quality reads. The input reads need to be in the bam format", action="store_true")
input_option_arguments.add_argument("--skipQC", help="skip entire QC step : filtering  low-quality, low-complexity and rRNA reads. The input reads need to be in the FASTA format", action="store_true")
input_option_arguments.add_argument("--skipPreliminary", '-s', help="skip the preliminary steps including (1) QC and (2) Remapping to human references (lost human reads). The input reads need to be in the FASTA format", action="store_true")



run_only_options = ap.add_argument_group('Run Options - select analysis (multi-selection possible)')
run_only_options.add_argument("--repeat", help = "Run lost repeat profiling ONLY", action = "store_true")
run_only_options.add_argument("--immune", help = "Run BCR/TCR profiling ONLY", action = "store_true")
run_only_options.add_argument("--metaphlan", help = "Run Metaphlan2 ONLY to obtain taxonomic profile of microbial communities", action = "store_true")
run_only_options.add_argument("--circRNA", help = "Run circular RNA profiling ONLY", action="store_true")
run_only_options.add_argument("--microbiome", help = "Run microbime profiling ONLY", action = "store_true")

misc_option_arguments = ap.add_argument_group('Miscellenous Options')

misc_option_arguments.add_argument("--outGz", help = "Beta Mode. Intermediate fasta files are stored as fasta.gz", action = "store_true")
misc_option_arguments.add_argument("--rezip", help = "rezip the fasta files after analysis", action = "store_true")
misc_option_arguments.add_argument("--clean", help = "clean all the intermediate files for maximum space efficiency - use with caution", action = "store_true")
misc_option_arguments.add_argument("--quiet", help = "Suppress progress report and warnings", action = "store_true")
misc_option_arguments.add_argument("--dev", help = "Keep intermediate files", action = "store_true")
misc_option_arguments.add_argument("--nonReductive", help = "non-reductive analysis - Dev mode - Please use with caution", action = "store_true")
misc_option_arguments.add_argument("--f", help = "force option to overwrite the analysis directory(provided as dir option) - Please use with caution  ", action = "store_true")


args = ap.parse_args()





# ONLY OPTION Configuration
# IF none of them are selected: make everything true
if not args.repeat and not args.immune and not args.circRNA and not args.microbiome and not args.metaphlan:
    args.repeat = True
    args.immune = True
    args.circRNA = True
    args.microbiome = True
    args.metaphlan = True
else:
    #It is gonna be non-reductive for now 
    args.nonReductive = True




#######################################################################
### MAIN CODE
#######################################################################

#relative path to absolute path
args.unmappedReads=os.path.abspath(args.unmappedReads)
args.dir=os.path.abspath(args.dir)




#check if args.dir exist
if os.path.exists(args.dir) and not args.f:
    message="ERROR ::::: The directory %s exist. Please choose different directory to save results of the analysis. Alternatively choose --f option to overwrite the results in %s" %(args.dir,args.dir)
    print message
    sys.exit(1)

if os.path.exists(args.dir) and args.f:
    cmd="rm -fr %s>temp 2>temp" %(args.dir)
    os.system(cmd)
    cmd="rm temp"
    os.system(cmd)

#basename
basename=os.path.splitext(os.path.basename(args.unmappedReads))[0].split(".fastq")[0]


#check if the input format is gzip and option --gzip provided
filename, file_extension = os.path.splitext(args.unmappedReads)
if file_extension == ".gz" and not args.gzip:
    print "ERROR ::: --gzip option is not selected with gzip format input. Please try with --gzip option"
    sys.exit(1)








#analysis directories
QCDir=args.dir+"/QC/"




lostHumanDir=args.dir+"/lostHumanReads/"
lostRepeatDir=args.dir+"/lostRepeatSequences/"

antibodyDir=args.dir+"/antibodyProfile/"
bcrDir=antibodyDir+"/BCR/"
tcrDir=antibodyDir+"/TCR/"

NCL_Dir=args.dir+"/NCL/"

ighDir=bcrDir+"/IGH/"
igkDir=bcrDir+"/IGK/"
iglDir=bcrDir+"/IGL/"

tcraDir=tcrDir+"/TCRA/"
tcrbDir=tcrDir+"/TCRB/"
tcrdDir=tcrDir+"/TCRD/"
tcrgDir=tcrDir+"/TCRG/"

microbiomeDir=args.dir+"/microbiomeProfile/"


metaphlanDir = microbiomeDir + "/metaphlan/"


bacteriaDir=args.dir+"/microbiomeProfile/bacteriaProfile/"
virusDir=args.dir+"/microbiomeProfile/viralProfile/"
eupathdbDir=args.dir+"/microbiomeProfile/eukaryoticPathogenProfile/"



if not os.path.exists(QCDir):
    os.makedirs(QCDir)
if not os.path.exists(lostHumanDir):
    os.makedirs(lostHumanDir)
if not os.path.exists(lostRepeatDir):
    os.makedirs(lostRepeatDir)
if not os.path.exists(antibodyDir):
    os.makedirs(antibodyDir)
if not os.path.exists(bcrDir):
    os.makedirs(bcrDir)
if not os.path.exists(tcrDir):
    os.makedirs(tcrDir)
if not os.path.exists(NCL_Dir):
    os.makedirs(NCL_Dir)
for i in [ighDir,igkDir,iglDir,tcraDir,tcrbDir,tcrdDir,tcrgDir]:
    if not os.path.exists(i):
        os.makedirs(i)
if not os.path.exists(microbiomeDir):
    os.makedirs(microbiomeDir)
if not os.path.exists(bacteriaDir):
    os.makedirs(bacteriaDir)
if not os.path.exists(virusDir):
    os.makedirs(virusDir)
if not os.path.exists(eupathdbDir):
    os.makedirs(eupathdbDir)
if not os.path.exists(metaphlanDir):
    os.makedirs(metaphlanDir)

#intermediate files
unmappedFastq=args.dir+"/unmapped_"+basename+".fastq"
lowQFile=QCDir+basename+"_lowQ.fastq"
lowQFileFasta=QCDir+basename+"_lowQ.fa"
lowQCFile=QCDir+basename+"_lowQC.fa"
rRNAFile=QCDir+basename+"_rRNA_blastFormat6.csv"
afterrRNAFasta=QCDir+basename+"_after_rRNA.fasta"

afterlostHumanFasta=lostHumanDir+basename+"_after_rRNA_lostHuman.fasta"
afterlostHumanFastaGzip=lostHumanDir+basename+"_after_rRNA_lostHuman.fasta.gz"
afterNCLFasta=NCL_Dir+basename+"_after_NCL.fasta"
afterImmuneFasta=bcrDir+basename+"_afterImmune.fasta"
afterBacteraFasta=bacteriaDir+basename+"_afterBacteria.fasta"
afterVirusFasta=virusDir+basename+"_afterVirus.fasta"


unaccountedReadsFasta=args.dir+"/"+basename+"_unaccountedReads.fasta"


metaphlan_intermediate_map = metaphlanDir + basename + "_metaphlan.map"
metaphlan_intermediate_bowtie2out = metaphlanDir + basename + "_bowtie2out.txt"
metaphlan_output = metaphlanDir + basename + "_metaphlan_output.tsv"


gBamFile=lostHumanDir+basename+"_genome.sam"
tBamFile=lostHumanDir+basename+"_transcriptome.sam"
repeatFile=lostRepeatDir+basename+"_lostRepeats_blastFormat6.tsv"
afterlostRepeatFasta=lostRepeatDir+basename+"_after_lostRepeat.fasta"
NCL_CIRI_file=NCL_Dir + basename + "_NCL_CIRI_after_bwa.sam"
after_NCL_CIRI_file_prefix = basename + "_circRNA.tsv"
ighFile=ighDir+basename+"_IGH_igblast.tsv"
igkFile=igkDir+basename+"_IGK_igblast.tsv"
iglFile=iglDir+basename+"_IGL_igblast.tsv"
tcraFile=tcraDir+basename+"_TCRA_igblast.tsv"
tcrbFile=tcrbDir+basename+"_TCRB_igblast.tsv"
tcrdFile=tcrdDir+basename+"_TCRD_igblast.tsv"
tcrgFile=tcrgDir+basename+"_TCRG_igblast.tsv"

#log files
logQC=QCDir+basename+"_QC.log"
logrRNA=QCDir + basename + "_rRNA.log"
logHuman=lostHumanDir + basename + "_lostHuman.log"

log_bowtieWG=lostHumanDir + basename + "_bowtieWG.log"
log_bowtieTR=lostHumanDir + basename + "_bowtieTR.log"




logLostRepeat=lostRepeatDir + basename + "_lostRepeat.log"

logNCL=NCL_Dir + basename + "_NCL.log"



bacteriaFile=bacteriaDir+basename+"_bacteria_blastFormat6.csv"
virusFile=virusDir+basename+"_virus_blastFormat6.csv"

bacteriaFileFiltered=bacteriaDir+basename+"_bacteriaFiltered_blastFormat6.csv"
virusFileFiltered=virusDir+basename+"_virusFiltered_blastFormat6.csv"

gLogfile=args.dir+"/"+basename+".log"






cmdLogfile=args.dir+"/"+"dev.log"

toolsLogfile=args.dir+"/"+"tools.log"




#runFiles
runLostHumanFile=lostHumanDir+"/runLostHuman_"+basename+".sh"
runLostRepeatFile=lostRepeatDir+"/runLostRepeat_"+basename+".sh"
runNCL_CIRIfile = NCL_Dir + "/run_NCL_CIRI" + basename + ".sh" 
runIGHFile=ighDir+"/runIGH_"+basename+".sh"
runIGKFile=igkDir+"/runIGK_"+basename+".sh"
runIGLFile=iglDir+"/runIGL_"+basename+".sh"
runTCRAFile=tcraDir+"/runTCRA_"+basename+".sh"
runTCRBFile=tcrbDir+"/runTCRB_"+basename+".sh"
runTCRDFile=tcrdDir+"/runTCRD_"+basename+".sh"
runTCRGFile=tcrgDir+"/runTCRG_"+basename+".sh"
runBacteriaFile=bacteriaDir +"/runBacteria_"+basename+".sh"
runVirusFile=virusDir +"/runVirus_"+basename+".sh"
run_metaphlan_file = metaphlanDir + "/run_metaphlan_" + basename + ".sh"

os.chdir(args.dir)



logIGH=ighDir + basename + "_igh.log"
logIGL=iglDir + basename + "_igl.log"
logIGK=igkDir + basename + "_igk.log"

logTCRA=tcraDir + basename + "_tcra.log"
logTCRB=tcrbDir + basename + "_tcrb.log"
logTCRG=tcrgDir + basename + "_tcrg.log"
logTCRD=tcrdDir + basename + "_tcrd.log"

logMetaphlan=metaphlanDir + basename + "_metaphlan.log"

logBacteria=bacteriaDir + basename + "_bacteria.log"
logVirus=virusDir + basename + "_virus.log"
logEukaryotes=eupathdbDir + basename + "_eukaryotes.log"



#######################################################################################################################################


readLength=0
n=0

### TODO : FIX WHEN INPUT IS BAM 
# The current version assumes the input is fastq or fasta 



# if the input is bam
# change the input file so the num_unmapped can be assessed
"""
unmapped file is used just to handle the number of reads

If input is bam: 
    if lowq bam file - then it is going to be converted into fasta and use it as unmapped_file
    else (if it is raw unmapped bam) - converted into fastq and used it as unmapped_file 
else (if fastq): 
    unmapped_file = input 
"""


if args.b:
    if args.skipLowq:
        bam2fasta(codeDir,args.unmappedReads,lowQFileFasta)
        unmapped_file = lowQFileFasta
    else:
        bam2fastq(codeDir,args.unmappedReads,unmappedFastq)
        unmapped_file = unmappedFastq
elif args.gzip and not args.skipPreliminary and not args.skipQC and not args.skipLowq:
    write2Log("--gzip option is selected. The input is in gzip format. It will be decompressed in the analysis dir.", gLogfile, args.quiet)
    unmappedFastq=args.dir+"/"+basename+".fastq"
    write_gzip_into_readable(args.unmappedReads, unmappedFastq)
    unmapped_file = unmappedFastq
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
write2Log(message,gLogfile,args.quiet)


nLowQReads=0
nLowCReads=0
n_rRNAReads=0
nlostHumanReads=0


if args.skipPreliminary:
    # afterrRNAFasta=args.unmappedReads
    write2Log("1. Quality Control is skipped",gLogfile,args.quiet)
    write2Log("2. Remapping to human references is skipped",gLogfile,args.quiet)
    
    filename, file_extension = os.path.splitext(args.unmappedReads)
    print file_extension
    if file_extension!=".fa" and file_extension!=".fasta" and not args.gzip :
        write2Log("ERROR ::: --skipPreliminary option is selected. Reads needs to be in FASTA format",gLogfile,args.quiet)
        sys.exit(1)
    elif file_extension == ".gz":
        temp=args.unmappedReads.split(".gz")[0]
        filename2, file_extension2 = os.path.splitext(temp)
        if file_extension2!=".fa" and file_extension2!=".fasta":
            write2Log("ERROR ::: --skipPreliminary option is selected. Reads needs to be in FASTA format",gLogfile,args.quiet)
            sys.exit(1)
        
        if not args.gzip:
            write2Log("ERROR ::: --gzip option is not selected with gzip format input. Please try with --gzip option", gLogfile, args.quiet)
    elif file_extension == ".gz" and args.gzip: 
        write2Log("SkipPreliminary option is selected. The input is in gzip format. It will be decompressed in the analysis dir.", gLogfile, args.quiet)


elif args.skipQC:
    filename, file_extension = os.path.splitext(args.unmappedReads)
    if file_extension!=".fa" and file_extension!=".fasta":
        write2Log("ERROR ::: --skipQC option is selected. Reads needs to be in FASTA format",gLogfile,args.quiet)
        sys.exit(1)
    # afterrRNAFasta = args.unmappedReads
    write2Log("1. Quality Control is skipped",gLogfile,args.quiet)
else:

    if not args.b and not args.gzip:
        unmappedFastq=args.unmappedReads







    valid = 'ACTG'

    if not args.skipLowq:
        #lowQ
        write2Log("1. Quality Control...",gLogfile,args.quiet)
        
        fastafile=open(lowQFileFasta, 'w')
        readLength=0
        nLowQReads=0
        nAfterLowQReads=0
        
        
        for record in SeqIO.parse(unmappedFastq, "fastq"):
            
            
            
            
            
            readLength=len(record) #assumes the same length, will not work for Ion Torrent or Pac Bio
            
            j=record.letter_annotations["phred_quality"]
            
            prc=len([i for i in j if i>=20])/float(len(j))
            if prc>0.75 and all(i in valid for i in record.seq):
                fastafile.write(str(">" + record.name) + "\n")
                fastafile.write(str(record.seq) + "\n")
                nAfterLowQReads+=1
        
            



        if args.b:
            os.remove(unmappedFastq)








        fastafile.close()
        nLowQReads=n-nAfterLowQReads
        write2Log("--filtered %s low quality reads" % (nLowQReads) ,gLogfile,args.quiet)
    else:
        write2Log("skipLowq option is selected - low quality filtering step is skipped" ,gLogfile, args.quiet)


    



    os.chdir(QCDir)
    #lowC
    cmd="export PATH=$PATH:%s/tools/seqclean-x86_64/bin" %(codeDir)
    os.system(cmd)

    cmd=codeDir+"/tools/seqclean-x86_64/seqclean %s -l 50 -M -o %s 2>>%s" %(lowQFileFasta, lowQCFile,logQC)
    write2Log(cmd,cmdLogfile,"False")
    os.system(cmd)

    cmd = "rm -rf %s/cleaning_1/ >temp 2 >temp; rm -f %s/*.cln >temp 2 >temp; rm -f %s/*.cidx>temp 2 >temp; rm -f %s/*.sort>temp 2 >temp" % (QCDir,QCDir,QCDir,QCDir)
    os.system(cmd)
    proc = subprocess.Popen(["grep trashed %s | awk -F \":\" '{print $2}'" %(logQC) ], stdout=subprocess.PIPE, shell=True)
    (nLowCReadsTemp, err) = proc.communicate()
    nLowCReads=int(nLowCReadsTemp.rstrip().strip())
    write2Log("--filtered %s low complexity reads (e.g. ACACACAC...)" %(nLowCReads) ,gLogfile,args.quiet)







    #rRNA
    cmd="%s/tools/blastn -task megablast -index_name %s/db/rRNA/rRNA -use_index true -query %s -db %s/db/rRNA/rRNA  -outfmt 6 -evalue 1e-05  >%s" %(codeDir,codeDir,lowQCFile,codeDir,rRNAFile)
    write2Log(cmd,cmdLogfile,"False")
    os.system(cmd)

    n_rRNATotal=0
    rRNAReads = set()
    with open(rRNAFile,'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
            n_rRNATotal+=1
            element=line[0]
            identity=float(line[2])
            alignmentLength=float(line[3])
            eValue=float(line[10])
            if eValue<1e-05 and alignmentLength==readLength and identity>=0.94*readLength:
                rRNAReads.add(element)

    excludeReadsFromFasta(lowQCFile,rRNAReads,afterrRNAFasta)
    n_rRNAReads=len(rRNAReads)
    write2Log("--filtered %s rRNA reads" %(n_rRNAReads) ,gLogfile,args.quiet)
    write2Log("In toto : %s reads failed QC and are filtered out" %(nLowQReads+nLowCReads+n_rRNAReads) ,gLogfile,args.quiet)
    write2File("done!",args.dir+"/step1_QC.done")

    message="Number of entries in %s is %s" %(rRNAFile,n_rRNATotal)
    write2Log(message,cmdLogfile,"False")

    if not args.dev:
        os.remove(lowQCFile)
        os.remove(lowQFileFasta)
        os.remove(rRNAFile)






#######################################################################################################################################
#2. Remaping to human references...
if not args.skipPreliminary:
    write2Log("2. Remapping to human references...",cmdLogfile,"False")
    write2Log("2. Remapping to human references...",gLogfile,args.quiet)
    # If input is afterQC fasta.gz
    if args.skipQC and args.outGz:
        
        
        write_gzip_into_readable(args.unmappedReads, afterrRNAFasta)
    elif args.skipQC and not args.outGz:
        afterrRNAFasta = args.unmappedReads
    cmdGenome="%s/tools/bowtie2 -k 1 -p 8 -f -x %s/db/bowtie2Index/genome -U %s 2>>%s | %s/tools/samtools view -SF4 -   >%s" %(codeDir,codeDir, afterrRNAFasta,log_bowtieWG,codeDir,gBamFile)

    #transcriptome
    cmdTranscriptome="%s/tools/bowtie2  -k 1 -f -p 8 -x %s/db/bowtie2Index/hg19KnownGene.exon_polya200 -U %s 2>>%s | %s/tools/samtools view -SF4 -  >  %s " %(codeDir,codeDir, afterrRNAFasta,log_bowtieTR, codeDir,tBamFile)
    write2Log(cmdGenome,cmdLogfile,"False")
    write2Log(cmdTranscriptome,cmdLogfile,"False")




    os.system(cmdTranscriptome)
    os.system(cmdGenome)

    lostHumanReads = set()
    lostHumanReads0= set()
    lostHumanReads1= set()
    lostHumanReads2= set()




    with open(tBamFile,'r') as f:
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


    with open(gBamFile,'r') as f:
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



    write2Log("Unmapped reads mapped to human genome and/or transcriptome (using bowtie2) are categorized as lost human reads and are excluded from the further analysis. This includes :  %s reads with 0 mismatches; %s reads with 1 mismatch; %s reads with 2 mismatches" %(len(lostHumanReads0),len(lostHumanReads1),len(lostHumanReads2)),logHuman,"False")
    write2Log("Complete list of lost human reads is available from sam files: %s,%s" %(gBamFile,tBamFile),logHuman,"False")


    if args.outGz:
        excludeReadsFromFastaGzip(afterrRNAFasta,lostHumanReads,afterlostHumanFastaGzip)
    else:
        excludeReadsFromFasta(afterrRNAFasta,lostHumanReads,afterlostHumanFasta)
    nlostHumanReads=len(lostHumanReads)
    write2Log("--identified %s lost human reads from unmapped reads. Among those: %s reads with 0 mismatches; %s reads with 1 mismatch; %s reads with 2 mismatches" %(len(lostHumanReads),len(lostHumanReads0),len(lostHumanReads1),len(lostHumanReads2)), gLogfile, args.quiet)
    write2Log("***Note: Complete list of lost human reads is available from sam files: %s,%s" %(gBamFile,tBamFile), gLogfile, args.quiet)
    write2File("done!",args.dir+"/step2_lostHumanReads.done")

    if not args.dev:
        if not args.skipQC and not args.skipPreliminary:
            os.remove(afterrRNAFasta)
    if args.clean:
        write2Log("Clean mode selected - removing analysis sam files", gLogfile, args.quiet)
        os.remove(gBamFile)
        os.remove(tBamFile)



### TODO
else:
    if args.outGz:
        write_gzip_into_readable(args.unmappedReads, afterlostHumanFasta)

    else:
        afterlostHumanFasta = args.unmappedReads

if args.nonReductive or args.qsub or args.qsubArray:
    branch_point_file = afterlostHumanFasta
    
    write2Log("*********************************",gLogfile, args.quiet)
    write2Log("Non-substractive mode is selected : Low quality, low complexity, rRNA reads and lost human reads are filtered out. Resulting high quality non-human reads are provided as input  for STEP3-STEP6",gLogfile, args.quiet)
    write2Log("*********************************",gLogfile, args.quiet)

    

#######################################################################################################################################
#3. Maping to repeat sequences...
if args.repeat:
    write2Log("3. Maping to repeat sequences...",cmdLogfile,"False")
    write2Log("3. Maping to repeat sequences...",gLogfile,args.quiet)

    #TO DO : make all fasta ->gzip
    #gzip -dc %s | , query -
    # CHANGED 
    # cmd="%s/tools/blastn -task megablast -index_name %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa -use_index true -query %s -db %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa  -outfmt 6 -evalue 1e-05  > %s" %(codeDir, codeDir, afterlostHumanFasta, codeDir, repeatFile)
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = afterlostHumanFasta
    cmd="%s/tools/blastn -task megablast -index_name %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa -use_index true -query %s -db %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa  -outfmt 6 -evalue 1e-05  > %s 2>%s" %(codeDir, codeDir, input_file, codeDir, repeatFile,logLostRepeat)




    if args.qsub or args.qsubArray:
        f = open(runLostRepeatFile,'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_lostRepeat.done \n" %(lostRepeatDir,basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runLostRepeatFile)
            else:
                cmdQsub="qsub -cwd -V -N lostRepeat -l h_data=16G,time=10:00:00 %s" %(runLostRepeatFile)
            os.system(cmdQsub)
            write2Log("Job for STEP3 has been submitted via qsub",gLogfile, args.quiet)


    else:
        os.system(cmd)
        write2Log(cmd,cmdLogfile,"False")



    if not args.qsub and not args.qsubArray:

        lostRepeatReads = set()
        
        with open(repeatFile,'r') as f:
            reader=csv.reader(f,delimiter='\t')
            for line in reader:
                element=line[0]
                identity=float(line[2])
                alignmentLength=float(line[3])
                eValue=float(line[10])
                if eValue<1e-05 and alignmentLength>=0.8*readLength and identity>=0.9*readLength:
                    lostRepeatReads.add(element)

        nRepeatReads=len(lostRepeatReads)
        write2Log("-- Identified %s lost repeat sequences from unmapped reads" %(nRepeatReads) ,gLogfile,args.quiet)
        write2Log("***Note : Repeat sequences classification into classes (e.g. LINE) and families (e.g. Alu) will be available in next release" ,gLogfile,args.quiet)
        excludeReadsFromFasta(afterlostHumanFasta,lostRepeatReads,afterlostRepeatFasta)
        write2File("done!",args.dir+"/step3_lostRepeatSequences.done")
        if not args.dev:
            if not args.skipPreliminary:
                os.remove(afterlostHumanFasta)


        write2Log("Lost repeat reads are mapped to the repeat sequences (using megablast)",logLostRepeat,"False")
        write2Log("Complete list of lost repeat reads is available from tsv file: %s" %(repeatFile),logLostRepeat,"False")

else:
    print "3. Maping to repeat sequences is skipped."

#######################################################################################################################################
#3. Non-co-linear RNA profiling
write2Log("4. Non-co-linear RNA profiling",cmdLogfile,"False")
write2Log("4. Non-co-linear RNA profiling",gLogfile,args.quiet)
write2Log("***Note : Trans-spicing and gene fusions  are currently not supported, but will be in the next release.",gLogfile,args.quiet)

os.chdir(NCL_Dir)
NCL_reads=set()
nNCLReads=0



if args.circRNA:
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = afterlostRepeatFasta
    cmd="%s/tools/bwa mem -T -S %s/db/BWAIndex/genome.fa %s > %s 2>%s \n" %(codeDir,codeDir,input_file,NCL_CIRI_file,logNCL)
    cmd = cmd + "perl %s/tools/CIRI_v1.2.pl -S -I %s -O %s -F %s/db/BWAIndex/genome.fa 1>>%s 2>>%s" %(codeDir,NCL_CIRI_file,after_NCL_CIRI_file_prefix,codeDir,logNCL,logNCL)


    if args.qsub or args.qsubArray:
        f = open(runNCL_CIRIfile,'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_NCL_CIRI.done \n" %(NCL_Dir,basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runNCL_CIRIfile)
            else:
                cmdQsub="qsub -cwd -V -N NCL_CIRI -l h_data=8G,time=10:00:00 %s" %(runNCL_CIRIfile)
            os.system(cmdQsub)
            write2Log("Job for STEP4 has been submitted via qsub",gLogfile, args.quiet)
    else:
        os.chdir(NCL_Dir)
        os.system(cmd)
        NCL_reads=nCirrcularReads(after_NCL_CIRI_file_prefix)
        nReadsNCL=len(NCL_reads)
        write2Log("--identified %s reads from circRNA" %(nReadsNCL) ,gLogfile,args.quiet)
        write2Log("***Note: circRNAs detected by CIRI are available here: %s" %(after_NCL_CIRI_file_prefix) ,gLogfile,args.quiet)
        excludeReadsFromFasta(input_file,NCL_reads,afterNCLFasta)
        nNCLReads=len(NCL_reads)
        write2File("done!",args.dir+"/step4_NCL.done")
        if not args.dev:
            os.remove(input_file)




else:
    print "4. Non-co-linear RNA profiling is skipped."

#######################################################################################################################################
#5. T and B lymphocytes profiling

if args.immune:
    immuneReads=set()

    write2Log("5a. B lymphocytes profiling...",cmdLogfile,"False")
    write2Log("5a. B lymphocytes profiling...",gLogfile,args.quiet)


    #IGH-------
    os.chdir(ighDir)
    cmd="ln -s %s//db/antibody/internal_data/ ./" %(codeDir)
    os.system(cmd)

    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = afterNCLFasta

    write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of immunoglobulin heavy locus (IGH)",logIGH,"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logIGH,"False")


    cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/IGHV.fa -germline_db_D %s/db/antibody/IGHD.fa  -germline_db_J %s/db/antibody/IGHJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05  2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,input_file, logIGH,ighFile)
    write2Log(cmd,cmdLogfile,"False")


    if args.qsub or args.qsubArray:
        f = open(runIGHFile,'w')
        f.write("ln -s %s//db/antibody/internal_data/ ./ \n" %(codeDir))
        f.write(cmd+"\n")
        f.write("echo \"done!\"> %s/%s_igh.done \n" %(ighDir,basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runIGHFile)
            else:
                cmdQsub="qsub -cwd -V -N igh -l h_data=16G,time=24:00:00 %s" %(runIGHFile)
            os.system(cmdQsub)
            write2Log("Job for STEP5a(IGK) has been submitted via qsub",gLogfile, args.quiet)
    else:
        os.chdir(ighDir)
        os.system(cmd)
        immuneReadsIGH=nReadsImmune(ighFile)
        nReadsImmuneIGH=len(immuneReadsIGH)
        write2Log("--identified %s reads mapped to immunoglobulin heavy (IGH) locus" %(nReadsImmuneIGH) ,gLogfile,args.quiet)


    #IGK---------
    os.chdir(igkDir)
    cmd="ln -s %s//db/antibody/internal_data/ ./" %(codeDir)
    os.system(cmd)

    write2Log("IgBlast	was	used	to	identify	reads	spanning	VJ recombinations of immunoglobulin kappa locus (IGK)",logIGK,"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logIGK,"False")

    

    cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/IGKV.fa -germline_db_D %s/db/antibody/IGHD.fa  -germline_db_J %s/db/antibody/IGKJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05  2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,input_file, logIGK,igkFile)


    write2Log(cmd,cmdLogfile,"False")

    if args.qsub or args.qsubArray:
        f = open(runIGKFile,'w')
        f.write("ln -s %s//db/antibody/internal_data/ ./ \n" %(codeDir))
        f.write(cmd+"\n")
        f.write("echo \"done!\"> %s/%s_igk.done \n" %(igkDir,basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runIGKFile)
            else:
                cmdQsub="qsub -cwd -V -N igk -l h_data=16G,time=24:00:00 %s" %(runIGKFile)
            os.system(cmdQsub)
            write2Log("Job for STEP5a(IGK) has been submitted via qsub",gLogfile, args.quiet)
    else:
        os.chdir(igkDir)
        os.system(cmd)
        immuneReadsIGK=nReadsImmune(igkFile)
        nReadsImmuneIGK=len(immuneReadsIGK)
        write2Log("--identified %s reads mapped to immunoglobulin kappa (IGK) locus " %(nReadsImmuneIGK) ,gLogfile,args.quiet)
                
                
    #IGL------------
    os.chdir(iglDir)
    cmd="ln -s %s//db/antibody/internal_data/ ./" %(codeDir)
    os.system(cmd)

    write2Log("IgBlast	was	used	to	identify	reads	spanning	VJ recombinations of immunoglobulin lambda locus (IGL)",logIGL,"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logIGL,"False")




    cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/IGLV.fa -germline_db_D %s/db/antibody/IGHD.fa  -germline_db_J %s/db/antibody/IGLJ.fa -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,input_file,logIGL,iglFile)
    write2Log(cmd,cmdLogfile,"False")

    if args.qsub or args.qsubArray:
        f = open(runIGLFile,'w')
        f.write("ln -s %s//db/antibody/internal_data/ ./ \n" %(codeDir))
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_igl.done \n" %(iglDir,basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runIGLFile)
            else:
                cmdQsub="qsub -cwd -V -N igl -l h_data=16G,time=24:00:00 %s" %(runIGLFile)
            os.system(cmdQsub)
            write2Log("Job for STEP5a(IGL) has been submitted via qsub",gLogfile, args.quiet)
    else:
        os.chdir(iglDir)
        os.system(cmd)
        immuneReadsIGL=nReadsImmune(iglFile)
        nReadsImmuneIGL=len(immuneReadsIGL)
        write2Log("--identified %s reads mapped to immunoglobulin lambda (IGL) locus" %(nReadsImmuneIGL) ,gLogfile,args.quiet)


    ##################
    ##################
    write2Log("5b. T lymphocytes profiling...",cmdLogfile,"False")
    write2Log("5b. T lymphocytes profiling...",gLogfile,args.quiet)

    #TCRA-----------------
    os.chdir(tcraDir)
    cmd="ln -s %s//db/antibody/internal_data/ ./" %(codeDir)
    os.system(cmd)

    write2Log("IgBlast	was	used	to	identify	reads	spanning	VJ recombinations of T cell receptor alpha locus (TCRA)",logTCRA,"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logTCRA,"False")
    
    cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/TRAV.fa -germline_db_D %s/db/antibody/TRBD.fa  -germline_db_J %s/db/antibody/TRAJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,input_file,logTCRA,tcraFile)
    write2Log(cmd,cmdLogfile,"False")
    if args.qsub or args.qsubArray:
        f = open(runTCRAFile,'w')
        f.write("ln -s %s//db/antibody/internal_data/ ./ \n" %(codeDir))
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_tcra.done \n"%(tcraDir,basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runTCRAFile)
            else:
                cmdQsub="qsub -cwd -V -N tcra -l h_data=16G,time=24:00:00 %s" %(runTCRAFile)
            os.system(cmdQsub)
            write2Log("Job for STEP5b(TCRA) has been submitted via qsub",gLogfile, args.quiet)
    else:
        os.chdir(tcraDir)
        os.system(cmd)
        immuneReadsTCRA=nReadsImmune(tcraFile)
        nReadsImmuneTCRA=len(immuneReadsTCRA)
        write2Log("--identified %s reads mapped to T cell receptor alpha (TCRA) locus" %(nReadsImmuneTCRA) ,gLogfile,args.quiet)
                
                

    #TCRB--------------
    os.chdir(tcrbDir)
    cmd="ln -s %s//db/antibody/internal_data/ ./" %(codeDir)
    os.system(cmd)

    write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of T cell receptor beta locus (TCRB)",logTCRB,"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logTCRB,"False")

    cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/TRBV.fa -germline_db_D %s/db/antibody/TRBD.fa  -germline_db_J %s/db/antibody/TRBJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,input_file,logTCRB,tcrbFile)
    write2Log(cmd,cmdLogfile,"False")
    if args.qsub or args.qsubArray:
        f = open(runTCRBFile,'w')
        f.write("ln -s %s//db/antibody/internal_data/ ./ \n" %(codeDir))
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_tcrb.done \n"%(tcrbDir,basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runTCRBFile)
            else:
                cmdQsub="qsub -cwd -V -N tcrb -l h_data=16G,time=24:00:00 %s" %(runTCRBFile)
            os.system(cmdQsub)
            write2Log("Job for STEP5b(TCRB) has been submitted via qsub",gLogfile, args.quiet)

    else:
        os.chdir(tcrbDir)
        os.system(cmd)
        immuneReadsTCRB=nReadsImmune(tcrbFile)
        nReadsImmuneTCRB=len(immuneReadsTCRB)
        write2Log("--identified %s reads mapped to T cell receptor beta (TCRB) locus" %(nReadsImmuneTCRB) ,gLogfile,args.quiet)


                

    #TCRD----------------
    os.chdir(tcrdDir)
    cmd="ln -s %s//db/antibody/internal_data/ ./" %(codeDir)
    os.system(cmd)

    write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of T cell receptor delta locus (TCRD)",logTCRD,"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logTCRD,"False")

    cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/TRDV.fa -germline_db_D %s/db/antibody/TRBD.fa  -germline_db_J %s/db/antibody/TRDJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,input_file,logTCRD,tcrdFile)
    write2Log(cmd,cmdLogfile,"False")
    if args.qsub or args.qsubArray:
        f = open(runTCRDFile,'w')
        f.write("ln -s %s//db/antibody/internal_data/ ./ \n" %(codeDir))
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_tcrd.done \n" %(tcrdDir, basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runTCRDFile)
            else:
                cmdQsub="qsub -cwd -V -N tcrd -l h_data=16G,time=24:00:00 %s" %(runTCRDFile)
            os.system(cmdQsub)
            write2Log("Job for STEP5b(TCRD) has been submitted via qsub",gLogfile, args.quiet)

    else:
        os.chdir(tcrdDir)
        os.system(cmd)
        immuneReadsTCRD=nReadsImmune(tcrdFile)
        nReadsImmuneTCRD=len(immuneReadsTCRD)
        write2Log("--identified %s reads mapped to T cell receptor delta (TCRD) locus" %(nReadsImmuneTCRD) ,gLogfile,args.quiet)


                
    #TCRG---------------------
    os.chdir(tcrgDir)
    cmd="ln -s %s//db/antibody/internal_data/ ./" %(codeDir)
    os.system(cmd)

    write2Log("IgBlast	was	used	to	identify	reads	spanning	VDJ recombinations of T cell receptor gamma locus (TCRG)",logTCRG,"False")
    write2Log("More details about IgBlast format are here : https://github.com/smangul1/rop/wiki/ROP-output-details ",logTCRG,"False")


    cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/TRGV.fa -germline_db_D %s/db/antibody/TRBD.fa  -germline_db_J %s/db/antibody/TRGJ.fa  -query %s -outfmt '7 std qseq sseq' -evalue 1e-05 2>>%s  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,input_file,logTCRG,tcrgFile)
    write2Log(cmd,cmdLogfile,"False")

    if args.qsub or args.qsubArray:
        f = open(runTCRGFile,'w')
        f.write("ln -s %s//db/antibody/internal_data/ ./ \n" %(codeDir))
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_tcrg.done \n" %(tcrgDir,basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runTCRGFile)
            else:
                cmdQsub="qsub -cwd -V -N tcrg -l h_data=16G,time=24:00:00 %s" %(runTCRGFile)
            os.system(cmdQsub)
            write2Log("Job for STEP5b(TCRG) has been submitted via qsub",gLogfile, args.quiet)

    else:
        os.chdir(tcrgDir)
        os.system(cmd)
        immuneReadsTCRG=nReadsImmune(tcrgFile)
        nReadsImmuneTCRG=len(immuneReadsTCRG)
        write2Log("--identified %s reads mapped to T cell receptor gamma locus (TCRG) locus" %(nReadsImmuneTCRG) ,gLogfile,args.quiet)

    nReadsImmuneTotal=0
    if not args.qsub and not args.qsubArray:
        nReadsImmuneTotal=nReadsImmuneIGH+nReadsImmuneIGL+nReadsImmuneIGK+nReadsImmuneTCRA+nReadsImmuneTCRB+nReadsImmuneTCRD+nReadsImmuneTCRG
        write2Log("In toto : %s reads mapped to antibody repertoire loci" %(nReadsImmuneTotal) ,gLogfile,args.quiet)
        write2Log("***Note : Combinatorial diversity of the antibody repertoire (recombinations of the of VJ gene segments)  will be available in the next release.",gLogfile,args.quiet)
        write2File("done!",args.dir+"/step5_antibodyProfile.done")
        immuneReads=set().union(immuneReadsTCRA,immuneReadsTCRB,immuneReadsTCRD,immuneReadsTCRG)
        excludeReadsFromFasta(input_file,immuneReads,afterImmuneFasta)
        if not args.dev:
            if args.circRNA:           
                os.remove(afterNCLFasta)
else:
    print "5a. B lymphocytes profiling is skipped."
    print "5b. T lymphocytes profiling is skipped."
#######################################################################################################################################
# 5. Metaphlan



if args.microbiome:
    write2Log("***Extra step.  Metaphlan profiling...",cmdLogfile,"False")
    write2Log("***Extra step.  Metaphlan profiling...",gLogfile,args.quiet)
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = afterImmuneFasta
    cmd = "python %s/tools/metaphlan2.py %s %s --mpa_pkl %s/db/metaphlan/mpa_v20_m200.pkl --bowtie2_exe %s/tools/bowtie2 --input_type multifasta --bowtie2db %s/db/metaphlan/mpa_v20_m200 -t reads_map --nproc 8 --bowtie2out %s>>%s 2>>%s" % (codeDir, input_file, metaphlan_intermediate_map, codeDir, codeDir, codeDir, metaphlan_intermediate_bowtie2out,logMetaphlan,logMetaphlan)
    cmd = cmd + "\n" + "python %s/tools/metaphlan2.py --mpa_pkl %s/db/metaphlan/mpa_v20_m200.pkl --bowtie2_exe %s/tools/bowtie2 --input_type bowtie2out %s -t rel_ab > %s 2>>%s" %(codeDir, codeDir, codeDir, metaphlan_intermediate_bowtie2out, metaphlan_output,logMetaphlan)
    write2Log(cmd,cmdLogfile,"False")

    if args.qsub or args.qsubArray:
        f= open(run_metaphlan_file, 'w')
        f.write(cmd + "\n")
        f.write("echo \"done!\" > %s/%s_metaphlan.done \n" % (metaphlanDir, basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(run_metaphlan_file)
            else:
                cmdQsub="qsub -cwd -V -N metaphlan -l h_data=16G,time=24:00:00 %s" %(run_metaphlan_file)
            os.system(cmdQsub)
            write2Log("Job for Metaphlan has been submitted via qsub",gLogfile, args.quiet)

    else:
        os.chdir(metaphlanDir)
        os.system(cmd)
        write2Log("***Microbiome profiling by Metaphlan2: taxonomic profile of microbial communities detected by Metaphlan2 is available here: %s" %(metaphlanDir) ,gLogfile,args.quiet)


else:
    print "Extra step.  Metaphlan profiling is skipped."


#######################################################################################################################################
# 6. Microbiome profiling...
if args.microbiome:
    write2Log("6.  Microbiome profiling...",cmdLogfile,"False")
    write2Log("6.  Microbiome profiling...",gLogfile,args.quiet)
    
    
    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = afterImmuneFasta






    #bacteria -------------------------------
    os.chdir(bacteriaDir)
    readLength=0
    if not args.outGz:
        fastafile = open(input_file, "rU")
    else:
        fastafile = gzip.open(input_file, "rU")

    for record in SeqIO.parse(fastafile,"fasta"):
        readLength=len(record) #assumes the same length, will not work for Ion Torrent or Pac Bio
        break;


    write2Log("Megablast was used to map the reads onto the	bacterial reference	genomes ",logBacteria,"False")
    write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides aligned are considered microbial reads and are saved into %s " %(bacteriaFileFiltered), logBacteria,"False")
    write2Log("---------------",logBacteria,"False")

    cmd="%s/tools/blastn -task megablast -index_name %s/db/bacteria/bacteria -use_index true -query %s -db %s/db/bacteria/bacteria  -outfmt 6 -evalue 1e-05  >%s 2>>%s" %(codeDir,codeDir,input_file,codeDir,bacteriaFile,logBacteria)
    write2Log(cmd,cmdLogfile,"False")

    if args.qsub or args.qsubArray:
        f = open(runBacteriaFile,'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_bacteria.done \n"%(bacteriaDir, basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runBacteriaFile)
            else:
                cmdQsub="qsub -cwd -V -N bacteria -l h_data=16G,time=24:00:00 %s" %(runBacteriaFile)
            os.system(cmdQsub)
            write2Log("Job for STEP6 (bacteriaProfile) has been submitted via qsub",gLogfile, args.quiet)

    else:
        os.chdir(bacteriaDir)
        os.system(cmd)
        bacteriaReads=nMicrobialReads(bacteriaFile,readLength,bacteriaFileFiltered)
        nReadsBacteria=len(bacteriaReads)
        write2Log("--identified %s reads mapped bacterial genomes" %(nReadsBacteria) ,gLogfile,args.quiet)
        excludeReadsFromFasta(input_file,bacteriaReads,afterBacteraFasta)
    



    #virus-----------------------------------------------------
    os.chdir(virusDir)

    if args.nonReductive or args.qsub or args.qsubArray:
        input_file = branch_point_file
    else:
        input_file = afterBacteraFasta
        if not args.dev:
            os.remove(afterImmuneFasta)


    write2Log("Megablast was used to map the reads onto the	viral reference	genomes ",logVirus,"False")
    write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides aligned are considered microbial reads and are saved into %s" %(virusFileFiltered),logVirus,"False")
    write2Log("----------------------",logVirus,"False")



    cmd="%s/tools/blastn -task megablast -index_name %s/db/virus/viruses -use_index true -query %s -db %s/db/virus/viruses  -outfmt 6 -evalue 1e-05  >%s 2>>%s" %(codeDir,codeDir,input_file,codeDir,virusFile,logVirus)
    write2Log(cmd,cmdLogfile,"False")

    if args.qsub or args.qsubArray:
        f = open(runVirusFile,'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_virus.done \n"%(bacteriaDir, basename))
        f.close()
        if args.qsub:
            if args.maui:
                cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runBacteriaFile)
            else:
                cmdQsub="qsub -cwd -V -N virus -l h_data=16G,time=24:00:00 %s" %(runBacteriaFile)
            os.system(cmdQsub)
            write2Log("Job for STEP6 (viralProfile) has been submitted via qsub",gLogfile, args.quiet)

    else:
        os.system(cmd)
        virusReads=nMicrobialReads(virusFile,readLength,virusFileFiltered)
        nReadsVirus=len(virusReads)
        write2Log("--identified %s reads mapped viral genomes" %(nReadsVirus) ,gLogfile,args.quiet)
        excludeReadsFromFasta(input_file,virusReads,afterVirusFasta)




    #eukaryotic pathogens-----------------------
    os.chdir(eupathdbDir)

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
        inFasta = afterVirusFasta
        if not args.dev:
            os.remove(afterBacteraFasta)


    nReadsEP=0

    for db in dbList:
        afterFasta=eupathdbDir+"%s_after_%s.fasta" %(basename,db)
        eupathdbFile=eupathdbDir+basename+"_"+db+"_blastFormat6.csv"
        eupathdbFileFiltered=eupathdbDir+basename+"_"+db+"Filtered_blastFormat6.csv"
        runEupathdbFile=eupathdbDir+"/run_"+basename+"_"+db+".sh"

        write2Log("Megablast was used to map the reads onto the	reference genomes of Eukaryotes (http://eupathdb.org/eupathdb/)",logEukaryotes,"False")
        write2Log("Reads with 0.9 identity and more than 0.8 of nucleotides aligned are considered microbial reads and are saved into %s" %(afterFasta), logEukaryotes,"False")
        write2Log("-------------",logEukaryotes,"False")
    
        cmd="%s/tools/blastn -task megablast -index_name %s/db/eupathdb/%s -use_index true -query %s -db %s/db/eupathdb/%s  -outfmt 6 -evalue 1e-05  >%s 2>>%s" %(codeDir,codeDir,db,inFasta,codeDir,db,eupathdbFile,logEukaryotes)
        write2Log(cmd,cmdLogfile,"False")
        
        if args.qsub or args.qsubArray:
            f = open(runEupathdbFile,'w')
            f.write(cmd+"\n")
            f.write("echo \"done!\">%s/%s.done" %(eupathdbDir, db)+ "\n")
            f.close()
            if args.qsub:
                if args.maui:
                    cmdQsub="qsub -d `pwd`  -l walltime=10:00:00 -l nodes=1:m16G:ppn=12 %s" %(runEupathdbFile)
                else:
                    cmdQsub="qsub -cwd -V -N %s -l h_data=16G,time=24:00:00 %s" %(db,runEupathdbFile)
                os.system(cmdQsub)
                write2Log("Job for STEP6 (%s) has been submitted via qsub" %(db),gLogfile, args.quiet)
        else:
            os.chdir(eupathdbDir)
            os.system(cmd)
            eupathdbReads=set()
            eupathdbReads=nMicrobialReads(eupathdbFile,readLength,eupathdbFileFiltered)
            nEupathdbReads=len(eupathdbReads)
            write2Log("--identified %s reads mapped %s genomes" %(nEupathdbReads,db) ,gLogfile,args.quiet)
           
            excludeReadsFromFasta(inFasta,eupathdbReads,afterFasta)
            inFasta=afterFasta
            nReadsEP+=nEupathdbReads

    if not args.qsubArray and not args.qsub:
        os.rename(eupathdbDir+"%s_after_tritryp.fasta" %(basename), unaccountedReadsFasta)

        if not args.dev:
            os.remove(afterVirusFasta)
        for db in ["ameoba",
               "crypto",
               "giardia",
               "microsporidia",
               "piroplasma",
               "plasmo",
               "toxo",
               "trich"]:
            os.remove(eupathdbDir+"%s_after_%s.fasta" %(basename,db))



    if not args.qsub and  not args.qsubArray:
        write2File("done!",args.dir+"/step6_microbiomeProfile.done")

    if not args.qsub and  not args.qsubArray and not args.nonReductive:
        write2Log("In toto : %s reads mapped to microbial genomes" %(nReadsBacteria+nReadsVirus+nReadsEP) ,gLogfile,args.quiet)
        
        nTotalReads=nLowQReads+nLowCReads+n_rRNAReads+nlostHumanReads+nRepeatReads+nNCLReads+nReadsImmuneTotal+nReadsBacteria+nReadsVirus+nReadsEP
        write2Log("Summary: The ROP protocol is able to account for %s reads" %(nTotalReads) ,gLogfile,args.quiet)
        write2Log("***Unaccounted reads (not explained by ROP) are saved to %s" %(unaccountedReadsFasta) ,gLogfile,args.quiet)
        
        write2Log("***Log file with all the commands used is available here: %s" %(cmdLogfile) ,gLogfile,args.quiet)

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


#tools
write2Log("The list of the tools used by ROP and the paramemers are provided below" ,toolsLogfile,"False")
write2Log("************" ,toolsLogfile,"False")
write2Log("**We have used FastQC (version 0.0.13, with the default parameters) downloaded from  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ to filter out low quality reads" ,toolsLogfile,"False")
write2Log("**We have used SEQLEAN (seqclean-x86_64, with the default parameters) downloaded from https://sourceforge.net/projects/seqclean/ to filter out low complexity reads , " ,toolsLogfile,"False")
write2Log("**We have used Megablast (BLAST+ version 2.2.30, with the following parameters: task	=	megablast,	use_index	=	 true; -outfmt 6 ;-evalue 1e-05; perc_identity	=	100) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to filter out reads mapped to rRNA	repeat sequence	(HSU13369 Human ribosomal DNA	complete repeating	unit) " ,toolsLogfile,"False")
write2Log("**We have used Bowtie2(version 2.0.5, with the following parameters: -k 1; -p 8; -f ) downloaded from  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml to identify lost human reads mapped to reference transcriptome and genome (Ensembl GRCh37GRCh37/hg19)" ,toolsLogfile,"False")
write2Log("**We have used Megablast (BLAST+ version 2.2.30, with the following options : task=megablast, use_index=true, -outfmt 6 -evalue 1e-05, perc_identity	=	90) downloaded from  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/  to identify lost repeat reads mapped to database of 	repeat	sequences (RepBase20.07) " ,toolsLogfile,"False")
write2Log("**We have used CIRI (version 1.2 with the following parameters : -S -I ) downloaded from https://sourceforge.net/projects/ciri/ to identify reads from circRNAs" ,toolsLogfile,"False")
write2Log("**We have used IgBLAST (version v 1.4.0 with the following parameters: -germline_db_V;	germline_db_D; -germline_db_J; -organism=human;	-outfmt 7 std qseq sseq; -evalue = 1e-05 )  downloaded from http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/igblast/release/1.4.0/  to identify immune reads spanningBCR/TCR	 receptor gene	 rearrangement	in	 the	variable	domain (V(D)J	recombinations)",toolsLogfile,"False")
write2Log("**We have used Megablast (BLAST+ version 2.2.30, with the following parameters: task=megablast, use_index=true; -outfmt 6 ;-evalue 1e-05) downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ to identify microbial reads mapped onto the microbial genomes (Bacteria, Viruses, and Eukaryotic Pathogens)" ,toolsLogfile,"False")
write2Log("**We have used Metaphlan2 (version 2.2.0, with the following parameters: --mpa_pkl ; --bowtie2_exe;  --input_type multifasta; --bowtie2db ; -t reads_map/rel_ab ;  --nproc 8) downloaded from https://bitbucket.org/biobakery/metaphlan2 to obtain taxonomic profile of microbial communities" ,toolsLogfile,"False")
write2Log("************" ,toolsLogfile,"False")
write2Log("For more information about the paramemers and databases used by ROP, please see the preprint : Dumpster diving in RNA-sequencing to find the source of every last read http://biorxiv.org/content/early/2016/05/13/053041" ,toolsLogfile,"False")


write2Log("********************",gLogfile,args.quiet)
write2Log("Important: ROP relies on  several open source tools that were developed by other groups. These components are (c) their respective developers and are redistributed with ROP to provide ease-of-use. The list of the tools used by ROP and the parameters/reference databases are provided here: %s " %(toolsLogfile) ,gLogfile,args.quiet)






if args.rezip:
    write_gzip(branch_point_file, afterlostHumanFastaGzip)
    os.remove(branch_point_file)
