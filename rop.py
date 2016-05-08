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
####################################################################
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
def write2Log(message,logFile,option):
    if not option:
        print message
    logFile.write(message)
    logFile.write("\n")

#######################################################################
def nReadsImmune(inFile):
    readsImmune=set()
    with open(inFile,'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
            print line
            read=line[1]
            eValue=float(line[10])
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
#1:25169311|25169341     1       25169311        25169341        2       1_2_0   6       0.400   n/a     /n/a    SRR1146076.13939638,SRR1146076.25457964,
def nCirrcularReads(inFile):
    reads=set()
    with open(inFile,'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
            for kv in line[9].split(","):
                reads.add(kv)
            
    return reads


print "*********************************************"
print "ROP is a computational protocol aimed to discover the source of all reads, originated from complex RNA molecules, recombinant antibodies and microbial communities. Written by Serghei Mangul (smangul@ucla.edu) and Harry Taegyun Yang (harry2416@gmail.com), University of California, Los Angeles (UCLA). (c) 2016. Released under the terms of the General Public License version 3.0 (GPLv3)"
print ""
print "For more details see:"
print "http://serghei.bioinformatics.ucla.edu/rop/"
print "https://github.com/smangul1/rop/wiki"
print "*********************************************"




ap = argparse.ArgumentParser('python rop.py')


ap.add_argument('unmappedReads', help='unmapped Reads in the fastq format')
ap.add_argument('dir', help='directory (absolute path) to save results of the analysis')

ap.add_argument("--qsub", help="submit qsub jobs on hoffman2 cluster",
                action="store_true")
ap.add_argument("--qsubArray", help="prepare qsub scripts to be run later using job array",
                action="store_true")
ap.add_argument("--b", help="if unmapped reads are in bam format",
                action="store_true")
ap.add_argument("--skipLowq", help="skip filtering step",
                action="store_true")
ap.add_argument("--skipQC", help="skip entire QC step : filtering low-quality, low-complexity and rRNA reads (reads matching rRNA repeat unit)",
                action="store_true")
#ap.add_argument("--circRNA", help="enable CIRI for circular RNA detection ", action="store_true")
ap.add_argument("--immune", help = "Only TCR/BCR immune gene analysis will be performed", action = "store_true")
#ap.add_argument("--gzip", help = "Gzip the fasta files after filtering step", action = "store_true")
ap.add_argument("--quiet", help = "suppress progress report and warnings", action = "store_true")
ap.add_argument("--dev", help = "keep intermediate files", action = "store_true")
ap.add_argument("--license", help= "Show ROP License Information", action = "store_true")
ap.add_argument("--version", help= "Show ROP version", action = "store_true")




args = ap.parse_args()


##################################
# VERSION
##################################
if args.version:
    print "ROP version 1.0"


##################################
# LICENSE
##################################

if args.license:
    print """
        Read Origin Protocol is a computational protocol for profiling the composition of unmapped reads, which failed to map to the human references. ROP profiles repeats, circRNAs, gene fusions, trans-splicing events, recombined B and T cell receptors and microbial communities.
        Copyright (C) 2016  Serghei Mangul and Harry Taegyun Yang
        
        Read Origin Protocol is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.
        
        Read Origin Protocol is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.
        
        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <http://www.gnu.org/licenses/>.
        """
    sys.exit(0)



##################################
#main code
##################################


#relative path to absolute path
args.unmappedReads=os.path.abspath(args.unmappedReads)
args.dir=os.path.abspath(args.dir)


#basename
basename=os.path.splitext(os.path.basename(args.unmappedReads))[0]





#analysis directories
QCDir=args.dir+"/QC/"



humanDir=args.dir+"/human/"

lostHumanDir=humanDir+"/lostHuman/"
lostRepeatDir=humanDir+"/lostRepeat/"
bcrDir=args.dir+"/BCR/"
tcrDir=args.dir+"/TCR/"

NCL_Dir=args.dir+"/NCL/"

ighDir=args.dir+"/BCR/IGH/"
igkDir=args.dir+"/BCR/IGK/"
iglDir=args.dir+"/BCR/IGL/"

tcraDir=args.dir+"/TCR/TCRA/"
tcrbDir=args.dir+"/TCR/TCRB/"
tcrdDir=args.dir+"/TCR/TCRD/"
tcrgDir=args.dir+"/TCR/TCRG/"

microbiomeDir=args.dir+"/microbiome/"

bacteriaDir=args.dir+"/microbiome/bacteria/"
virusDir=args.dir+"/microbiome/virus/"
eupathdbDir=args.dir+"/microbiome/eupathdb/"



if not os.path.exists(QCDir):
    os.makedirs(QCDir)
if not os.path.exists(lostHumanDir):
    os.makedirs(lostHumanDir)
if not os.path.exists(lostRepeatDir):
    os.makedirs(lostRepeatDir)
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

#intermediate files
unmappedFastq=args.dir+"/unmapped_"+basename+".fastq"
lowQFile=QCDir+basename+"_lowQ.fastq"
lowQFileFasta=QCDir+basename+"_lowQ.fa"
lowQCFile=QCDir+basename+"_lowQC.fa"
rRNAFile=QCDir+basename+"_rRNA_blastFormat6.csv"
afterrRNAFasta=QCDir+basename+"_after_rRNA.fasta"
afterlostHumanFasta=lostHumanDir+basename+"_after_rRNA_lostHuman.fasta"
afterImmuneFasta=bcrDir+basename+"_afterImmune.fasta"
afterBacteraFasta=bacteriaDir+basename+"_afterBacteria.fasta"
afterVirusFasta=virusDir+basename+"_afterVirus.fasta"


gBamFile=lostHumanDir+basename+"_genome.sam"
tBamFile=lostHumanDir+basename+"_transcriptome.sam"
repeatFile=lostRepeatDir+basename+"_lostRepeats_blastFormat6.csv"
afterlostRepeatFasta=lostRepeatDir+basename+"_after_lostRepeat.fasta"

afterNCLFasta=NCL_Dir+basename+"_after_NCL.fasta"


NCL_CIRI_file=NCL_Dir + basename + "_NCL_CIRI_after_bwa.sam"
after_NCL_CIRI_file_prefix = NCL_Dir + "/"+basename + "_circRNA.csv"
if os.path.exists(after_NCL_CIRI_file_prefix):
    os.remove(after_NCL_CIRI_file_prefix)

ighFile=ighDir+basename+"_IGH_igblast.csv"
igkFile=igkDir+basename+"_IGK_igblast.csv"
iglFile=iglDir+basename+"_IGL_igblast.csv"
tcraFile=tcraDir+basename+"_TCRA_igblast.csv"
tcrbFile=tcrbDir+basename+"_TCRB_igblast.csv"
tcrdFile=tcrdDir+basename+"_TCRD_igblast.csv"
tcrgFile=tcrgDir+basename+"_TCRG_igblast.csv"

#log files
logQC=QCDir+basename+"_QC.log"
logrRNA=QCDir + basename + "_rRNA.log"
logHuman=lostHumanDir + basename + "_lostHuman.log"
logNCL=lostHumanDir + basename + "_lostHuman.log"


bacteriaFile=bacteriaDir+basename+"_bacteria_blastFormat6.csv"
virusFile=virusDir+basename+"_virus_blastFormat6.csv"

bacteriaFileFiltered=bacteriaDir+basename+"_bacteriaFiltered_blastFormat6.csv"
virusFileFiltered=virusDir+basename+"_virusFiltered_blastFormat6.csv"

gLog=args.dir+"/"+basename+".log"
gLogfile=open(gLog,'w')

tLog=args.dir+"/"+"numberReads_"+basename+".log"
tLogfile=open(tLog,'w')



cmdLog=args.dir+"/"+"dev.log"
cmdLogfile=open(cmdLog,'w')


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


os.chdir(args.dir)

#######################################################################################################################################
if args.skipQC:
    afterrRNAFasta=args.unmappedReads
else:
    if args.b:
        if args.skipLowq:
            bam2fasta(codeDir,args.unmappedReads,lowQFileFasta)
        else:
            bam2fastq(codeDir,args.unmappedReads,unmappedFastq)
    else :
        unmappedFastq=args.unmappedReads


    #######################################################################################################################################
    if args.skipLowq==False:
        
        #number of reads
        with open(unmappedFastq) as f:
            for i, l in enumerate(f):
                pass
        n=(i + 1)/4
        
        message="Processing %s unmapped reads" %(n)
        write2Log(message,gLogfile,args.quiet)
        
        
        
        
        
        #lowQ
        write2Log("1. Quality Control...",gLogfile,args.quiet)
        cmd=codeDir+"/tools/fastq_quality_filter -v -Q 33 -q 20 -p 75 -i %s -o %s > %s \n" %(unmappedFastq,lowQFile,logQC)
        write2Log(cmd,cmdLogfile,"False")
        os.system(cmd)
        if args.b:
            os.remove(unmappedFastq)
        
        
        readLength=0
        #Convert from fastq to fasta
        fastafile=open(lowQFileFasta,'w')
        fastqfile = open(lowQFile, "rU")
        nAfterLowQReads=0
        for record in SeqIO.parse(fastqfile,"fastq"):
            readLength=len(record) #assumes the same length, will not work for Ion Torrent or Pac Bio
            fastafile.write(str(">"+record.name)+"\n")
            fastafile.write(str(record.seq)+"\n")
            nAfterLowQReads+=1
        fastafile.close()
        nLowQReads=n-nAfterLowQReads
        write2Log("--filtered %s low quality reads" %(nLowQReads) ,gLogfile,args.quiet)




    os.chdir(QCDir)
    #lowC
    cmd="export PATH=$PATH:%s/tools/seqclean-x86_64/bin" %(codeDir)
    os.system(cmd)
    
    cmd=codeDir+"/tools/seqclean-x86_64/seqclean %s -l 50 -M -o %s 2>>%s" %(lowQFileFasta, lowQCFile,logQC)
    write2Log(cmd,cmdLogfile,"False")
    os.system(cmd)
    
    cmd = "rm -rf %s/cleaning_1/ ; rm -f %s/*.cln ; rm -f %s/*.cidx; rm -f %s/*.sort" % (QCDir,QCDir,QCDir,QCDir)
    os.system(cmd)
    proc = subprocess.Popen(["grep trashed %s | awk -F \":\" '{print $2}'" %(logQC) ], stdout=subprocess.PIPE, shell=True)
    (nLowCReadsTemp, err) = proc.communicate()
    nLowCReads=int(nLowCReadsTemp.rstrip().strip())
    write2Log("--filtered %s low complexity reads (e.g. ACACACAC...)" %(nLowCReads) ,gLogfile,args.quiet)
    
    
    
    
    
    
    
    #rRNA
    cmd="%s/tools/blastn -task megablast -index_name %s/db/rRNA/rRNA -use_index true -query %s -db %s/db/rRNA/rRNA  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,lowQCFile,codeDir,rRNAFile)
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
    
    
    message="Number of entries in %s is %s" %(rRNAFile,n_rRNATotal)
    write2Log(message,cmdLogfile,"False")
    
    if not args.dev:
        os.remove(lowQFile)
        os.remove(lowQCFile)
        os.remove(lowQFileFasta)
        os.remove(rRNAFile)






#######################################################################################################################################
#2. Remaping to human references...
write2Log("2. Remaping to human references...",cmdLogfile,"False")
write2Log("2. Remaping to human references...",gLogfile,args.quiet)

cmdGenome="%s/tools/bowtie2 -k 1 -p 8 -f -x %s/db/human/Bowtie2Index/genome -U %s 2>%s | %s/tools/samtools view -SF4 -   >%s" %(codeDir,codeDir, afterrRNAFasta,logHuman,codeDir,gBamFile)

#transcriptome
cmdTranscriptome="%s/tools/bowtie2  -k 1 -f -p 8 -x %s/db/human/Bowtie2Index/hg19KnownGene.exon_polya200 -U %s 2>%s | %s/tools/samtools view -SF4 -  >  %s " %(codeDir,codeDir, afterrRNAFasta,logHuman, codeDir,tBamFile)
write2Log(cmdGenome,cmdLogfile,"False")
write2Log(cmdTranscriptome,cmdLogfile,"False")






os.system(cmdGenome)
os.system(cmdTranscriptome)




nlostHumanReads_10=0

lostHumanReads = set()


with open(gBamFile,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        if int(line[16].split(':')[2])<3:
            lostHumanReads.add(line[0])

with open(tBamFile,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        if int(line[16].split(':')[2])<3:
            lostHumanReads.add(line[0])



excludeReadsFromFasta(afterrRNAFasta,lostHumanReads,afterlostHumanFasta)
nlostHumanReads=len(lostHumanReads)
write2Log("--identified %s lost human reads from unmapped reads " %(len(lostHumanReads)) ,gLogfile,args.quiet)


if not args.dev:
    os.remove(afterrRNAFasta)
    os.remove(gBamFile)
    os.remove(tBamFile)





#######################################################################################################################################
#3. Maping to repeat sequences...
write2Log("3. Maping to repeat sequences...",cmdLogfile,"False")
write2Log("3. Maping to repeat sequences...",gLogfile,args.quiet)

#TO DO : make all fasta ->gzip
#gzip -dc %s | , query -

cmd="%s/tools/blastn -task megablast -index_name %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa -use_index true -query %s -db %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,afterlostHumanFasta,codeDir,repeatFile)


if args.qsub or args.qsubArray:
    f = open(runLostRepeatFile,'w')
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_lostRepeat.done \n" %(lostRepeatDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N lostRepeat -l h_data=16G,time=10:00:00 %s" %(runLostRepeatFile)
        os.system(cmdQsub)
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
    write2Log("-Identify %s lost repeat sequences from unmapped reads" %(nRepeatReads) ,gLogfile,args.quiet)
    write2Log("***Note : Repeat sequences classification into classes (e.g. LINE) and families (e.g. Alu) will be available in next release" ,gLogfile,args.quiet)
    
    
    excludeReadsFromFasta(afterlostHumanFasta,lostRepeatReads,afterlostRepeatFasta)
    
    if not args.dev:
        os.remove(afterlostHumanFasta)

#######################################################################################################################################
#3. Non-co-linear RNA profiling
write2Log("3. Non-co-linear RNA profiling",cmdLogfile,"False")
write2Log("3. Non-co-linear RNA profiling",gLogfile,args.quiet)
write2Log("***Note : Trans-spicing and gene fusions  are currently not supported, but will be in the next release.",gLogfile,args.quiet)

os.chdir(NCL_Dir)
NCL_reads=set()

cmd="%s/tools/bwa mem -T -S %s/db/human/BWAIndex/genome.fa %s > %s 2>%s \n" %(codeDir,codeDir,afterlostRepeatFasta,NCL_CIRI_file,logNCL)
cmd = cmd + "perl %s/tools/CIRI_v1.2.pl -S -I %s -O %s -F %s/db/human/BWAIndex/genome.fa 1>>%s 2>>%s" %(codeDir,NCL_CIRI_file,after_NCL_CIRI_file_prefix,codeDir,logNCL,logNCL)
write2Log(cmd,cmdLogfile,"False")

if args.qsub or args.qsubArray:
    f = open(runNCL_CIRIfile,'w')
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_NCL_CIRI.done \n" %(NCL_Dir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N NCL_CIRI -l h_data=8G,time=10:00:00 %s" %(runNCL_CIRIfile)
        os.system(cmdQsub)
else:
    os.system(cmd)
    NCL_reads=nCirrcularReads(after_NCL_CIRI_file_prefix)
    nReadsNCL=len(NCL_reads)
    write2Log("--identified %s reads from circRNA" %(nReadsNCL) ,gLogfile,args.quiet)

excludeReadsFromFasta(afterlostRepeatFasta,NCL_reads,afterNCLFasta)
    
if not args.dev:
        os.remove(afterlostRepeatFasta)



#######################################################################################################################################
#4. T and B lymphocytes profiling
immuneReads=set()

write2Log("4a. B lymphocytes profiling...",cmdLogfile,"False")
write2Log("4a. B lymphocytes profiling...",gLogfile,args.quiet)

#IGH-------
os.chdir(ighDir)
cmd="ln -s %s//db/antibody/internal_data/ ./" %(codeDir)
os.system(cmd)


cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/IGHV.fa -germline_db_D %s/db/antibody//IGHD.fa  -germline_db_J %s/db/antibody//IGHJ.fa -query %s -outfmt 7 -evalue 1e-05  2>temp.txt | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterNCLFasta,ighFile)
write2Log(cmd,cmdLogfile,"False")


if args.qsub or args.qsubArray:
    f = open(runIGHFile,'w')
    f.write("ln -s %s//db/antibody//internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\"> %s/%s_igh.done \n" %(ighDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N igh -l h_data=16G,time=24:00:00 %s" %(runIGHFile)
        os.system(cmdQsub)
else:
    os.chdir(ighDir)
    os.system(cmd)
    print ighFile
    
    immuneReadsIGH=nReadsImmune(ighFile)
    nReadsImmuneIGH=len(immuneReadsIGH)
    write2Log("--identified %s reads mapped to immunoglobulin heavy (IGH) locus" %(nReadsImmuneIGH) ,gLogfile,args.quiet)


#IGK---------
os.chdir(igkDir)
cmd="ln -s %s//db/antibody//internal_data/ ./" %(codeDir)

cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/IGKV.fa -germline_db_D %s/db/antibody//IGHD.fa  -germline_db_J %s/db/antibody//IGKJ.fa -query %s -outfmt 7 -evalue 1e-05 2>temp.txt | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterNCLFasta,igkFile)
write2Log(cmd,cmdLogfile,"False")

if args.qsub or args.qsubArray:
    f = open(runIGKFile,'w')
    f.write("ln -s %s//db/antibody//internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\"> %s/%s_igk.done \n" %(igkDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N igk -l h_data=16G,time=24:00:00 %s" %(runIGKFile)
        os.system(cmdQsub)
else:
    os.chdir(igkDir)
    os.system(cmd)
    immuneReadsIGK=nReadsImmune(igkFile)
    nReadsImmuneIGK=len(immuneReadsIGK)
    write2Log("--identified %s reads mapped to immunoglobulin kappa (IGK) locus " %(nReadsImmuneIGK) ,gLogfile,args.quiet)


#IGL------------
os.chdir(iglDir)
cmd="ln -s %s//db/antibody//internal_data/ ./" %(codeDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/IGLV.fa -germline_db_D %s/db/antibody//IGHD.fa  -germline_db_J %s/db/antibody//IGLJ.fa -query %s -outfmt 7 -evalue 1e-05 2>temp.txt  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterNCLFasta,iglFile)
write2Log(cmd,cmdLogfile,"False")

if args.qsub or args.qsubArray:
    f = open(runIGLFile,'w')
    f.write("ln -s %s//db/antibody//internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_igl.done \n" %(iglDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N igl -l h_data=16G,time=24:00:00 %s" %(runIGLFile)
        os.system(cmdQsub)
else:
    os.chdir(iglDir)
    os.system(cmd)
    immuneReadsIGL=nReadsImmune(iglFile)
    nReadsImmuneIGL=len(immuneReadsIGL)
    write2Log("--identified %s reads mapped to immunoglobulin lambda (IGL) locus" %(nReadsImmuneIGL) ,gLogfile,args.quiet)


##################
##################
write2Log("4b. T lymphocytes profiling...",cmdLogfile,"False")
write2Log("4b. T lymphocytes profiling...",gLogfile,args.quiet)

#TCRA-----------------
os.chdir(tcraDir)
cmd="ln -s %s//db/antibody//internal_data/ ./" %(codeDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/TRAV.fa -germline_db_D %s/db/antibody//TRBD.fa  -germline_db_J %s/db/antibody//TRAJ.fa -query %s -outfmt 7 -evalue 1e-05 2>temp.txt | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterNCLFasta,tcraFile)
write2Log(cmd,cmdLogfile,"False")
if args.qsub or args.qsubArray:
    f = open(runTCRAFile,'w')
    f.write("ln -s %s//db/antibody//internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_tcra.done \n"%(tcraDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N tcra -l h_data=16G,time=24:00:00 %s" %(runTCRAFile)
        os.system(cmdQsub)
else:
    os.chdir(tcraDir)
    os.system(cmd)
    immuneReadsTCRA=nReadsImmune(tcraFile)
    nReadsImmuneTCRA=len(immuneReadsTCRA)
    write2Log("--identified %s reads mapped to T cell receptor alpha (TCRA) locus" %(nReadsImmuneTCRA) ,gLogfile,args.quiet)



#TCRB--------------
os.chdir(tcrbDir)
cmd="ln -s %s//db/antibody//internal_data/ ./" %(codeDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/TRBV.fa -germline_db_D %s/db/antibody//TRBD.fa  -germline_db_J %s/db/antibody//TRBJ.fa -query %s -outfmt 7 -evalue 1e-05 2>temp.txt  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterNCLFasta,tcrbFile)
write2Log(cmd,cmdLogfile,"False")
if args.qsub or args.qsubArray:
    f = open(runTCRBFile,'w')
    f.write("ln -s %s//db/antibody//internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_tcrb.done \n"%(tcrbDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N tcrb -l h_data=16G,time=24:00:00 %s" %(runTCRBFile)
        os.system(cmdQsub)
else:
    os.chdir(tcrbDir)
    os.system(cmd)
    immuneReadsTCRB=nReadsImmune(tcrbFile)
    nReadsImmuneTCRB=len(immuneReadsTCRB)
    write2Log("--identified %s reads mapped to T cell receptor beta (TCRB) locus" %(nReadsImmuneTCRB) ,gLogfile,args.quiet)




#TCRD----------------
os.chdir(tcrdDir)
cmd="ln -s %s//db/antibody//internal_data/ ./" %(codeDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/TRDV.fa -germline_db_D %s/db/antibody//TRBD.fa  -germline_db_J %s/db/antibody//TRDJ.fa -query %s -outfmt 7 -evalue 1e-05 2>temp.txt  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterNCLFasta,tcrdFile)
write2Log(cmd,cmdLogfile,"False")
if args.qsub or args.qsubArray:
    f = open(runTCRDFile,'w')
    f.write("ln -s %s//db/antibody//internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_tcrd.done \n" %(tcrdDir, basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N tcrd -l h_data=16G,time=24:00:00 %s" %(runTCRDFile)
        os.system(cmdQsub)
else:
    os.chdir(tcrdDir)
    os.system(cmd)
    immuneReadsTCRD=nReadsImmune(tcrdFile)
    nReadsImmuneTCRD=len(immuneReadsTCRD)
    write2Log("--identified %s reads mapped to T cell receptor delta (TCRD) locus" %(nReadsImmuneTCRD) ,gLogfile,args.quiet)



#TCRG---------------------
os.chdir(tcrgDir)
cmd="ln -s %s//db/antibody//internal_data/ ./" %(codeDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/antibody/TRGV.fa -germline_db_D %s/db/antibody/TRBD.fa  -germline_db_J %s/db/antibody/TRGJ.fa -query %s -outfmt 7 -evalue 1e-05 2>temp.txt  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterNCLFasta,tcrgFile)
write2Log(cmd,cmdLogfile,"False")

if args.qsub or args.qsubArray:
    f = open(runTCRGFile,'w')
    f.write("ln -s %s//db/antibody/internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_tcrg.done \n" %(tcrgDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N tcrg -l h_data=16G,time=24:00:00 %s" %(runTCRGFile)
        os.system(cmdQsub)
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
    
    immuneReads=set().union(immuneReadsTCRA,immuneReadsTCRB,immuneReadsTCRD,immuneReadsTCRG)
    excludeReadsFromFasta(afterNCLFasta,immuneReads,afterImmuneFasta)
    if not args.dev:
        os.remove(afterNCLFasta)




#######################################################################################################################################
#5. Microbiome profiling...
if not args.immune:
    write2Log("5.  Microbiome profiling...",cmdLogfile,"False")
    write2Log("5.  Microbiome profiling...",gLogfile,args.quiet)
    
    #bacteria ----------
    cmd="%s/tools/blastn -task megablast -index_name %s/db/microbiome/bacteria/bacteria -use_index true -query %s -db %s/db/microbiome/bacteria/bacteria  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s 2>temp.txt" %(codeDir,codeDir,afterImmuneFasta,codeDir,bacteriaFile)
    write2Log(cmd,cmdLogfile,"False")
    
    if args.qsub or args.qsubArray:
        f = open(runBacteriaFile,'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_bacteria.done \n"%(bacteriaDir, basename))
        f.close()
        if args.qsub:
            cmdQsub="qsub -cwd -V -N bacteria -l h_data=16G,time=24:00:00 %s" %(runBacteriaFile)
            os.system(cmdQsub)
    else:
        os.chdir(bacteriaDir)
        os.system(cmd)
        bacteriaReads=nMicrobialReads(bacteriaFile,readLength,bacteriaFileFiltered)
        nReadsBacteria=len(bacteriaReads)
        write2Log("--identified %s reads mapped bacterial genomes" %(nReadsBacteria) ,gLogfile,args.quiet)
        excludeReadsFromFasta(afterImmuneFasta,bacteriaReads,afterBacteraFasta)
    
    
    #MetaPhlAn
    #    cmd=" python %s/tools/metaphlan.py
    
    
    #metaphlan2.py metagenome.fastq --input_type fastq
    
    
    #virus-----------
    cmd="%s/tools/blastn -task megablast -index_name %s/db/microbiome/virus/viruses -use_index true -query %s -db %s/db/microbiome/virus/viruses  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s 2>temp" %(codeDir,codeDir,afterBacteraFasta,codeDir,virusFile)
    write2Log(cmd,cmdLogfile,"False")
    if args.qsub or args.qsubArray:
        f = open(runVirusFile,'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_virus.done \n" %(virusDir, basename))
        f.close()
        if args.qsub:
            cmdQsub="qsub -cwd -V -N virus -l h_data=16G,time=24:00:00 %s" %(runVirusFile)
            os.system(cmdQsub)
    else:
        os.chdir(virusDir)
        os.system(cmd)
        virusReads=nMicrobialReads(virusFile,readLength,virusFileFiltered)
        nReadsVirus=len(virusReads)
        write2Log("--identified %s reads mapped viral genomes" %(nReadsVirus) ,gLogfile,args.quiet)
        excludeReadsFromFasta(afterBacteraFasta,virusReads,afterVirusFasta)


    #eukaryotic pathogens----------------
    dbList=["ameoba",
            "crypto",
            "giardia",
            "microsporidia",
            "piroplasma",
            "plasmo",
            "toxo",
            "trich",
            "tritryp"]




    inFasta=afterVirusFasta
    nReadsEP=0
    
    for db in dbList:
        eupathdbFile=eupathdbDir+basename+"_"+db+"_blastFormat6.csv"
        eupathdbFileFiltered=eupathdbDir+basename+"_"+db+"Filtered_blastFormat6.csv"
        runEupathdbFile=eupathdbDir+"/run_"+basename+"_"+db+".sh"
        
        
        cmd="%s/tools/blastn -task megablast -index_name %s/db/microbiome/eupathdb/%s -use_index true -query %s -db %s/db/microbiome/eupathdb/%s  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s 2>temp.txt" %(codeDir,codeDir,db,inFasta,codeDir,db,eupathdbFile)
        write2Log(cmd,cmdLogfile,"False")
        if args.qsub or args.qsubArray:
            f = open(runEupathdbFile,'w')
            f.write(cmd+"\n")
            f.write("echo \"done!\">%s/%s.done" %(eupathdbDir, db)+ "\n")
            f.close()
            if args.qsub:
                cmdQsub="qsub -cwd -V -N %s -l h_data=16G,time=24:00:00 %s" %(db,runEupathdbFile)
                os.system(cmdQsub)
        else:
            os.chdir(eupathdbDir)
            os.system(cmd)
            eupathdbReads=set()
            eupathdbReads=nMicrobialReads(eupathdbFile,readLength,eupathdbFileFiltered)
            nEupathdbReads=len(eupathdbReads)
            write2Log("--identified %s reads mapped %s genomes" %(nEupathdbReads,db) ,gLogfile,args.quiet)
            afterFasta=eupathdbDir+"%s_after_%s.fasta" %(basename,db)
            excludeReadsFromFasta(inFasta,eupathdbReads,afterFasta)
            inFasta=afterFasta
            nReadsEP+=nEupathdbReads


if not args.qsub and  not args.qsubArray:
    write2Log("In toto : %s reads mapped to microbial genomes" %(nReadsBacteria+nReadsVirus+nReadsEP) ,gLogfile,args.quiet)
    nTotalReads=nLowQReads+nLowCReads+n_rRNAReads+nlostHumanReads+nRepeatReads+nReadsImmuneTotal+nReadsBacteria+nReadsVirus+nReadsEP
    write2Log("Summary:   The ROP protocol is able to account for %s reads" %(nTotalReads) ,gLogfile,args.quiet)
    
    message=basename+","+str(n)+","+str(nLowQReads)+","+str(nLowCReads)+","+str(n_rRNAReads)+","+str(nlostHumanReads)+","+str(nRepeatReads)+","+str(nReadsImmuneTotal)+","+str(nReadsBacteria+nReadsVirus+nReadsEP)
    tLogfile.write(message)
    tLogfile.write("\n")




tLogfile.close()
gLogfile.close()
cmdLogfile.close()


