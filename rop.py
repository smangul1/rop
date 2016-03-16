import sys
import csv
import os
import argparse
from Bio import SeqIO # module needed for sequence input
import pysam




####################################################################
def excludeReadsFromFasta(inFasta,reads,outFasta):

    fasta_sequences = SeqIO.parse(open(inFasta),'fasta')
    with open(outFasta, "w") as f:
        for seq in fasta_sequences:
            name = seq.name
            if name not in reads:
                SeqIO.write([seq], f, "fasta")





ap = argparse.ArgumentParser()
ap.add_argument('unmappedReads', help='unmapped Reads in the fastq format')
ap.add_argument('dir', help='directory to save results of the analysis')

ap.add_argument("--qsub", help="submit qsub jobs on hoffman2 cluster",
                    action="store_true")

args = ap.parse_args()



#if args.qsub:
#    print "weklwjerklwjerlkjw





#number of reads
with open(args.unmappedReads) as f:
    for i, l in enumerate(f):
        pass
n=(i + 1)/4


basename=os.path.splitext(os.path.basename(args.unmappedReads))[0]


print "Number of unmapped reads",n


#codeDir
codeDir=os.path.dirname(os.path.realpath(__file__))

#run file to save all commands
runFile=args.dir+"/commands_"+basename+".sh"




#analysis directories
QCDir=args.dir+"/QC/"
lostHumanDir=args.dir+"/lostHuman/"
lostRepeatDir=args.dir+"/lostRepeat/"
bcrDir=args.dir+"/BCR/"
tcrDir=args.dir+"/TCR/"


ighDir=args.dir+"/BCR/IGH/"
igkDir=args.dir+"/BCR/IGK/"
iglDir=args.dir+"/BCR/IGL/"

tcraDir=args.dir+"/TCR/TRA/"
tcrbDir=args.dir+"/TCR/TRB/"
tcrdDir=args.dir+"/TCR/TRD/"
tcrgDir=args.dir+"/TCR/TRG/"

microbiomeDir=args.dir+"/microbiome/"

bacteriaDir=args.dir+"/microbiome/bacteria/"
virusDir=args.dir+"/microbiome/virus/"



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
for i in [ighDir,igkDir,iglDir,tcraDir,tcrbDir,tcrdDir,tcrgDir]:
    if not os.path.exists(i):
        os.makedirs(i)
if not os.path.exists(microbiomeDir):
    os.makedirs(microbiomeDir)
if not os.path.exists(bacteriaDir):
    os.makedirs(bacteriaDir)
if not os.path.exists(virusDir):
    os.makedirs(virusDir)

#intermediate files
lowQFile=QCDir+basename+"_lowQ.fastq"
lowQFileFasta=QCDir+basename+"_lowQ.fa"
lowQCFile=QCDir+basename+"_lowQC.fa"
rRNAFile=QCDir+basename+"_rRNA_blastFormat6.csv"
afterrRNAFasta=QCDir+basename+"_after_rRNA.fasta"
afterlostHumanFasta=lostHumanDir+basename+"_after_lostHuman.fasta"
gBamFile=lostHumanDir+basename+"_genome.bam"
tBamFile=lostHumanDir+basename+"_transcriptome.bam"
repeatFile=lostRepeatDir+basename+"_lostRepeats_blastFormat6.csv"
afterlostRepeatFasta=lostRepeatDir+basename+"_after_lostRepeat.fasta"
ighFile=ighDir+basename+"_IGH_igblast.csv"
igkFile=igkDir+basename+"_IGK_igblast.csv"
iglFile=iglDir+basename+"_IGL_igblast.csv"
tcraFile=tcraDir+basename+"_TCRA_igblast.csv"
tcrbFile=tcrbDir+basename+"_TCRB_igblast.csv"
tcrdFile=tcrdDir+basename+"_TCRD_igblast.csv"
tcrgFile=tcrgDir+basename+"_TCRG_igblast.csv"

bacteriaFile=bacteriaDir+basename+"_bacteria_blastFormat6.csv"
virusFile=virusDir+basename+"_virus_blastFormat6.csv"


#runFiles
runFile=args.dir+"/commands_"+basename+".sh"
runLostHumanFile=lostHumanDir+"/runLostHuman_"+basename+".sh"
runLosRepeatFile=lostRepeatDir+"/runLostRepeat_"+basename+".sh"
runIGHFile=ighDir+"/runIGH_"+basename+".sh"
runIGKFile=igkDir+"/runIGK_"+basename+".sh"
runIGLFile=iglDir+"/runIGL_"+basename+".sh"
runTCRAFile=tcraDir+"/runTCRA_"+basename+".sh"
runTCRBFile=tcrbDir+"/runTCRB_"+basename+".sh"
runTCRDFile=tcrdDir+"/runTCRD_"+basename+".sh"
runTCRGFile=tcrgDir+"/runTCRG_"+basename+".sh"
runBacteriaFile=bacteriaDir +"/runBacteria_"+basename+".sh"
runVirusFile=virusDir +"/runVirus_"+basename+".sh"



f = open(runFile,'w')


#######################################################################################################################################

#lowQ
print "*****************************Running FASTX to filter low quality reads******************************"
cmd=codeDir+"/tools/fastq_quality_filter -v -Q 33 -q 20 -p 75 -i %s -o %s \n" %(args.unmappedReads,lowQFile)
print "Run ", cmd
os.system(cmd)
print "Save reads after filtering low quality reads to ", lowQFile
f.write(cmd+"\n" )


#Convert from fastq to fasta
print "**********************************Convert ",lowQFile,"to ",lowQFileFasta,"****************************"

fastafile=open(lowQFileFasta,'w')

fastqfile = open(lowQFile, "rU")
for record in SeqIO.parse(fastqfile,"fastq"):
    fastafile.write(str(">"+record.name)+"\n")
    fastafile.write(str(record.seq)+"\n")
fastafile.close()

#need to add command to commands_....sh


#lowC
print "*****************************Running FASTX to filter low complexity reads******************************"
cmd=codeDir+"/tools/seqclean-x86_64/seqclean %s -l 50 -M -o %s" %(lowQFileFasta, lowQCFile)
print "Run ", cmd
os.system(cmd)
print "Save reads after filtering low complexity (e.g. ACACACAC...) reads to ", lowQCFile
f.write(cmd+"\n" )

#megablast rRNA
print "*****************************Identify reads from RNA repeat unit******************************"
cmd="%s/tools/blastn -task megablast -index_name %s/db/rRNA/rRNA -use_index true -query %s -db %s/db/rRNA/rRNA  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,lowQCFile,codeDir,rRNAFile)
print "Run :", cmd
os.system(cmd)
f.write(cmd+"\n")


rRNAReads = set()

with open(rRNAFile,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        element=line[0]
        identity=float(line[2])
        alignmentLength=float(line[3])
        eValue=float(line[10])
        if eValue<1e-05 and alignmentLength==100 and identity>=94:
            rRNAReads.add(element)

excludeReadsFromFasta(lowQCFile,rRNAReads,afterrRNAFasta)


#######################################################################################################################################
print "*****************************Identify lost human reads******************************"





os.system(". /u/local/Modules/default/init/modules.sh \n")
os.system("module load bowtie2/2.1.0 \n")
os.system("module load samtools \n")


#genome
cmdGenome="bowtie2 -k 1 --very-sensitive -p 8 -f -x %s/db/human/Bowtie2Index/genome -U %s | samtools view -bSF4 - | samtools sort -  %s" %(codeDir, afterrRNAFasta,os.path.splitext(gBamFile)[0])
print "Run: ",cmd
#transcriptome
cmdTranscriptome="bowtie2 -k 1 --very-sensitive -f -p 8 -x %s/db/human/Bowtie2Index/hg19KnownGene.exon_polya200 -U %s | samtools view -bSF4 - | samtools sort -  %s" %(codeDir, afterrRNAFasta,os.path.splitext(tBamFile)[0])
print "Run: ",cmd

if args.qsub:
    f = open(runLostHumanFile,'w')
    f.write(". /u/local/Modules/default/init/modules.sh \n")
    f.write("module load bowtie2/2.1.0 \n")
    f.write("module load samtools \n")
    f.write(cmdGenome+"\n")
    f.write(cmdGenome+"\n")
    f.write("samtools index %s \n" %gBamFile)
    f.write("samtools index %s \n" %tBamFile)

    f.close()
else:
    os.system(cmdGenome)
    os.system(cmdTranscriptome)
    os.system("samtools index %s \n" %gBamFile)



if args.qsub:
    f.write(". /u/local/Modules/default/init/modules.sh \n")
    f.write("module load bowtie2/2.1.0 \n")
    f.write("module load samtools \n")
    f.write(cmd+"\n")
    f.write("samtools index %s \n" %tBamFile)
    f.close()
else:
    os.system(cmd)
    os.system("samtools index %s \n" %tBamFile)



if not args.qsub:

    lostHumanReads = set()


    samfile = pysam.AlignmentFile(gBamFile, "rb")
    for r in samfile.fetch():
        for tag in r.tags:
            if tag[0] == 'NM':
                if int(tag[1])<=6:
                    lostHumanReads.add(r.query_name)
    samfile = pysam.AlignmentFile(tBamFile, "rb")
    for r in samfile.fetch():
        for tag in r.tags:
            if tag[0] == 'NM':
                if int(tag[1])<=6:
                    lostHumanReads.add(r.query_name)

    excludeReadsFromFasta(afterrRNAFasta,lostHumanReads,afterlostHumanFasta)

#######################################################################################################################################
print "*****************************Identify lost repeat reads******************************"
cmd="%s/tools/blastn -task megablast -index_name %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa -use_index true -query %s -db %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,afterlostHumanFasta,codeDir,repeatFile)
print "Run :", cmd
os.system(cmd)

lostRepeatReads = set()

with open(repeatFile,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        element=line[0]
        identity=float(line[2])
        alignmentLength=float(line[3])
        eValue=float(line[10])
        if eValue<1e-05 and alignmentLength>=80 and identity>=90:
            lostRepeatReads.add(element)

excludeReadsFromFasta(afterlostHumanFasta,lostRepeatReads,afterlostRepeatFasta)

#######################################################################################################################################
print "*****************************Identify NCL events******************************"
#TO DO!!!!!!!!!

#######################################################################################################################################

print "*****************************Identify VDJ recombinations from BCR and TCR******************************"
cmd="ln -s %s//db/BCRTCR/internal_data/ ./" %(codeDir)
os.system(cmd)

#IGH
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/IGHV.fa -germline_db_D %s/db/BCRTCR/IGHD.fa  -germline_db_J %s/db/BCRTCR/IGHJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterlostRepeatFasta,ighFile)
os.system(cmd)
print "Run: ",cmd

#IGK
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/IGKV.fa -germline_db_D %s/db/BCRTCR/IGHD.fa  -germline_db_J %s/db/BCRTCR/IGKJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterlostRepeatFasta,igkFile)
os.system(cmd)
print "Run: ",cmd

#IGL
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/IGLV.fa -germline_db_D %s/db/BCRTCR/IGHD.fa  -germline_db_J %s/db/BCRTCR/IGLJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterlostRepeatFasta,iglFile)
os.system(cmd)
print "Run: ",cmd

#TCRA
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/TRAV.fa -germline_db_D %s/db/BCRTCR/TRBD.fa  -germline_db_J %s/db/BCRTCR/TRAJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterlostRepeatFasta,tcraFile)
os.system(cmd)
print "Run: ",cmd


#TCRB
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/TRBV.fa -germline_db_D %s/db/BCRTCR/TRBD.fa  -germline_db_J %s/db/BCRTCR/TRBJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterlostRepeatFasta,tcrbFile)
os.system(cmd)
print "Run: ",cmd


#TCRD
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/TRDV.fa -germline_db_D %s/db/BCRTCR/TRBD.fa  -germline_db_J %s/db/BCRTCR/TRDJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterlostRepeatFasta,tcrdFile)
os.system(cmd)
print "Run: ",cmd

#TCRG
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/TRGV.fa -germline_db_D %s/db/BCRTCR/TRBD.fa  -germline_db_J %s/db/BCRTCR/TRGJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterlostRepeatFasta,tcrgFile)
os.system(cmd)
print "Run: ",cmd

#######################################################################################################################################
print "*****************************Identify microbial reads**********************************************************"
#bacteria
cmd="%s/tools/blastn -task megablast -index_name %s/db/microbiome/bacteria/bacteria -use_index true -query %s -db %s/db/microbiome/bacteria/bacteria  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,afterrRNAFasta,codeDir,bacteriaFile)
print "Run :", cmd
os.system(cmd)

#virus
cmd="%s/tools/blastn -task megablast -index_name %s/db/microbiome/virus/viruses -use_index true -query %s -db %s/db/microbiome/virus/viruses  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,afterrRNAFasta,codeDir,bacteriaFile)
print "Run :", cmd
os.system(cmd)




f.close()