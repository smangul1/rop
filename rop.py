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
ap.add_argument('dir', help='directory to save results')
args = ap.parse_args()


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

if not os.path.exists(QCDir):
    os.makedirs(QCDir)
if not os.path.exists(lostHumanDir):
    os.makedirs(lostHumanDir)
if not os.path.exists(lostRepeatDir):
    os.makedirs(lostRepeatDir)

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


runFile=args.dir+"/commands_"+basename+".sh"





f = open(runFile,'w')

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



print "*****************************Identify lost human reads******************************"





os.system(". /u/local/Modules/default/init/modules.sh \n")
os.system("module load bowtie2/2.1.0 \n")
os.system("module load samtools \n")


#genome
cmd="bowtie2 -k 1 --very-sensitive -p 8 -f -x %s/db/human/Bowtie2Index/genome -U %s | samtools view -bSF4 - | samtools sort -  %s" %(codeDir, afterrRNAFasta,os.path.splitext(gBamFile)[0])
os.system(cmd)
print "Run: ",cmd
os.system("samtools index %s \n" %gBamFile)

#transcriptome
cmd="bowtie2 -k 1 --very-sensitive -f -p 8 -x %s/db/human/Bowtie2Index/hg19KnownGene.exon_polya200 -U %s | samtools view -bSF4 - | samtools sort -  %s" %(codeDir, afterrRNAFasta,os.path.splitext(tBamFile)[0])
os.system(cmd)
print "Run: ",cmd
os.system("samtools index %s \n" %tBamFile)


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

print "*****************************Identify lost repeat reads******************************"
cmd="%s/tools/blastn -task megablast -index_name %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa -use_index true -query %s -db %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,afterlostHumanFasta,codeDir,repeatFile)
print "Run :", cmd
os.system(cmd)





f.close()