import sys
import csv
import os
import argparse
from Bio import SeqIO # module needed for sequence input




####################################################################
excludeReadsFromFasta(in,reads,out):

    fasta_sequences = SeqIO.parse(open(in),'fasta')
    with open(out, "w") as f:
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


#intermediate files
lowQFileFasta=QCDir+basename+"_lowQ.fa"
lowQCFile=QCDir+basename+"_lowQC.fa"
rRNAFile=QCDir+basename+"_rRNA_blastFormat6.csv"
afterrRNAFasta=QCDir+basename+"_after_rRNA.fasta"

runFile=args.dir+"/commands_"+basename+".sh"


#analysis directories
QCDir=args.dir+"/QC/"
lostHumanDir=args.dir+"/lostHuman/"
if not os.path.exists(QCDir):
    os.makedirs(QCDir)
if not os.path.exists(lostHumanDir):
    os.makedirs(lostHumanDir)


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



f.close()



sys.exit(1)

#bowtie2

bowtie2Dir=args.dir+"/lostHuman/"
if not os.path.exists(bowtie2Dir):
    os.makedirs(bowtie2Dir)
runBowtie2=args.dir+"/runBowtie2_"+basename+".sh"
gBamFile=bowtie2Dir+basename+"_genome.bam"
tBamFile=bowtie2Dir+basename+"_transcriptome.bam"



f = open(runBowtie2,'w')
f.write(". /u/local/Modules/default/init/modules.sh \n")
f.write("module load bowtie2/2.1.0 \n")
f.write("module load samtools \n")
f.write("bowtie2 -k 10 --very-sensitive -p 8 -x ~/project/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome -U %s | samtools view -bSF4 - | samtools sort -  %s \n" %(lowQFile,os.path.splitext(gBamFile)[0]) )
f.write("samtools index %s \n" %gBamFile)
f.write("bowtie2 -k 10 --very-sensitive -p 8 -x /u/home/s/serghei/project/Homo_sapiens/hg19KnownIsoforms/Bowtie2Index/hg19KnownGene.exon_polya200 -U %s | samtools view -bSF4 -  | samtools sort -  %s \n" %(lowQFile,os.path.splitext(tBamFile)[0]) )
f.write("samtools index %s \n" %tBamFile)

f.close()


