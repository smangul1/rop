import sys
import csv
import os
import argparse




#codeDir
codeDir=os.path.dirname(os.path.realpath(__file__))



sys.path.append('%s/tools/biopython/biopython-1.66/' %(codeDir))


import Bio
from Bio import SeqIO # module needed for sequence input

sys.path.append('/u/home/s/serghei/project/code/import/pysam-master/')




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


####################################################################
def bam2fasta(codeDir,inFile,outFile):
    print "*****************************Convert unmapped bam to fasta ******************************"
    cmdConvertBam2Fastq="%s/tools/bamtools convert -in %s -format fasta >%s" %(codeDir,inFile,outFile)
    print "Run:",cmdConvertBam2Fastq
    os.system(cmdConvertBam2Fastq)

####################################################################
def bam2fastq(codeDir,inFile,outFile):
    print "*****************************Convert unmapped bam to fastq ******************************"
    #os.system(". /u/local/Modules/default/init/modules.sh")
    #os.system("module load bamtools")
    cmdConvertBam2Fastq="%s/tools/bamtools convert -in %s -format fastq >%s" %(codeDir,inFile,outFile)
    print "Run:",cmdConvertBam2Fastq
    os.system(cmdConvertBam2Fastq)




ap = argparse.ArgumentParser('python rop.py')
ap.add_argument('unmappedReads', help='unmapped Reads in the fastq format')
ap.add_argument('dir', help='directory (absolute path) to save results of the analysis')

ap.add_argument("--qsub", help="submit qsub jobs on hoffman2 cluster",
                    action="store_true")
ap.add_argument("--qsubArray", help="prepare qsub scripts to be run latetr using job array",
                action="store_true")
ap.add_argument("--b", help="unmapped reads in bam format",
                action="store_true")
ap.add_argument("--skipLowq", help="skip step filtering ",
                action="store_true")
ap.add_argument("--skipQC", help="skip entire QC step : filtering  low-quality, low-complexity and rRNA reads (reads mathing rRNA repeat unit)",
                action="store_true")
ap.add_argument("--NCL_CIRI", help="enable CIRI for non-co-linear RNA sequence analysis", action="store_true")
ap.add_argument("--immune", help = "Only TCR/BCR immune gene analysis will be performed", action = "store_true")



args = ap.parse_args()



#move to halp message
print "For more details see:"
print "<https://github.com/smangul1/rop>"
print "<https://github.com/smangul1/rop/ROPmanual.pdf>"








#basename
basename=os.path.splitext(os.path.basename(args.unmappedReads))[0]





#analysis directories
QCDir=args.dir+"/QC/"
humanDir=args.dir+"/human/"

lostHumanDir=humanDir+"/lostHuman/"
lostRepeatDir=humanDir+"/lostRepeat/"
bcrDir=args.dir+"/BCR/"
tcrDir=args.dir+"/TCR/"

NCL_CIRI_Dir=args.dir+"/NCL_CIRI/" 

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
if not os.path.exists(NCL_CIRI_Dir):
    os.makedirs(NCL_CIRI_Dir)
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
afterlostHumanFasta=lostHumanDir+basename+"_after_lostHuman.fasta"
gBamFile=lostHumanDir+basename+"_genome.bam"
tBamFile=lostHumanDir+basename+"_transcriptome.bam"
repeatFile=lostRepeatDir+basename+"_lostRepeats_blastFormat6.csv"
afterlostRepeatFasta=lostRepeatDir+basename+"_after_lostRepeat.fasta"
NCL_CIRI_file=NCL_CIRI_Dir + basename + "_NCL_CIRI_after_bwa.sam"
after_NCL_CIRI_file_prefix = basename + "NCL_CIRI_AFTER"
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
runLostHumanFile=lostHumanDir+"/runLostHuman_"+basename+".sh"
runLostRepeatFile=lostRepeatDir+"/runLostRepeat_"+basename+".sh"
runNCL_CIRIfile = NCL_CIRI_Dir + "/run_NCL_CIRI" + basename + ".sh" 
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
        print "Number of unmapped reads",n




        #lowQ
        print "*****************************Running FASTX to filter low quality reads******************************"
        cmd=codeDir+"/tools/fastq_quality_filter -v -Q 33 -q 20 -p 75 -i %s -o %s \n" %(unmappedFastq,lowQFile)
        print "Run ", cmd
        os.system(cmd)
        if args.b:
            os.remove(unmappedFastq)
        print "Save reads after filtering low quality reads to ", lowQFile


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
    cmd="export PATH=$PATH:%s/tools/seqclean-x86_64/bin" %(codeDir)
    print "Run ", cmd
    os.system(cmd)

    cmd=codeDir+"/tools/seqclean-x86_64/seqclean %s -l 50 -M -o %s" %(lowQFileFasta, lowQCFile)
    print "Run ", cmd
    os.system(cmd)

    cmd="export PATH=$PATH:%s/tools/seqclean-x86_64/bin" %(codeDir)


    print "Save reads after filtering low complexity (e.g. ACACACAC...) reads to ", lowQCFile

    #megablast rRNA
    print "*****************************Identify reads from RNA repeat unit******************************"
    cmd="%s/tools/blastn -task megablast -index_name %s/db/rRNA/rRNA -use_index true -query %s -db %s/db/rRNA/rRNA  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,lowQCFile,codeDir,rRNAFile)
    #print "Run :", cmd
    os.system(cmd)


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


    os.remove(lowQFile)
    os.remove(lowQCFile)
    os.remove(lowQFileFasta)



#######################################################################################################################################
print "*****************************Identify lost human reads******************************"





os.system(". /u/local/Modules/default/init/modules.sh \n")
os.system("module load bowtie2/2.1.0 \n")
os.system("module load samtools \n")


#genome
cmdGenome="bowtie2 -k 1 --very-sensitive -p 8 -f -x %s/db/human/Bowtie2Index/genome -U %s | samtools view -bSF4 - | samtools sort -  %s" %(codeDir, afterrRNAFasta,os.path.splitext(gBamFile)[0])
#print "Run: ",cmdGenome
#transcriptome
cmdTranscriptome="bowtie2 -k 1 --very-sensitive -f -p 8 -x %s/db/human/Bowtie2Index/hg19KnownGene.exon_polya200 -U %s | samtools view -bSF4 - | samtools sort -  %s" %(codeDir, afterrRNAFasta,os.path.splitext(tBamFile)[0])
#print "Run: ",cmdGenome

if args.qsub or args.qsubArray:
    f = open(runLostHumanFile,'w')
    f.write(". /u/local/Modules/default/init/modules.sh \n")
    f.write("module load bowtie2/2.1.0 \n")
    f.write("module load samtools \n")
    f.write(cmdGenome+"\n")
    f.write(cmdTranscriptome+"\n")
    f.write("samtools index %s \n" %gBamFile)
    f.write("samtools index %s \n" %tBamFile)
    f.write("echo \"done!\">%s/%s_lostHuman.done \n" %(lostHumanDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N lostHuman -l h_data=16G,time=10:00:00 %s" %(runLostHumanFile)
        os.system(cmdQsub)

else:
    os.system(cmdGenome)
    os.system(cmdTranscriptome)
    os.system("samtools index %s \n" %gBamFile)




if not args.qsub and not args.qsubArray:

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
cmd="%s/tools/blastn -task megablast -index_name %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa -use_index true -query %s -db %s/db/repeats/human_repbase_20_07/human_repbase_20_07.fa  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,afterrRNAFasta,codeDir,repeatFile)
#print "Run :", cmd

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


if not args.qsub and not args.qsubArray:

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

if args.NCL_CIRI and not args.immune:
    print "*****************************Identify NCL events******************************"
    cmd="%s/tools/bwa mem -T -S %s/db/human/BWAIndex/genome.fa %s > %s \n" %(codeDir,codeDir,afterrRNAFasta,NCL_CIRI_file)
    cmd = cmd + "pearl %s/tools/CIRI_v1.2.pl -S -I %s -O %s -F %s/db/human/BWAIndex/genome.fa" %(codeDir,NCL_CIRI_file,after_NCL_CIRI_file_prefix,codeDir)
    if args.qsub or args.qsubArray:
        f = open(runNCL_CIRIfile,'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_NCL_CIRI.done \n" %(NCL_CIRI_Dir,basename))
        f.close()
        if args.qsub:
            cmdQsub="qsub -cwd -V -N NCL_CIRI -l h_data=8G,time=10:00:00 %s" %(runNCL_CIRIfile)
            os.system(cmdQsub)
    else:
        os.system(cmd)


# if not args.qsub and not args.qsubArray:

#     lostRepeatReads = set()

#     with open(repeatFile,'r') as f:
#         reader=csv.reader(f,delimiter='\t')
#         for line in reader:
#             element=line[0]
#             identity=float(line[2])
#             alignmentLength=float(line[3])
#             eValue=float(line[10])
#             if eValue<1e-05 and alignmentLength>=80 and identity>=90:
#                 lostRepeatReads.add(element)
"""
Should we exclude??? FIX IT
"""
    # excludeReadsFromFasta(afterlostHumanFasta,lostRepeatReads,afterlostRepeatFasta)















#######################################################################################################################################

print "*****************************Identify VDJ recombinations from BCR and TCR******************************"

#IGH
cmd="ln -s %s//db/BCRTCR/internal_data/ %s" %(codeDir,ighDir)
os.system(cmd)

cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/IGHV.fa -germline_db_D %s/db/BCRTCR/IGHD.fa  -germline_db_J %s/db/BCRTCR/IGHJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"D\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterrRNAFasta,ighFile)
#print "Run: ",cmd
            
if args.qsub or args.qsubArray:
    f = open(runIGHFile,'w')
    f.write("ln -s %s//db/BCRTCR/internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\"> %s/%s_igh.done \n" %(ighDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N igh -l h_data=16G,time=24:00:00 %s" %(runIGHFile)
        os.system(cmdQsub)
else:
    os.system(cmd)
                    

            
#IGK
cmd="ln -s %s//db/BCRTCR/internal_data/ %s" %(codeDir,igkDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/IGKV.fa -germline_db_D %s/db/BCRTCR/IGHD.fa  -germline_db_J %s/db/BCRTCR/IGKJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterrRNAFasta,igkFile)
#print "Run: ",cmd
            
if args.qsub or args.qsubArray:
    f = open(runIGKFile,'w')
    f.write("ln -s %s//db/BCRTCR/internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\"> %s/%s_igk.done \n" %(igkDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N igk -l h_data=16G,time=24:00:00 %s" %(runIGKFile)
        os.system(cmdQsub)
else:
    os.system(cmd)
            
            
#IGL
cmd="ln -s %s//db/BCRTCR/internal_data/ %s" %(codeDir,iglDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/IGLV.fa -germline_db_D %s/db/BCRTCR/IGHD.fa  -germline_db_J %s/db/BCRTCR/IGLJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterrRNAFasta,iglFile)
#print "Run: ",cmd
if args.qsub or args.qsubArray:
    f = open(runIGLFile,'w')
    f.write("ln -s %s//db/BCRTCR/internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_igl.done \n" %(iglDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N igl -l h_data=16G,time=24:00:00 %s" %(runIGLFile)
        os.system(cmdQsub)
else:
    os.system(cmd)
            
#TCRA
cmd="ln -s %s//db/BCRTCR/internal_data/ %s" %(codeDir,tcraDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/TRAV.fa -germline_db_D %s/db/BCRTCR/TRBD.fa  -germline_db_J %s/db/BCRTCR/TRAJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterrRNAFasta,tcraFile)
#print "Run: ",cmd
if args.qsub or args.qsubArray:
    f = open(runTCRAFile,'w')
    f.write("ln -s %s//db/BCRTCR/internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_tcra.done \n"%(tcraDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N tcra -l h_data=16G,time=24:00:00 %s" %(runTCRAFile)
        os.system(cmdQsub)
else:
    os.system(cmd)
            
            

#TCRB
cmd="ln -s %s//db/BCRTCR/internal_data/ %s" %(codeDir,tcrbDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/TRBV.fa -germline_db_D %s/db/BCRTCR/TRBD.fa  -germline_db_J %s/db/BCRTCR/TRBJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterrRNAFasta,tcrbFile)
#print "Run: ",cmd
if args.qsub or args.qsubArray:
    f = open(runTCRBFile,'w')
    f.write("ln -s %s//db/BCRTCR/internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_tcrb.done \n"%(tcrbDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N tcrb -l h_data=16G,time=24:00:00 %s" %(runTCRBFile)
        os.system(cmdQsub)
else:
    os.system(cmd)
            
            

#TCRD
cmd="ln -s %s//db/BCRTCR/internal_data/ %s" %(codeDir,tcrdDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/TRDV.fa -germline_db_D %s/db/BCRTCR/TRBD.fa  -germline_db_J %s/db/BCRTCR/TRDJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterrRNAFasta,tcrdFile)
#print "Run: ",cmd
if args.qsub or args.qsubArray:
    f = open(runTCRDFile,'w')
    f.write("ln -s %s//db/BCRTCR/internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_tcrd.done \n" %(tcrdDir, basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N tcrd -l h_data=16G,time=24:00:00 %s" %(runTCRDFile)
        os.system(cmdQsub)
else:
    os.system(cmd)
            
            
#TCRG
cmd="ln -s %s//db/BCRTCR/internal_data/ %s" %(codeDir,tcrgDir)
os.system(cmd)
cmd="%s/tools/igblastn -germline_db_V %s/db/BCRTCR/TRGV.fa -germline_db_D %s/db/BCRTCR/TRBD.fa  -germline_db_J %s/db/BCRTCR/TRGJ.fa -query %s -outfmt 7 -evalue 1e-05  | awk '{if($13<1e-05 && ($1==\"V\" || $1==\"J\")) print }' >%s" %(codeDir,codeDir,codeDir,codeDir,afterrRNAFasta,tcrgFile)
#print "Run: ",cmd
            
if args.qsub or args.qsubArray:
    f = open(runTCRGFile,'w')
    f.write("ln -s %s//db/BCRTCR/internal_data/ ./ \n" %(codeDir))
    f.write(cmd+"\n")
    f.write("echo \"done!\">%s/%s_tcrg.done \n" %(tcrgDir,basename))
    f.close()
    if args.qsub:
        cmdQsub="qsub -cwd -V -N tcrg -l h_data=16G,time=24:00:00 %s" %(runTCRGFile)
        os.system(cmdQsub)
else:
    os.system(cmd)
            
            
#######################################################################################################################################
if not args.immune:
    print "*****************************Identify microbial reads**********************************************************"
    #bacteria
    cmd="%s/tools/blastn -task megablast -index_name %s/db/microbiome/bacteria/bacteria -use_index true -query %s -db %s/db/microbiome/bacteria/bacteria  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,afterrRNAFasta,codeDir,bacteriaFile)
    #print "Run :", cmd
    if args.qsub or args.qsubArray:
        f = open(runBacteriaFile,'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_bacteria.done \n"%(bacteriaDir, basename))
        f.close()
        if args.qsub:
            cmdQsub="qsub -cwd -V -N bacteria -l h_data=16G,time=24:00:00 %s" %(runBacteriaFile)
            os.system(cmdQsub)
    else:
        os.system(cmd)


    #virus
    cmd="%s/tools/blastn -task megablast -index_name %s/db/microbiome/virus/viruses -use_index true -query %s -db %s/db/microbiome/virus/viruses  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,afterrRNAFasta,codeDir,virusFile)
    #print "Run :", cmd
    if args.qsub or args.qsubArray:
        f = open(runVirusFile,'w')
        f.write(cmd+"\n")
        f.write("echo \"done!\">%s/%s_virus.done \n" %(virusDir, basename))
        f.close()
        if args.qsub:
            cmdQsub="qsub -cwd -V -N virus -l h_data=16G,time=24:00:00 %s" %(runVirusFile)
            os.system(cmdQsub)
    else:
        os.system(cmd)


    #http://eupathdb.org/eupathdb/
    #eukaryotic pathogens

    dbList=["ameoba",
            "crypto",
            "fungi",
            "giardia",
            "microsporidia",
            "piroplasma",
            "plasmo",
            "toxo",
            "trich",
            "tritryp"]

    print dbList


    for db in dbList:
        print db
        eupathdbFile=eupathdbDir+basename+"_"+db+"_blastFormat6.csv"
        runEupathdbFile=eupathdbDir+"/run_"+basename+"_"+db+".sh"


        cmd="%s/tools/blastn -task megablast -index_name %s/db/microbiome/eupathdb/%s -use_index true -query %s -db %s/db/microbiome/eupathdb/%s  -outfmt 6 -evalue 1e-05 -max_target_seqs 1 >%s" %(codeDir,codeDir,db,afterrRNAFasta,codeDir,db,eupathdbFile)
        ##print "Run :", cmd
        if args.qsub or args.qsubArray:
            f = open(runEupathdbFile,'w')
            f.write(cmd+"\n")
            f.write("echo \"done!\">%s/%s.done" %(eupathdbDir, db)+ "\n")
            f.close()
            if args.qsub:
                cmdQsub="qsub -cwd -V -N %s -l h_data=16G,time=24:00:00 %s" %(db,runEupathdbFile)
                os.system(cmdQsub)
        else:
            os.system(cmd)






