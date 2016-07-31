import pysam
import sys
import csv
import os
import argparse



from quicksect import IntervalNode
from random import randint, seed


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))




def is_junction(read):
    for c in read.cigartuples:
        if c[0]==3:
            return True
    return False





def whichFeature(read,chr):
    find_list_cds=find(read.reference_start, read.reference_end , tree_cds[chr])
    find_list_utr3=find(read.reference_start, read.reference_end , tree_utr3[chr])
    find_list_utr5=find(read.reference_start, read.reference_end , tree_utr5[chr])
    find_list_intron=find(read.reference_start, read.reference_end , tree_geneCoordinates[chr])
    find_list_intergenic=find(read.reference_start, read.reference_end , tree_intergenic[chr])
    

    
    threshold=len(read.query_sequence)*0.75
    
    
    tag_cds=0
    tag_utr3=0
    tag_utr5=0
    max_cds=0
    max_utr3=0
    max_utr5=0
    
    cds_overlap=[]
    utr3_overlap=[]
    utr5_overlap=[]
    
    for i in find_list_cds:
        overlap=getOverlap((read.reference_start,read.reference_end),i)
        if overlap>threshold:
            cds_overlap.append(overlap)
    
    
    for i in find_list_utr3:
        overlap=getOverlap((read.reference_start,read.reference_end),i)
        if overlap>threshold:
            utr3_overlap.append(overlap)
            
    for i in find_list_utr5:
        overlap=getOverlap((read.reference_start,read.reference_end),i)
        if overlap>threshold:
            utr5_overlap.append(overlap)
    

 
    
    if len(cds_overlap)>0:
        tag_cds=1
        max_cds=max(cds_overlap)
    if len(utr3_overlap)>0:
        tag_utr3=1
        max_cds=max(utr3_overlap)
    if len(utr5_overlap)>0:
        tag_utr5=1
        max_cds=max(utr5_overlap)

    
    
    if tag_cds>1 and tag_utr3+tag_utr5>1:
        print "-------"
        print "MIXED"
        print tag_cds,tag_utr3,tag_utr5
        print max_cds,tag_utr3,tag_utr5
        print find_list_cds
        print find_list_utr3
        print find_list_utr5
        print read
        print "-------"
        return 'MIXED'
    elif tag_utr3+tag_utr5>1:
        x=find_list_utr3[0][0]
        y=find_list_utr3[0][1]
        return ('UTR_',geneUTR3[(x,y)])
    elif tag_cds==1:
        
        x=find_list_cds[0][0]
        y=find_list_cds[0][1]
        return ('CDS',geneCDS[(x,y)])
    elif tag_utr3==1:
        x=find_list_utr3[0][0]
        y=find_list_utr3[0][1]
        return ('UTR3',geneUTR3[(x,y)])
    elif tag_utr5==1:
        x=find_list_utr5[0][0]
        y=find_list_utr5[0][1]
        return ('UTR5',geneUTR5[(x,y)])
    elif tag_cds+tag_utr3+tag_utr5==0:
        if len(find_list_intron)>0:
            x=find_list_intron[0][0]
            y=find_list_intron[0][1]
            
            
            return ('INTRON',dictGeneNames[x,y])
        elif len(find_list_intergenic)>0:
            return ('INTERGENIC',("NA,NA"))
        else:
            return ('DEEP',("NA","NA"))

        
        

#------

def find(start, end, tree):
    "Returns a list with the overlapping intervals"
    out = []
    tree.intersect( start, end, lambda x: out.append(x) )
    return [ (x.start, x.end) for x in out ]

'''
    tree = IntervalNode( 5, 20 )
    
    
    
    overlap = find(27, 28 , tree)
    if overlap==[]:
    print "----"
    
    '''


ap = argparse.ArgumentParser()
ap.add_argument('bam', help='sorted bam file with mapped reads')
ap.add_argument('out', help='file to save the number of reads per genome category')
ap.add_argument("--perCategory", help="reports the assigment for each read. A separate file per chromosome will be created",action="store_true")
ap.add_argument("--mouse", help="Use mouse genome annotations (NCBIM37). Default is human",action="store_true")
args = ap.parse_args()








##About Interval trees  https://www.biostars.org/p/99/



outDir=os.path.dirname(args.out)
if outDir=="":
    outDir=os.getcwd()

chr_list=[]


                
                
if not args.mouse:
    for i in range(1,23):
        chr_list.append(str(i))
    chr_list.append('X')
    chr_list.append('Y')

    utr3_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/human/bedPrepared/UTR3_GRCh37_prepared.bed'
    utr5_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/human/bedPrepared/UTR5_GRCh37_prepared.bed'
    cds_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/human/bedPrepared/CDS_GRCh37_prepared.bed'
    geneCoordinates_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/human/bedPrepared/geneCoordinates_GRCh37.bed'

elif args.mouse:
    for i in range(1,20):
        chr_list.append(str(i))
    chr_list.append('X')
    chr_list.append('Y')
    utr3_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/mouse/bedPrepared/UTR3_NCBIM37_prepared.bed'
    utr5_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/mouse/bedPrepared/UTR5_NCBIM37_prepared.bed'
    cds_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/mouse/bedPrepared/CDS_NCBIM37_prepared.bed'
    geneCoordinates_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/mouse/bedPrepared/geneCoordinates_NCBIM37.bed'


base=os.path.basename(args.bam)
prefix=os.path.splitext(base)[0]

dirOutPerCategory=""


#-------perCategory-------
if args.perCategory:
    dirOutPerCategory=outDir+"/"+prefix+"_perCategory/"
    if not os.path.exists(dirOutPerCategory):
        os.makedirs(dirOutPerCategory)
    
    outFile={}
    for chr in chr_list:
        f_file=dirOutPerCategory+prefix+"."+chr+".genomicFeature"
        outfile = open(f_file, 'w' )
        outFile[chr]=open(f_file, 'w' )
    #MT
    f_file=dirOutPerCategory+prefix+"."+'MT'+".genomicFeature"
    outfile = open(f_file, 'w' )
    outFile['MT']=open(f_file, 'w' )














#DATA STRUCTURE
tree_utr3={}
tree_utr5={}
tree_cds={}
tree_geneCoordinates={}
tree_rRNA={}
tree_intergenic={} # +10,000


for chr in chr_list:
    tree_utr3[chr]=IntervalNode(0,0)
    tree_utr5[chr]=IntervalNode(0,0)
    tree_cds[chr]=IntervalNode(0,0)
    tree_geneCoordinates[chr]=IntervalNode(0,0)
    tree_rRNA[chr]=IntervalNode(0,0)
    tree_intergenic[chr]=IntervalNode(0,0)










print "Load gene annotations ..."


geneUTR3={}

#UTR3
print "Load",utr3_file
with open(utr3_file,'r') as f:
    
    reader=csv.reader(f)
    for line in reader:
        chr=line[0]
        if chr in chr_list:
            x=int(line[1])
            y=int(line[2])
            geneID=line[3]
            geneName=line[4]
            geneUTR3[(x,y)]=(geneID,geneName)
            tree_utr3[chr]=tree_utr3[chr].insert( x, y )








#find_list=find(67208778, 67210057 , tree_utr3[chr])

geneUTR5={}

#UTR5
print "Load",utr5_file
with open(utr5_file,'r') as f:
    
    reader=csv.reader(f)
    for line in reader:
        chr=line[0]
        if chr in chr_list:
            x=int(line[1])
            y=int(line[2])
            geneID=line[3]
            geneName=line[4]
            geneUTR5[(x,y)]=(geneID,geneName)
            tree_utr5[chr]=tree_utr5[chr].insert( x, y )


#CDS
geneCDS={}


print "Load",cds_file
with open(cds_file,'r') as f:
    
    reader=csv.reader(f)
    for line in reader:
        chr=line[0]
        if chr in chr_list:
            x=int(line[1])
            y=int(line[2])
            geneID=line[3]
            geneName=line[4]
            geneCDS[(x,y)]=(geneID,geneName)
            tree_cds[chr]=tree_cds[chr].insert( x, y )






#gene coordinates
nGenes_non_rRNA=0
nGenes_rRNA=0


#['1', 'non-rRNA', 'ENSG00000000457', 'SCYL3', '169818772', '169863408']


print "Load",geneCoordinates_file

#1,non-rRNA,ENSMUSG00000000544,Gpa33,168060369,168096

dictGeneNames={}

with open(geneCoordinates_file,'r') as f:
    
    reader=csv.reader(f)
    for line in reader:
        chr=line[0]
        if chr in chr_list:
            x=int(line[4])
            y=int(line[5])
            geneID=line[2]
            geneName=line[3]
            
            if line[1]=='non-rRNA':
                nGenes_non_rRNA+=1
                tree_geneCoordinates[chr]=tree_geneCoordinates[chr].insert( x, y )
                dictGeneNames[x,y]=(geneID,geneName)
                
                
                x_10K=x-10000
                y_10K=y+10000
                if x_10K<0:
                    x_10K=0
                tree_intergenic[chr]=tree_intergenic[chr].insert( x_10K, y_10K )
            elif line[1]=='rRNA':
                tree_rRNA[chr]=tree_rRNA[chr].insert( x, y )
                nGenes_rRNA+=1







#
#======================================================================
#BAM

if args.perCategory:
    outFile={}
    for chr in chr_list:
        f_file=dirOutPerCategory+prefix+"."+chr+".genomicFeature"
        outfile = open(f_file, 'w' )
        outFile[chr]=open(f_file, 'w' )


    #MT
    f_file=dirOutPerCategory+prefix+"."+'MT'+".genomicFeature"
    outfile = open(f_file, 'w' )
    outFile['MT']=open(f_file, 'w' )



print "Open bam file",args.bam
bamfile = pysam.Samfile(args.bam, "rb")


#list for read categories
multiMappedReads=[]
fusionReads=[]

#counts
nrRNA=0
nDeep=0
nIntergenic=0
nIntron=0
nCDS=0
nUTR3=0
nUTR5=0
nUTR_=0
nJunction=0
nMultiMapped=0
nIntron=0
nMT=0


singleton=[]

for chr in chr_list:
    
    print "Process chr",chr
    for read in bamfile.fetch(chr):
        readName=read.query_name
             
        if read.mapq!=50:
            multiMappedReads.append(readName)
        elif is_junction(read):
            feature=whichFeature(read,chr)
            
            if args.perCategory:
                outFile[chr].write( readName+','+chr + ',' + 'junction' + ',' + feature[1][0] + ',' + feature[1][1] + '\n' )
            nJunction+=1
        else:
            feature=whichFeature(read,chr)
            if args.perCategory:
                outFile[chr].write( readName+','+chr + ',' + feature[0] + ',' + feature[1][0] + ',' + feature[1][1] + '\n' )
            if feature[0]=='CDS':
                nCDS+=1
            elif feature[0]=='INTRON':
                nIntron+=1
            elif feature[0]=='UTR3':
                nUTR3+=1
            elif feature[0]=='UTR5':
                nUTR5+=1
            elif feature[0]=='UTR_':
                 nUTR_+=1
            elif feature[0]=='INTERGENIC':
                nIntergenic+=1
            elif feature[0]=='DEEP':
                nDeep+=1
                    

for read in bamfile.fetch('MT'):
    flag_multiMapped=0
    if read.is_read1:
        readName=read.query_name
    else:
        readName=read.query_name
    
    if read.mapq!=50:
            multiMappedReads.append(readName)
    else:
        if args.perCategory:
            outFile['MT'].write( readName+','+'MT' + ',' + 'MT' + ",NA,NA"+'\n' )
        nMT+=1

multiMappedReads=set(multiMappedReads)





if args.perCategory:
    #multiMappedReads
    f_multiMappedReads=outDir+"/"+prefix+'__multiMappedReads.reads'
    outfile = open(f_multiMappedReads, 'w' )
    for i in multiMappedReads:
        outfile.write(i)
        outfile.write("\n")



nMultiMapped=len(multiMappedReads)

nTotalMapped=nrRNA+nDeep+nIntergenic+nIntron+nCDS+nUTR3+nUTR5+nUTR_+nMultiMapped+nJunction+nMT


header=[]


header.append('sampleName')
header.append('nTotalMapped')
header.append('nJunction')
header.append('nCDS')
header.append('nUTR3')
header.append('nUTR5')
header.append('nUTR_')
header.append('nIntron')
header.append('nIntergenic')
header.append('nDeep')
header.append('nMT')
header.append('nMultiMapped')


gf=[]


gf.append(prefix)
gf.append(nTotalMapped)
gf.append(nJunction)
gf.append(nCDS)
gf.append(nUTR3)
gf.append(nUTR5)
gf.append(nUTR_)
gf.append(nIntron)
gf.append(nIntergenic)
gf.append(nDeep)
gf.append(nMT)
gf.append(nMultiMapped)



c = csv.writer(open(args.out, "w"))
c.writerow(header)
c.writerow(gf)



















