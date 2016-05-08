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
ap.add_argument('bam', help='sorted bam file')
ap.add_argument('outPrefix', help='file to save number of reads per genome category')
args = ap.parse_args()









print os.path.dirname(os.path.realpath(__file__))


chr_list=[]


for i in range(1,23):
    chr_list.append(str(i))

chr_list.append('X')
chr_list.append('Y')






repeat_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/human/bedPrepared/hg19_rmsk_TE_prepared_noCDS.bed'
    


base=os.path.basename(args.bam)
prefix=os.path.splitext(base)[0]








#DATA STRUCTURE
tree_repeat={}



for chr in chr_list:
    tree_repeat[chr]=IntervalNode(0,0)




print "Load repeat annotations from ",repeat_file



dictClass={}
dictFamily={}
dictGene={}

dictClassCount={}
dictFamilyCount={}
dictGeneCount={}


class_list=[]
family_list=[]

#1,AluSp,Alu,SINE,16777161,16777470


k=0
with open(repeat_file,'r') as f:
    
    reader=csv.reader(f)
    for line in reader:
        chr=line[0]
        if chr in chr_list:
            Class=line[3]
            Family=line[3]+"____"+line[2]
            Gene=line[2]+"____"+line[3]+"____"+line[1]
            
            x=int(line[4])
            y=int(line[5])
            tree_repeat[chr]=tree_repeat[chr].insert( x, y )
            
            dictClass[(x,y,chr)]=Class
            dictFamily[(x,y,chr)]=Family
            dictGene[x,y,chr]=Gene
            dictClassCount[Class]=[0]
            dictFamilyCount[Family]=[0]
            dictGeneCount[Gene]=[0]

        k+=1
        if k%100000==0:
            print k, "repeat elements were loaded "
	    





#
#======================================================================
#BAM



kOverlapClass=0
kOverlapFamily=0
kOverlapGene=0


print "Open bam file",args.bam
bamfile = pysam.Samfile(args.bam, "rb")



for chr in chr_list:
    print "Process chr",chr
    for read in bamfile.fetch(chr):
        readName=read.query_name
        
        if read.mapq==50 and not is_junction(read):
            #feature=whichFeature(read,chr)
            #outFile[chr].write( readName+','+chr + ',' + feature + '\n' )
            find_list_repeat=find(read.reference_start, read.reference_end , tree_repeat[chr])
           
            if len(find_list_repeat)>1:
                overlap_list=[]
                overlap_list[:]=[]
                for f in find_list_repeat:
                    overlap=getOverlap((read.reference_start,read.reference_end),f)
                    if overlap==len(read.query_sequence):
                        overlap_list.append(overlap)
                    
                if len(overlap_list)>1:
                    print find_list_repeat
                    
                    classT=[]
                    familyT=[]
                    for f in find_list_repeat:
                        classT.append(dictClass[((f[0],f[1],chr))])
                        familyT.append(dictFamily[((f[0],f[1],chr))])
                        print "overlapping-elements",dictGene[((f[0],f[1],chr))]
                    classT=set(classT)
                    familyT=set(familyT)
                    kOverlapGene+=1
                    
                    if (len(classT)==1):
                        print classT
                        #dictClassCount[classT[0]][0]+=1
                    else:
                        kOverlapClass+=1
                
                    if (len(familyT)==1):
                        print familyT
                        #dictClassCount[familyT[0]][0]+=1
                    else:
                        kOverlapFamily+=1
                    
                            
                       
            
            if len(find_list_repeat)!=0:
                for f in find_list_repeat:
                    overlap=getOverlap((read.reference_start,read.reference_end),f)
                    if overlap==len(read.query_sequence):
                        #print overlap
                        #print read
                        #print find_list_repeat
                        Class=dictClass[((f[0],f[1],chr))]
                        Family=dictFamily[((f[0],f[1],chr))]
                        Gene=dictGene[((f[0],f[1],chr))]
                        
                        
                        #Class
                        dictClassCount[Class][0]+=1
                        
                        #Family
                        dictFamilyCount[Family][0]+=1
                         
                        #Gene
                        dictGeneCount[Gene][0]+=1
                        
                        
                        



#Class
class_names=['sample']
class_counts=[prefix]


for key,val in dictClassCount.items():
        class_names.append(key)
        class_counts.append(val[0])
    

c = csv.writer(open(args.outPrefix+"repeatClass.csv", "w"))
c.writerow(class_names)
c.writerow(class_counts)


#Family
family_names=['sample']
family_counts=[prefix]




for key,val in dictFamilyCount.items():
    family_names.append(key)
    family_counts.append(val[0])


c = csv.writer(open(args.outPrefix+"repeatFamily.csv", "w"))
c.writerow(family_names)
c.writerow(family_counts)


#Gene

gene_names=['sample']
gene_counts=[prefix]

gene_names=['sample']
gene_counts=[prefix]
for key,val in dictGeneCount.items():
    gene_names.append(key)
    gene_counts.append(val[0])


c = csv.writer(open(args.outPrefix+"repeatGene.csv", "w"))
c.writerow(gene_names)
c.writerow(gene_counts)




print 'kOverlapClass', 'kOverlapFamily','kOverlapGene'
print kOverlapClass, kOverlapFamily,kOverlapGene


print "DONE!"






