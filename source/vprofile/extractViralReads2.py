import sys
import csv
import os

import argparse
from collections import Counter


ap = argparse.ArgumentParser()
ap.add_argument('virusMegablastSample', help='virusMegablastSample')
ap.add_argument('readRength', help='')

args = ap.parse_args()


rLength=int(args.readRength)







base=os.path.basename(args.virusMegablastSample)
prefix=os.path.splitext(base)[0]


print "Read ",args.virusMegablastSample

out=prefix+"_filtered.csv"

#print "Save to",out

dict={}
genomesList=[]
genomesList2=set()
genomeReads={}
genomeReads2={}


genomes={}

genomeL={}
genomeR={}

genomeCoverage={}
genomeStartPosition={}
genomeNumberReads={}

k=0

#HWI-ST1148:179:C4BAKACXX:5:1101:11944:1949/1    gi|8486122|ref|NC_002016.1|     94.68   94      5       0       7 100     750     843     7e-35    147


#get list with genomes
with open(args.virusMegablastSample,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        genome=line[1]
        genomesList.append(genome)

for g in genomesList:
    genomeReads[g]=0


#get reads from each genome
with open(args.virusMegablastSample,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        genome=line[1]
        genomeReads[genome]+=1

genomesList=list(set(genomesList))



#if there are more then 10 reads assigned onto the genome prepare the coordinates of LEFT and RIGHT coordinates of the covered region
for g in genomesList:
    if genomeReads[g]>100:
        genomesList2.add(g)
        k+=1
        genomeL[g]=-1
        genomeR[g]=-1

print k



# get LEFT and RIGHT coordinates of the covered region
with open(args.virusMegablastSample,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        genome=line[1]
        left=int(line[8])
        right=int(line[9])
        
        if right<left:
            temp=left
            left=right
            right=left
    
        if genome in genomesList2:
            if genomeL[genome]==-1 or left<genomeL[genome]:
                genomeL[genome]=left-2
            
            if genomeR[genome]==-1 or right>genomeR[genome]:
                genomeR[genome]=right+2











import array

for g in genomesList2:
    
    genomeCoverage[g] =[0] * genomeR[g]
    
    
    #array.array('i',(0,)*genomeR[g])
    genomeStartPosition[g]=[0] * genomeR[g]
    genomeNumberReads[g]=0
    
    





print "---------"

with open(args.virusMegablastSample,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        genome=line[1]
        left=int(line[8])
        right=int(line[9])
        if right<left:
            temp=left
            left=right
            right=left
        
        if genome in genomesList2:
            #print genome,genomeL[genome],genomeR[genome]
            
            genomeNumberReads[genome]+=1
            genomeStartPosition[genome][left]+=1
            
            
            for i in range(left,right):
                genomeCoverage[genome][i]+=1


for g in genomesList2:
    
    print g,genomeNumberReads[g],genomeR[g]-genomeL[g]




g="gi|548558394|ref|NC_022518.1|"
l=int(genomeL[g])
r=int(genomeR[g])
print r-l
#print genomeCoverage[g][l:r]
#print genomeStartPosition[g][l:r]


length=r-l
nPosCovered=length-genomeCoverage[g][l:r].count(0)
nPosStart=length-genomeStartPosition[g][l:r].count(0)

print nPosStart,nPosCovered
print nPosStart/float(nPosCovered)


print genomeCoverage[g][l:r]
print "------"
print genomeStartPosition[g][l:r]

from itertools import groupby
a=genomeCoverage[g][l:r]
b = range(len(a))
for group in groupby(iter(b), lambda x: a[x]):
    if group[0]==0:
        lis=list(group[1])
        if max(lis)-min(lis)>rLength/2:
            print max(lis)-min(lis)



sys.exit(1)

#print genomeCoverage["gi|109255272|ref|NC_008168.1|"]


print "=================================="

for g in genomesList2:
    print g,genomeL[g],genomeR[g]
    print "len(genomeCoverage[g])",len(genomeCoverage[g])
    print "Number of reads", genomeNumberReads[g]
    
    
    
    for i in range(0,len(genomeCoverage[g])):
        print i, genomeCoverage[g][i], genomeStartPosition[g][i]
    sys.exit(1)







sys.exit(1)

f=open(out,'w')


print k1,k2

for key, value in dict.iteritems():
    f.write(key+","+str(value[0])+","+str(value[1])+","+str(value[2]))
    f.write("\n")
f.close()


#-------------------------------------------
#viral profiling

#genomes



#H0A9DADXX130405:2:2202:19273:35203/2    gi|548558394|ref|NC_022518.1|   95.89   73      3       0       1       73      248     176     1e-26    119


with open(args.virusMegablastSample,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        genomes[line[1]].append(line)
        element=line[0]
        identity=float(line[2])
        alignmentLength=float(line[7])-float(line[6])
        eValue=float(line[10])
        genomesList.append(line[1])
        dict[element]=[alignmentLength,identity,eValue]



for key,value in genomes.iteritems():
    print key,len(value)


print "Number of genomes",len(genomesList)



