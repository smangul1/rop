#Note :
#1) make sure the same read is not counted multiple times to the same genome


import sys
import csv
import os
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt


import argparse
from collections import Counter


ap = argparse.ArgumentParser()
ap.add_argument('virusMegablastSample', help='virusMegablastSample')
ap.add_argument('readRength', help='')

args = ap.parse_args()


rLength=int(args.readRength)

#NR - number of reads
NR=50





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
    genomeReads[g]=set()


#get reads from each genome
with open(args.virusMegablastSample,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        genome=line[1]
        genomeReads[genome].add(line[0])

genomesList=list(set(genomesList))



#if there are more then 10 reads assigned onto the genome prepare the coordinates of LEFT and RIGHT coordinates of the covered region
for g in genomesList:
    
        k+=1
        genomeL[g]=[]
        genomeR[g]=[]

print k

dict={}


for g in genomesList:
    for r in genomeReads[g]:
        dict[r]=[]


# make dict for each read
with open(args.virusMegablastSample,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        read=line[0]
        genome=line[1]
        dict[read].append(line)







# get LEFT and RIGHT coordinates of the covered region

dict2={}
for g in genomesList:
    dict2[g]=[]


nMultimapped=0
for key, value in dict.iteritems():
    if len(value)>1:
        
        gList=set()
        for v in value:
            if v[1] not in gList:
                dict2[v[1]].append(v)
        
        
        nMultimapped+=1

    else :
        dict2[value[0][1]].append(value[0])


print "Number of multimapped reads", nMultimapped

dictGenome={}
for key,value in dict2.items():
    dictGenome[key]=[0,0,0,[],[],""]


#['HWI-ST1148:179:C4BAKACXX:5:1101:17968:1893/1', 'gi|752901102|ref|NC_026427.1|', '95.65', '92', '4', '0', '7', '98', '798', '889', '2e-35', ' 148']

for key,value in dict2.items():
    
    
    lGlobal=[]
    rGlobal=[]
    
    if len(value)>NR:
        
        for v in value:
            
            
            l=int(v[8])
            r=int(v[9])
            
            
            
            

            
            if l>r:
                l, r = r, l
            
            
            lGlobal.append(int(l))
            rGlobal.append(int(r))



        dictGenome[key][0]=len(value)
        dictGenome[key][5]=line[1]
        dictGenome[key][1]=min(lGlobal)
        dictGenome[key][2]=max(rGlobal)
        dictGenome[key][3]=[0]*(dictGenome[key][2]+1)
        dictGenome[key][4]=[0]*(dictGenome[key][2]+1)







for key,value in dict2.items():
    
    if len(value)>NR:
        for v in value:
            l=int(v[8])
            r=int(v[9])
            if l>r:
                l, r = r, l
            

            
            
            dictGenome[key][3][l]+=1
            
            for i in range(l,r+1):
                dictGenome[key][4][i]+=1

    plt.figure(1)
    plt.subplot(211)
    plt.plot( dictGenome[key][3], 'ro')


    plt.subplot(212)
    plt.plot( dictGenome[key][4], 'ro')




    plt.savefig("%s.png" %key)







