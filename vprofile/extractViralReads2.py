import sys
import csv
import os

import argparse
from collections import Counter


ap = argparse.ArgumentParser()
ap.add_argument('virusMegablastSample', help='virusMegablastSample')
args = ap.parse_args()



base=os.path.basename(args.virusMegablastSample)
prefix=os.path.splitext(base)[0]


print "Read ",args.virusMegablastSample

out=prefix+"_filtered.csv"

print "Save to",out

dict={}
genomesList=[]
genomes={}



k1=0
k2=0

with open(args.virusMegablastSample,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        genomesList.append(line[1])
        genomes[line[1]]=[]
        k1+=1
        element=line[0]
        identity=float(line[2])
        alignmentLength=float(line[7])-float(line[6])
        eValue=float(line[10])
        if eValue<1e-05 and alignmentLength>=80 and identity>=90:
            k2+=1
            dict[element]=[alignmentLength,identity,eValue]

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



