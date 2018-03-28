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


