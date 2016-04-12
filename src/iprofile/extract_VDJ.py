import sys
import csv
import os
import argparse


ap = argparse.ArgumentParser()
ap.add_argument('file', help='file from IgBlast')
ap.add_argument('dir', help='dir to save results')
ap.add_argument('chainType', help='type of the chain to use : IGH, IGK, ...')
ap.add_argument('eValue', help='e value threshold, 1e-20 is recomended')
args = ap.parse_args()



try:
    os.stat(args.dir)
except:
    os.mkdir(args.dir)

base=os.path.basename(args.file)
prefix=os.path.splitext(base)[0]

#J       reversed|HWI-D00108:215:C41Y3ACXX:3:2113:4836:9536/1    IGHJ2P*01       95.24   21      1       0       0       57      77      21      41      3e-06   34.2


currentDir=os.path.dirname(os.path.realpath(__file__))

ighvFile=currentDir+"/../../db/BCRTCR/gnames/%sV.gNames" %(args.chainType)
ighjFile=currentDir+"/../../db/BCRTCR/gnames/%sJ.gNames" %(args.chainType)

ighvNames = [line.strip() for line in open(ighvFile, 'r')]
ighjNames = [line.strip() for line in open(ighjFile, 'r')]


dict={}

reads=[]


with open(args.file,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        read=line[1].replace("reversed|", "")
        if len(line)==14: #ignore lines with problem in format
            geneAllele=line[2]
            if geneAllele in ighvNames or geneAllele in ighjNames:
                reads.append(read)




reads=set(reads)
readsV=[]
readsD=[]
readsJ=[]




for r in reads:
    dict[r]=[0,(0,'',0,0,0.0,0.0),(0,'',0,0,0.0,0.0),(0,'',0,0,0.0,0.0)]

print "Number of immune reads detected by IgBlast", len(reads)

print 'Open', args.file

with open(args.file,'r') as f:
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        if len(line)==14:
            read=line[1].replace("reversed|", "")
            geneAllele=line[2]
            identity=line[3]
            eValue=float(line[12])
        
            if geneAllele in ighvNames or geneAllele in ighjNames:
                x=line[8]
                y=line[9]


                if geneAllele in ighvNames and dict[read][1][0]==0 and eValue<=float(args.eValue):
                    readsV.append(read)
                    dict[read][0]=1
                    dict[read][1]=(1,geneAllele,x,y,identity, eValue)
            
                    
                elif geneAllele in ighjNames and dict[read][3][0]==0: 
                    readsJ.append(read)
                    dict[read][0]=1
                    dict[read][3]=(1,geneAllele,x,y,identity, eValue)
   

comb=[]




for key,val in dict.items():
    if val[0]==1:
        s=val[1][0]+val[2][0]+val[3][0]
        if s>=2:
            if val[1][0]+val[3][0]==2:
                comb.append(val[1][1]+";"+val[3][1])







from collections import Counter

outFile=args.dir+"/"+prefix+".%srecomb" %(args.chainType)
print "Save to ",outFile

out = open(outFile, "w")

header="recombination,relativeAbundance"
out.write(header+"\n")



combCounter=Counter(comb)



for key,val in combCounter.items():
    out.write(key+","+str(val)+"\n")




out.close()

