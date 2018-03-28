import sys
import csv
import os
import argparse


ap = argparse.ArgumentParser()
ap.add_argument('bed', help='File in bed format donwloaded from UCSC')
ap.add_argument('gt', help='File with gene Ids as a first column, and transcript Ids as a second one. Please see the instruction how to prepare file https://github.com/smangul1/rop/wiki/How-to-prepare-gene-annotations-for-gprofile')
ap.add_argument('gc', help='File prepare from gtf by extractGeneBoundaries.py')
ap.add_argument('out', help='clean bed fileto be created')
ap.add_argument('organism', help='h- human, m - for mouse')


#ap.add_argument('--testN', type=int,
#                help='Run a test using only the first N features, and then '
#                'print out some example feature IDs and their attributes')
#ap.add_argument('--force', action='store_true',
#                help='Overwrite an existing database')


#You can output the UTRs in BED format using the Table Browser
#(http://genome.ucsc.edu/cgi-bin/hgTables). Once you've selected your
#assembly of interest select:

#group: Genes and Gene Prediction Tracks
#track: UCSC Genes
#table: knownGene
#region: you can retrieve data for the whole genome or for a specified region

#output format: BED - browser extensible data




args = ap.parse_args()


#chr1    67208778        67210057        ENST00000237247_utr3_26_0_chr1_67208779_f



chr_list=[]


#human or mouse
if args.organism=='m':
    for i in range(1,20):
        chr_list.append(str(i))
    
    chr_list.append('X')
    chr_list.append('Y')
    chr_list.append('MT')
elif args.organism=='h':
    for i in range(1,23):
        chr_list.append(str(i))
    
    chr_list.append('X')
    chr_list.append('Y')
    chr_list.append('MT')
else:
    print "ERROR"
    sys.exit(1)



bed=[]


trIDSet=set()

#-------------------------
dictTrID_GeneID={}
#ENSMUSG00000090025,ENSMUST00000160944
with open(args.gt,'r') as f:
    reader=csv.reader(f)
    for line in reader:
        trID=line[1]
        geneID=line[0]
        dictTrID_GeneID[trID]=geneID
        trIDSet.add(trID)

#-------------------------
dictGeneID_GeneName={}
#1,non-rRNA,ENSMUSG00000000544,Gpa33,168060369,168096641
with open(args.gc,'r') as f:
    reader=csv.reader(f)
    for line in reader:
        geneID=line[2]
        geneName=line[3]
        dictGeneID_GeneName[geneID]=geneName






with open(args.bed,'r') as f:
    
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        chr=line[0].split('chr')[1]
        if chr in chr_list:
            
            x=line[1]
            y=line[2]
            trID=line[3]
            if trID in trIDSet:
                geneID=dictTrID_GeneID[trID]
                geneName=dictGeneID_GeneName[geneID]
                bed.append((chr,x,y,geneID,geneName))



print "Save to",args.out


outFile = open(args.out,'wb')
wr = csv.writer(outFile)
wr.writerows(sorted(bed))


