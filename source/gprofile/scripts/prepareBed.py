import sys
import csv
import os
import argparse


ap = argparse.ArgumentParser()
ap.add_argument('bed', help='File in bed format donwloaded from UCSC')
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

with open(args.bed,'r') as f:
    
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        chr=line[0].split('chr')[1]
        if chr in chr_list:
            x=line[1]
            y=line[2]
            bed.append((chr,x,y))

print "Number of element initially",len(bed)
bed=set(bed)
print "Number of element after removing dublicated elements",len(bed)


print "Save to",args.out


outFile = open(args.out,'wb')
wr = csv.writer(outFile)
wr.writerows(sorted(bed))


