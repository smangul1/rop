import sys
import csv
import os
import argparse


ap = argparse.ArgumentParser()
ap.add_argument('gtf', help='File in bed format donwloaded from UCSC')
ap.add_argument('out', help='file with the gene boundaries')
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
    print "ERROR . Please choose m/h"
    sys.exit(1)


geneNames=[]
gtf=[]
gene={}
gene_rRNA={}
geneNames_rRNA=[]
stuff = dict()
dictRead=dict()
geneChromosome={}


#1       processed_transcript    exon    3195982 3197398 .       -       .       exon_number "2"; gene_biotype "protein_coding"; gene_id "ENSMUSG00000051951"; gene_name "Xkr4"; transcript_id "ENSMUST00000162897"; transcript_name "Xkr4-003"; tss_id "TSS49891";


dictIdName={}

with open(args.gtf,'r') as f:
    
    reader=csv.reader(f,delimiter='\t')
    for line in reader:
        if line[0] in chr_list:
            gtf.append(line)
            gn = line[8].split("gene_id \"")[1].split("\";")[0]
            gName=line[8].split("gene_name \"")[1].split("\";")[0]
            dictIdName[gn]=gName
            geneChromosome[gn]=line[0]
            if line[1]=="rRNA":
                geneNames_rRNA.append(gn)
            else:
                geneNames.append(gn)



geneNames=set(geneNames)
print "Number of genes",len(geneNames)


geneNames_rRNA=set(geneNames_rRNA)
print "Number of rRNA genes",len(geneNames_rRNA)
#create dictionary to acess the elements of gene by geneName




for i in geneNames:
    gene[i]=[]
    
for i in geneNames_rRNA:
    gene_rRNA[i]=[]


for i in gtf:
    gn=i[8].split("gene_id \"")[1].split("\";")[0]
    if i[2]!="transcript":
        if i[1]=="rRNA":
            gene_rRNA[gn].append(int(i[3]))
            gene_rRNA[gn].append(int(i[4]))
        else:
            gene[gn].append(int(i[3]))
            gene[gn].append(int(i[4]))


filename, file_extension = os.path.splitext(args.out)





bed=[]

for i in geneNames:
    bed.append((geneChromosome[i],'non-rRNA',i,dictIdName[i],min(gene[i]),max(gene[i])))
     
for i in geneNames_rRNA:
    bed.append((geneChromosome[i],'rRNA',i,dictIdName[i],min(gene_rRNA[i]),max(gene_rRNA[i])))

                            



print "Save to",args.out


outFile = open(args.out,'wb')
wr = csv.writer(outFile)
wr.writerows(sorted(bed))


