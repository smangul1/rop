import pysam
import sys
import csv
import os
import argparse
import re




ap = argparse.ArgumentParser()
ap.add_argument('dir', help='directory where the tsv  are located')
ap.add_argument('sampleNames', help='File with sample names. Assumes the files have .MetaPhlan.tsv extension')
ap.add_argument('prefixOut', help='fileto save the stat')
ap.add_argument('taxa', help='taxa rank e.g. phylum')






#ap.add_argument('--testN', type=int,
#                help='Run a test using only the first N features, and then '
#                'print out some example feature IDs and their attributes')
#ap.add_argument('--force', action='store_true',
#                help='Overwrite an existing database')

#cmd https://gist.github.com/daler/ec481811a44b3aa469f3

args = ap.parse_args()



#=======================================



header=['sample','taxa','relativeFrequency']
c = csv.writer(open(args.prefixOut + "_tProfile.txt", "w"))
c.writerow(header)




#k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Oxalobacteraceae|g__Candidatus_Zinderia|s__Candidatus_Zinderia_unclassified   80.76285


alphaList=[]


taxaList=[]
richness=[]
taxonomicProfile={}


    


print "Read",args.sampleNames
with open(args.sampleNames,'r') as f:
    
    
    
    reader=csv.reader(f)
    for line in reader:
        pFile=args.dir+"/"+line[0]+'.MetaPhlan.tsv'
        
        tProfile=[]
        tProfile[:]=[]
        rFreq=[]
        rFreq[:]=[]
        alpha=0
        
        
        with open(pFile) as f2:
            reader=csv.reader(f2)
            for line2 in reader:
                
                result = re.search(args.taxa+'__([a-zA-Z_]*)\t', line2[0])
                if result:
                    
                    taxa= result.group(1)
                    taxaList.append(taxa)
                    freq=float(line2[0].split('\t')[1])
                    print taxa,freq
                    
                    
                    
                    sample=line[0].split('.')[0]
                    
                    if freq>0:
                        tProfile.append((sample,tissue,phenotype,phenotype1,phenotype2,phenotype3,taxa,freq))
        

        
        if len(tProfile)>0:
            
            tProfile = [(i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]/sum(zip(*tProfile)[7])) for i in tProfile]
            print tProfile
            
            for t in tProfile:
                c.writerow(t)
                taxonomicProfile[t[1]].append(t[2])
                rFreq.append(t[7])

            
            for r in rFreq:
                alpha+=float(r*r)
            print "alpha",alpha
            if alpha==0:
                    alphaList.append((sample,phenotype,phenotype1,phenotype2,phenotype3,0,len(tProfile)))
            else:
                alphaList.append((sample,phenotype,phenotype1,phenotype2,phenotype3,1/alpha,len(tProfile)))

            richness.append(len(tProfile))







print "Save alpha to ", args.prefixOut+"_summary.txt"

c = csv.writer(open(args.prefixOut+"_summary.txt", "w"))
header=['sample','alpha','richness']
c.writerow(header)
for a in alphaList:
    c.writerow(a)






import pandas as pd
import numpy as np

df = pd.read_csv(args.prefixOut + "_tProfile.txt")

df.groupby(['taxa']).mean().to_csv(args.prefixOut + "_tProfile_mean.txt")
df.groupby(['taxa']).std().to_csv(args.prefixOut + "_tProfile_std.txt", mode='a')

print "richness",np.mean(richness),"+-",np.std(richness)

sys.exit(1)

taxaList=set(taxaList)
for t in taxaList:
    print t,taxonomicProfile['CSF'].count(t), taxonomicProfile['ALS'].count(t),taxonomicProfile['BPD'].count(t),taxonomicProfile['SCZ'].count(t)



print "Gamma diversity"

print len(set(taxonomicProfile['CSF'])),len(set(taxonomicProfile['ALS'])),len(set(taxonomicProfile['BPD'])),len(set(taxonomicProfile['SCZ']))


import matplotlib as mpl
