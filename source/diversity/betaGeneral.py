import sys
import csv
import os
import argparse
import collections


ap = argparse.ArgumentParser()
ap.add_argument('dir', help='dir')
ap.add_argument('outDir', help='dir to save teh results')
ap.add_argument('extension', help='extension of the individulas files in the dir')
ap.add_argument('column1', help='number of column1 with items')
ap.add_argument('column2', help='number of column2 with relative frequency')
ap.add_argument('header', help='0/1 if the header is in the individual file')


args = ap.parse_args()


# (A+B-2*J)/(A+B)

def Soerenson(dict1,dict2):
    shared=set(dict1.keys()) & set(dict2.keys())
    return 1-len(shared)*2.0/(len(set(dict1.keys()))+len(set(dict2.keys())))




def BrayCurtis(dict1,dict2):
    # Given a hash { 'species': count } for each sample, returns Compositional dissimilarity
    #BrayCurtis({'a': 0.2, 'b': 0.5, 'c': 0.3,},{'a': .4, 'b': 0.1, 'c':0.5 ,})
    #0.4
    #    # confirmed by vegan R package
    #http://www.inside-r.org/packages/cran/vegan/docs/vegdist
    # this formula is equivalent to formula with differnece
    
    s=0
    s1=0.0
    s2=0.0
    
    shared=set(dict1.keys()) & set(dict2.keys())
    for i in shared:
        s+=min(dict1[i],dict2[i])
    return 1-2.0*s/(sum(dict1.values())+sum(dict2.values()))




def Jaccard(dict1,dict2):
    ## Given a hash { 'species': count } for each sample, returns Compositional dissimilarity
    #Jaccard({'a': 0.2, 'b': 0.5, 'c': 0.3,},{'a': .4, 'b': 0.1, 'c':0.5 ,})
    #0.571428571429
    # confirmed by vegan R package
    #http://www.inside-r.org/packages/cran/vegan/docs/vegdist
    s=0
    s1=0.0
    s2=0.0
    
    shared=set(dict1.keys()) & set(dict2.keys())
    for i in shared:
        s+=min(dict1[i],dict2[i])
    b=1-2.0*s/(sum(dict1.values())+sum(dict2.values()))
    return 2*b/(1+b)


#===================================================================================
#===================================================================================



if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)

iProfile={}

from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(args.dir) if isfile(join(args.dir, f)) and f.endswith(args.extension)]




for sample in onlyfiles:
    dict={}
    
    with open(args.dir+"/"+sample,'r') as f:
        print args.dir+"/"+sample
        reader=csv.reader(f)
        print args.header
        if args.header=="1":
            headers = reader.next()
            print "rewerwerwer"
        for line in reader:
            print line
            column1=int(args.column1)
            column2=int(args.column2)

            dict[line[column1]]=float(line[column2])


    sample=sample.split(".")[0]
    
    iProfile[sample]=dict


out=open(args.outDir+'/beta_'+args.extension+'.csv',"w")
out.write("sample1,sample2,itemsShared,Soerenson,BrayCurtis,Jaccard"+"\n")


for sample1,dict1 in iProfile.items():
    for sample2,dict2 in iProfile.items():
        if sample1<sample2 and sample1!="" and sample2!="":
            if (len(dict1)!=0 and len(dict2)!=0):
                shared=set(dict1.keys()) & set(dict2.keys())
                #soerenson=1-len(shared)*2.0/(len(set(dict1.keys()))+len(set(dict2.keys())))
                out.write(sample1+","+sample2+","+str(len(shared))+","+str(Soerenson(dict1,dict2))+","+str(BrayCurtis(dict1,dict2))+","+str(Jaccard(dict1,dict2))+"\n")


#gama diversity
if not os.path.exists(args.outDir+"/QC"):
    os.makedirs(args.outDir+"QC")

events=[]
for sample,dict in iProfile.items():
    events+=dict.keys()


from collections import Counter


print len(set(events))

out=open(args.outDir+'/QC/eventsPerSample_'+args.extension+'.csv',"w")
out.write("VDJ recombination,Number of RNA-Seq samples\n")

for key,val in Counter(events).items():
    out.write(str(key)+","+str(val)+"\n")




print "DONE!"

#beta=1-2J/A+B
#where J is the sum of the lesser values for the shared taxa, A and B are the sum of the total values for all taxa for each sample respectively

#Compositional dissimilarity


