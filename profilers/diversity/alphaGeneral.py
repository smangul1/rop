import sys
import csv
import os
import argparse
import collections
from math import log as ln


def p(n, N):
    """ Relative abundance """
    if n is  0:
        return 0
    else:
        return (float(n)/N) * ln(float(n)/N)

def sdi(data):
    N = sum(data.values())
    
    return -sum(p(n, N) for n in data.values() if n is not 0)


#Looks like it doesn't have relative abundance
#def InverseSimpson(data):
#    return 1/sum(float(n*n)for n in data.values() if n is not 0)


#===================================================================================
#===================================================================================

ap = argparse.ArgumentParser()
ap.add_argument('dir', help='dir')
ap.add_argument('outDir', help='dir to save the results')
ap.add_argument('prefix', help='for example igh,tcra,circular')
ap.add_argument('column1', help='number of column1')
ap.add_argument('column2', help='number of column2')
ap.add_argument('fileExtension', help='suffix - extension of th efiles in directory dir')



args = ap.parse_args()



if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)


iProfile={}


from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(args.dir) if isfile(join(args.dir, f))]

c1=int(args.column1)
c2=int(args.column2)


print "Number of samples to proccess:", len(onlyfiles)





for sample in onlyfiles:
    if sample.endswith(args.fileExtension):
    
        dict={}
        
        with open(args.dir+"/"+sample,'r') as f:
            reader=csv.reader(f)
            headers = reader.next()
            for line in reader:
                dict[line[c1]]=float(line[c2])
        ##USE this one
        #sample=sample.split(args.fileExtension)[0]
        #TEMP for GTEX
        sample=sample.split(".")[0]+"."+sample.split(".")[1]+"."+sample.split(".")[2]
        iProfile[sample]=dict

out=open(args.outDir+'/alpha'+args.prefix+'.csv',"w")
out.write("File_Name,load__%s,richness__%s,Shannon__%s" %(args.prefix,args.prefix,args.prefix)+"\n")






for sample,dict in iProfile.items():
    sdiTemp=0.0
    if len(dict)==1:
        sdiTemp=0.0
    else:
        sdiTemp=sdi(dict)


    if len(dict)!=0:
        out.write(sample+","+str(sum(dict.values()))+","+str(len(dict))+","+str(sdiTemp)+"\n")
    else:
        out.write(sample+",0,0,0.0\n")


print "Calculated alpha for %i samples" %len(iProfile)


print "DONE!"
