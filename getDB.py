#!/u/home/s/serghei/collab/code/rop/tools/Miniconda-Install/YourApplicationFolder/bin/python
import sys
import os
import argparse
import re
from subprocess import call, Popen, PIPE
import subprocess

#md5 values
humandict={}
humandict['viral_vipr']='9dce447328dfbc3a62cc7dd5b052242f'
humandict['fungi']='9f2d304fd5c49981682b2bb7a900a30e'
humandict['BWAindex']='4f009e3732d9f513e7b19b58edc41c13'
humandict['protozoa']='23e12115a5e9d526553c901e772731f5'
humandict['metaphlan']='3c9b9d6414d86a0c3d5018aefa5aaec4'
humandict['repeats']='109a97423f505b73a7e55995b827e2fd'
humandict['ribosomal.DNA']='9663a0e1121a0b122c8e23f41c558083'
humandict['viral']='7ce95144827603a64dc5996aa0112cc0'


dict={}
dict["human"]=humandict


#codeDir: the directory which contains this file
codeDir=os.path.dirname(os.path.realpath(__file__))

class SubprocessError(Exception):
	def __init__(self, message=""):
		print ("ERROR: Subprocess crashed. " + message)

def download(dbName,dbType,dirDB):
    """change directory to dirDB, and download dbName to it. Check that the  md5sum matches. dbName must be found on the server as a .tar.gz, but this command will extract and then delete the tar."""
    tarName=dirDB+"/"+"%s" %(dbName)
    md5sum=dict[dbType][dbName]
    os.chdir(dirDB)
    print ("-----------------------------")
    print humandict[dbName]
    md5sum=dict[dbType][dbName]
    os.chdir(dirDB)
    print ("-----------------------------")

    
    
    print tarName
    cmd=codeDir+"/db.scripts/download."+dbName+".sh " + tarName+".tar.gz"
    
    print cmd
    
    
    if subprocess.Popen([cmd], shell=True).wait(): raise SubprocessError()
    
    print ("Checking md5sum of %s database." %(dbName))
    cmd="md5sum "+tarName+".tar.gz"
    
    tarFile=tarName+".tar.gz"
    
    downloaded = Popen(["md5sum",tarFile], stdout=PIPE)
    
    
    checksum_downloaded = downloaded.communicate()[0].split()[0]

    if checksum_downloaded != md5sum:
        print ("DOWNLOAD of %s failed. Checksum does not match. Please re-run the script." %(dbName))
        sys.exit(23)
    else:
        print ("DOWNLOAD of %s was successful. Checksum matches." %(dbName))
        
    untar = Popen(["tar", "-zxvf", tarFile])
    if untar.wait():
        print ("Problem extracting downloaded DB in TAR format.")
    removetar = Popen(["rm", tarFile])
    if removetar.wait():
        print ("Problem deleting downloaded DB in TAR format.")
    
    

    
    

    
    

def symlink(dirDB):
	os.chdir(codeDir)
	symlinks = Popen(["ln", "-s", dirDB, "./"])
	if symlinks.wait():
		print ("problem creating link.")
	
####################################################################

print ("*********************************************")
print ("ROP (version 1.0.8) is a computational protocol aimed to discover the source of all reads, originated from complex RNA molecules, recombinant T and B cell receptors and microbial communities. Written by Serghei Mangul (smangul@ucla.edu) and Harry Taegyun Yang (harry2416@gmail.com), University of California, Los Angeles (UCLA). (c) 2016. Released under the terms of the General Public License version 3.0 (GPLv3)")
print ("")
print ("For more details see:")
print ("https://sergheimangul.wordpress.com/rop/")
print ("*********************************************")

#######################################################################
### Arguments
#######################################################################





ap = argparse.ArgumentParser('python getDB.py')


necessary_arguments = ap.add_argument_group('Necessary Inputs')
necessary_arguments.add_argument('dirDB', help='directory where the reference databases will be downloaded')


input_option_arguments = ap.add_argument_group('Select Database (multi-selection is possible)')
input_option_arguments.add_argument("--organism", help="downloads the (mouse) version of the database files", default="human")
input_option_arguments.add_argument("--repeat", help="Set up database for repeat sequences ONLY (lost repeat reads) ", action="store_true")
input_option_arguments.add_argument("--immune", help="Set up database for VDJ gene segments of B cell and T cell receptors ONLY (immune reads)", action="store_true")
input_option_arguments.add_argument("--metaphlan", help=" Set up database for Metaphlan2 ONLY", action="store_true")
input_option_arguments.add_argument("--microbiome", help="Set up database for microbiome ONLY", action="store_true")
input_option_arguments.add_argument("--viral", help="Set up database for viruses ONLY", action="store_true")
input_option_arguments.add_argument("--fungi", help="Set up database for fungi ONLY", action="store_true")
input_option_arguments.add_argument("--protozoa", help="Set up database for protozoa ONLY", action="store_true")


release = ap.add_argument_group('Connect database with the new release')
release.add_argument("--link2db", help="Connect the existing reference database with ROP", action="store_true")
release.add_argument("--f", help="Reconnects the ROP to the new database. Please note the existing link will removed", action="store_true")



args = ap.parse_args()




#=====================================================================================

#relative path to absolute path
args.dirDB=os.path.abspath(args.dirDB)


targeted=False

# IF none of them are selected: make everything true
if not args.repeat and not args.immune and not args.microbiome and not args.metaphlan and not args.viral and not args.fungi and not args.protozoa :
    args.repeat = True
    args.fungi = True
    args.microbiome = True
    args.metaphlan = True
    args.viral = True
    args.fungi = True
    args.protozoa = True
else:
	#It is gonna be non-reductive for now
	targeted = True



#set organism parameters
dbType=args.organism
dirDB=args.dirDB+"/db_"+dbType
currentDB=codeDir+"/db_"+dbType



if os.path.exists(currentDB):
	if args.f:
		removecurrDB = Popen(["rm", "-f", currentDB])
		if removecurrDB.wait():
			print ("Problem removing currentDB.")
	else:
		print ("Error : ROP was already linked to the reference database. To remove the current link please use --f option")
		sys.exit(23)

		
if args.link2db:
	chdir2 = Popen(["cd", codeDir])
	if chdir2.wait():
		print ("Problem changing directories (the second time).")
	if os.path.exists(dirDB):
		dblink = Popen(["ln", "-s", dirDB, "/db_"+dbType])
		if dblink.wait():
			print ("Error linking database to ROP.")
		print ("Reference databases are connected. Please use rop.py")
	else:
		print ("Error : Directory with reference databases doesn't exist! Please download the reference databases first. ")
else:
	if not os.path.exists(dirDB):
		os.makedirs(dirDB)
	else:
		print ("Error: Directory %s already exists. Please use  --link2db to connect the ROP with the reference databases." %(dirDB))
		sys.exit(23)






if not targeted:
	print ("The complete reference database will be downloaded.")
	print ("Downloading the database (~17GB) takes up to 15 minutes.")
	print ("Please wait until the downloading is completed.")


	for dbName in ['metaphlan','ribosomal.DNA','repeats','viral_vipr','viral','fungi','protozoa','BWAindex']:
		print ("Downloading %s database files ..." %(dbName))
		download(dbName,dbType,dirDB)

	symlink(dirDB)
	print ("Reference databases are ready. Please use rop.py")

else:
    

    if args.repeat:
        print "repeat"
        for dbName in ['ribosomal.DNA','BWAindex','repeats']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
        symlink(dirDB)
        print ("Reference databases are ready. Please use rop.py")

    if args.immune:
		for dbName in ['ribosomal.DNA','BWAindex']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")
		
    if args.metaphlan:
		for dbName in ['metaphlan','ribosomal.DNA','BWAindex']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")

    if args.microbiome:
		print ("--microbiome options was selected. Refrence database for microbiome will be downloaded")
		for dbName in ['ribosomal.DNA','BWAindex','viral','viral_vipr','fungi','protozoa','metaphlan']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")
    if args.viral:
		print ("--viral options was selected. Refrence database for microbiome will be downloaded")
		for dbName in ['ribosomal.DNA','BWAindex','viral','viral_vipr']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")

    if args.fungi:
		print ("--fungi options was selected. Refrence database for microbiome will be downloaded")
		for dbName in ['ribosomal.DNA','BWAindex','fungi']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")
    if args.protozoa:
		print ("--protozoa options was selected. Refrence database for microbiome will be downloaded")
		for dbName in ['ribosomal.DNA','BWAindex','protozoa']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")
