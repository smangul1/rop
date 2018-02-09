import sys
import os
import argparse
import re
from subprocess import call, Popen, PIPE

#md5 values
humandict={}
humandict['antibody']='37dd2b878d69734fa96fe6f4e182b703'
humandict['bacteria']='dae3d5d9220b1c40bac72cbf9f710d5c'
humandict['BWAIndex']='f9f61c97917adf37a738b31d9888a1a6'
humandict['eupathdb']='e38ac558fdcbe4ae85745016712b0e8b' 
humandict['metaphlan_db']='83e0be082014811189bcb57856d3213d'  
humandict['repeats']='132b9419746ea542eab966c06b090587'  
humandict['rRNA']='93756de9ed7a6f0488cdebf7e5a3eb93'
humandict['virus']='02a7c899924ce6bac64472abeade2e6d'  

mousedict={}
mousedict['antibody']='5fcbd7d0bfcf534317b2643b1aa7ee3d'
mousedict['bacteria']='5dd33ee0456bbbca077f13338d1fdbbd'
mousedict['bowtie2Index']='e8437ca0f382ca41c86f3b8c8c290999'
mousedict['BWAIndex']='915b3673f970f3045a2410860589839f'
mousedict['eupathdb']='ea35b8e758df4939540c8540bed4e7c4'
mousedict['metaphlan_db']='2ecd61ee68065903011446e37fecc186'
mousedict['repeats']='fc4c643475ce8b6c74136a6b2018b445'
mousedict['rRNA']='c6f70df0fcd70ae14e1963301644e134'
mousedict['virus']='02a7c899924ce6bac64472abeade2e6d'

dict={}
dict["mouse"]=mousedict
dict["human"]=humandict

humanlinks={}
#humanlinks['bacteria']=["https://drive.google.com/uc?export=download","&id=0Bx1fyWeQo3cOdWt6NzU4TGZONjA"]
humanlinks['BWAIndex']=["https://drive.google.com/uc?export=download","&id=1g1An0PEAkFisZojrygxovwKCbmtC4CVi"]
humanlinks['metaphlan_db']=["https://drive.google.com/uc?export=download","&id=0Bx1fyWeQo3cOVUVkcFNjdU5BbXM"]
humanlinks['repeats']=["https://drive.google.com/uc?export=download","&id=1C03gbi7xXiBQ2_0cIYarlggXXlq7QEO4"]
humanlinks['rRNA']=["https://drive.google.com/uc?export=download","&id=1_lkqJXh98jsk5wpqhPMiPnjeY1BizDP1"]
humanlinks['viral']=["https://drive.google.com/uc?export=download","&id=152FCuyWkY2oDhcqyIGImu4aE4lIu1AxN"]
humanlinks['fungi']=["https://drive.google.com/uc?export=download","&id=1lcR_2E-QWThFe0rdW80iPo96JVg4I3-F"]
humanlinks['protozoa']=["https://drive.google.com/uc?export=download","&id=1KUUSyDBnGjl4oXUf9ZmLiXCVbyHB06lP"]
humanlinks['viral_vipr']=["https://drive.google.com/uc?export=download","&id=1DCFMU2paZwn6rA3SpGlW5G-v52_sqr1B"]

mouselinks={}

mouselinks['antibody']=["https://drive.google.com/uc?export=download","&id=0Bx1fyWeQo3cOSzJxWm5obXBtc3c"]
mouselinks['BWAIndex']=["https://drive.google.com/uc?export=download","&id=0Bx1fyWeQo3cOR083NnNackYzdWM"]
mouselinks['eupathdb']=["https://drive.google.com/uc?export=download","&id=0Bx1fyWeQo3cOUHJHWUZFYzZoM1k"]
mouselinks['metaphlan_db']=["https://drive.google.com/uc?export=download","&id=0Bx1fyWeQo3cOeFA2bi02RlRCMGc"]
mouselinks['repeats']=["https://drive.google.com/uc?export=download","&id=0Bx1fyWeQo3cOWVFPZGhDd0g2SUk"]
mouselinks['rRNA']=["https://drive.google.com/uc?export=download","&id=0Bx1fyWeQo3cOR0kzeWVtS1Q4MkU"]
mouselinks['virus']=["https://drive.google.com/uc?export=download","&id=0Bx1fyWeQo3cObXVfMDd2aENNbms"]

links={}	
links["human"]=humanlinks
links["mouse"]=mouselinks

#codeDir: the directory which contains this file
codeDir=os.path.dirname(os.path.realpath(__file__))

def download(dbName,dbType,dirDB):
	"""change directory to dirDB, and download dbName to it. Check that the 
	md5sum matches. dbName must be found on the server as a .tar.gz, but this
	command will extract and then delete the tar."""
	
	tarName=dirDB+"/"+"%s" %(dbName)
	link=links[dbType][dbName][0] + links[dbType][dbName][1]

	link="https://drive.google.com/open?id=1_lkqJXh98jsk5wpqhPMiPnjeY1BizDP1"

	md5sum=dict[dbType][dbName]
	os.chdir(dirDB)
	print ("-----------------------------")
	print ("Attempting to download %s." %(link))
	if dbName in ['rRNA', 'antibody', 'repeats']:
		getting = Popen(["wget", link, "--output-document="+tarName, "--no-check-certificate"])
		if getting.wait():
			print ("problem downloading online database.")

	else:
		htmlName=dirDB+"/confirm"+dbName+".html"
		getting = Popen(["wget", link, "--output-document="+htmlName, "-nv", "--keep-session-cookies", "--save-cookies", "./cookies.txt", "--no-check-certificate"])
		if getting.wait():
			print ("problem downloading gdrive html file.")
		with open(htmlName, "r") as cfile:
			code = cfile.read()
			newlink = links[dbType][dbName][0] + "&confirm=" + re.search("confirm=(.{4})", code).group(1) + links[dbType][dbName][1]
			getting = Popen(["wget", newlink, "--output-document="+tarName, "-nv", "--load-cookies", "./cookies.txt", "--no-check-certificate"])
			if getting.wait():
				print ("problem downloading online database.")
		removehtml = Popen(["rm", "-f", htmlName])
		if removehtml.wait():
			print ("problem deleting gdrive html page.")
	
	print ("Checking md5sum of %s database." %(dbName))
	downloaded = Popen(["md5sum",tarName], stdout=PIPE)
	checksum_downloaded = downloaded.communicate()[0].split()[0]
	print (checksum_downloaded)
	if checksum_downloaded != md5sum:
		print ("DOWNLOAD of %s failed. Checksum does not match. Please re-run the script." %(link))
		sys.exit(23)

	untar = Popen(["tar", "-zxvf", tarName])
	if untar.wait():
		print ("problem extracting downloaded DB in TAR format.")
	removetar = Popen(["rm", tarName])
	if removetar.wait():
		print ("problem deleting downloaded DB in TAR format.")
	
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
if not args.repeat and not args.immune and not args.microbiome and not args.metaphlan and not args.viral and not args.fungi and not not args.protozoa :
	args.repeat = True
	args.immune = True
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
			print ("problem removing currentDB.")
	else:
		print ("Error : ROP was already linked to the reference database. To remove the current link please use --f option")
		sys.exit(23)

		
if args.link2db:
	chdir2 = Popen(["cd", codeDir])
	if chdir2.wait():
		print ("problem changing directories (the second time).")
	if os.path.exists(dirDB):
		dblink = Popen(["ln", "-s", dirDB, "/db_"+dbType])
		if dblink.wait():
			print ("error linking database to ROP.")
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
	print ("Downloading the database (~65GB) takes up to 45 minutes.")
	print ("Please wait until the installation is completed.")


	for dbName in ['rRNA','repeats','BWAIndex','metaphlan_db','antibody','viral','fungi','protozoa']:
		print ("Downloading %s database files ..." %(dbName))
		download(dbName,dbType,dirDB)

	symlink(dirDB)
	print ("Reference databases are ready. Please use rop.py")

else:


	if args.repeat:
		for dbName in ['rRNA','BWAIndex','repeats']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")

	if args.immune:
		for dbName in ['rRNA','BWAIndex']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")
		
	if args.metaphlan:
		for dbName in ['rRNA','BWAIndex','metaphlan_db']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")

	if args.microbiome:
		print ("--microbiome options was selected. Refrence database for microbiome will be downloaded")
		for dbName in ['rRNA','BWAIndex','viral','fungi','protozoa','metahplan']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")
	if args.viral:
		print ("--viral options was selected. Refrence database for microbiome will be downloaded")
		for dbName in ['rRNA','BWAIndex','viral']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")

	if args.fungi:
		print ("--fungi options was selected. Refrence database for microbiome will be downloaded")
		for dbName in ['rRNA','BWAIndex','fungi']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")
	if args.protozoa:
		print ("--protozoa options was selected. Refrence database for microbiome will be downloaded")
		for dbName in ['rRNA','BWAIndex','protozoa']:
			print ("Downloading %s database files ..." %(dbName))
			download(dbName,dbType,dirDB)
		symlink(dirDB)
		print ("Reference databases are ready. Please use rop.py")