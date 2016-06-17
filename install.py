import sys
import os 
import argparse
from subprocess import call, Popen, PIPE


"""
option requirements - TCR/BCR only option - in terms of database
"""

"""
make environment option 
env path name = ROPDB

"""

ap = argparse.ArgumentParser()
ap.add_argument("rop_db", help = "Directory ")
ap.add_argument("--immune", help="Set up database for immune reads only (i.e. TCR/BCR) ", action="store_true")
ap.add_argument("--standard", help="standard installation option", action="store_true")
args = ap.parse_args()

try:
	rop_db_path = os.environ['ROPDB']
	print "ROP DB is found: new links will be generated."
	# TODO : MAKE THE LINK IF DB IS FOUND 
except KeyError:
	print "Database for ROP will be unzipped at : %s " %(args.rop_db)
	# TODO : install DB and make link

if args.immune:
	print "Immune only option selected"
	print "Unzipping Immune Databases (i.e. BCR/TCR)"
	os.chdir('./db/')
	call(["tar","-xvf", './immune.tar'])

elif args.standard: 
	checksum_original="1052af7849f099ed77fa6e668278a4ec" 
	print "Standard installation option selected"
	print "Downloading the database (~65GB) takes up to 45 minute"
	print "Please wait until the installation is completed."
	print "Downloading the database files"
	os.chdir('./db/')
	call(["wget", "https://googledrive.com/host/0B_NUyiE86yDwaUxoVjhlSjN5SkE/database.tar", "--no-check-certificate"]) 
	
	print "Checking md5sum of the databases" 
	downloaded = Popen(["md5sum",'./database.tar'], stdout=PIPE)
	checksum_downloaded = downloaded.communicate()[0].split()[0]

	if checksum_downloaded != checksum_original:
		print "DOWNLOAD failed. Please re-run the script"
		sys.exit(23)
	else:
		print "MD5 Checksum matches"
	print "Downloading the metaphlan database"
	call(["wget", "https://googledrive.com/host/0B_NUyiE86yDwaUxoVjhlSjN5SkE/metaphlan_db.tar", "--no-check-certificate"]) 

	print "Unzipping the databases"
	call(["tar","-xvf", 'database.tar'])
	os.remove('database.tar')
	call(["tar", "-xvf", 'metaphlan_db.tar'])
	os.remove('metaphlan_db.tar')
	print "Installation Completed! Please use rop.py"

else:
	print "No option is selected. Please use -h option to see available options"
	sys.exit(233)


