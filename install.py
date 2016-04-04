import sys
import os 
import argparse
from subprocess import call, Popen, PIPE


"""
option requirements - TCR/BCR only option - in terms of database


"""

ap = argparse.ArgumentParser()
ap.add_argument("--immune", help="Set up database for immune reads only (i.e. TCR/BCR", action="store_true")
ap.add_argument("--standard", help="standard installation option", action="store_true")
args = ap.parse_args()

if args.immune:
	print "Immune only option selected"
	print "Unzipping Immune Databases (i.e. BCR/TCR)"
	call(["tar","-zxvf", './db/immune.tar'])
elif args.standard: 
	checksum_original="" # update checksum 
	print "Standard installation option selected"
	
	print "Downloading the gtex files"
	os.chdir('./db/')
	call(["wget", ""]) # FIX IT - ADD THE HTTPS ADDRESS
	
	print "Checking md5sum of the databases" 
	downlaoded = Popen(["md5sum",'./db/database.tar'],stdout=PIPE)
	checksum_downloaded = checsum_test.communicate()[0].split()[0]

	if checksum_downloaded != checksum_original:
		return "DOWNLOAD failed. Please re-run the script"
		sys.exit(23)
	else:
		return "MD5 Checksum matches"

	print "Unzipping the databases"
	call(["tar","-zxvf", 'database.tar'])


	print "Unzipping the databases"
	call(["tar", "-zxvf", './db/'])
else:
	print "No option is selected. Please use -h option to see available options"



