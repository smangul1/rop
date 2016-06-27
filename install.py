import sys
import os 
import argparse
from subprocess import call, Popen, PIPE



ap = argparse.ArgumentParser('python install.py')

necessary_arguments = ap.add_argument_group('Necessary Inputs')
necessary_arguments.add_argument("--standard", help="standard installation option", action="store_true")\


input_option_arguments = ap.add_argument_group('Select Database')
input_option_arguments.add_argument("--repeat", help="Set up database for repeat sequences ONLY (lost repeat reads) ", action="store_true")
input_option_arguments.add_argument("--immune", help="Set up database for VDJ gene segments of B cell and T cell receptors (immune reads)", action="store_true")

release = ap.add_argument_group('Connect database with the new release')
release.add_argument("--link2db", help="Set up database for VDJ gene segments of B cell and T cell receptors (immune reads)", action="store_true")



args = ap.parse_args()

#=====================================================================================

if args.standard:
	checksum_original="1052af7849f099ed77fa6e668278a4ec" 
	print "Standard installation option selected"
	print "Downloading the database (~65GB) takes up to 45 minute"
	print "Please wait until the installation is completed."
	print "Downloading the database files"
	os.chdir('./db/')
	call(["wget", "https://googledrive.com/host/0Bx1fyWeQo3cOMjFNMzBrcWZfXzA/rRNA.tar.gz", "--no-check-certificate"])
	
	print "Checking md5sum of the databases" 
	downloaded = Popen(["md5sum",'./database.tar'], stdout=PIPE)
	checksum_downloaded = downloaded.communicate()[0].split()[0]
    
    
    sys.exit(1)

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


