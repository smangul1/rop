import sys
import os 
import argparse
from subprocess import call	


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
	call(["tar","-zxvf", "./db/immune.tar"])
elif args.standard: 
	print "Standard installation option selected"
	print "Downloading the gtex files"
	call(["wget", ""]) # FIX IT - ADD THE HTTPS ADDRESS
	print "Unzipping the databases"
	call(["tar", "-zxvf", '../db'])
else:
	print "wrong"



