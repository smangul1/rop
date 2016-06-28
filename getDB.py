import sys
import os 
import argparse
from subprocess import call, Popen, PIPE


#codeDir
codeDir=os.path.dirname(os.path.realpath(__file__))

#Gets MD5 from file
def getmd5(filename):
    return m.hexdigest()


####################################################################

print "*********************************************"
print "ROP (version 1.0.3) is a computational protocol aimed to discover the source of all reads, originated from complex RNA molecules, recombinant antibodies and microbial communities. Written by Serghei Mangul (smangul@ucla.edu) and Harry Taegyun Yang (harry2416@gmail.com), University of California, Los Angeles (UCLA). (c) 2016. Released under the terms of the General Public License version 3.0 (GPLv3)"
print ""
print "For more details see:"
print "https://sergheimangul.wordpress.com/rop/"
print "*********************************************"

#######################################################################
### Arguments
#######################################################################


ap = argparse.ArgumentParser('python rop.py')


necessary_arguments = ap.add_argument_group('Necessary Inputs')
necessary_arguments.add_argument('dirDB', help='directory where the reference databases will be downloaded')


input_option_arguments = ap.add_argument_group('Select Database')
input_option_arguments.add_argument("--repeat", help="Set up database for repeat sequences ONLY (lost repeat reads) ", action="store_true")
input_option_arguments.add_argument("--immune", help="Set up database for VDJ gene segments of B cell and T cell receptors (immune reads)", action="store_true")

release = ap.add_argument_group('Connect database with the new release')
release.add_argument("--link2db", help="Connect the reference database with ROP", action="store_true")
release.add_argument("--f", help="Reconnect the ROP to the new database provided. Please note the existing link will romoved", action="store_true")



args = ap.parse_args()

#=====================================================================================

dirDB=args.dirDB+"/db/"

if not os.path.exists(dirDB):
    os.makedirs(dirDB)
else:
    print "Error: Directory %s already exists. Please use  --link2db to connect the ROP with the reference databases " %(dirDB)
    sys.exit(23)


currentDB=codeDir+'/db'



if os.path.exists(currentDB):
    if args.f:
        cmd='rm -f %s' %(currentDB)
        os.system(cmd)
    else:
        print "Error : ROP was already linked to the reference database. To remove the current link please use --f option"
        sys.exit(23)







if not args.repeat and not args.immune:
    print "Standard installation option selected"
    print "Downloading the database (~65GB) takes up to 45 minute"
    print "Please wait until the installation is completed."
    print "Downloading the database files"


    for dbName in ['rRNA','bowtie2Index','repeats']:
    
        print "Downloading "
        os.chdir(dirDB)
        tarName=dirDB+"/"+"%s.tar.gz" %(dbName)
        
        
        link="https://googledrive.com/host/0Bx1fyWeQo3cOMjFNMzBrcWZfXzA/%s.tar.gz" %(dbName)
        print "-----------------------------"
        print "Downloading %s" %(link)

        call(["wget", link, "--no-check-certificate"])
        



        cmd="tar -zxvf %s" %(tarName)
        
        
        
        os.system(cmd)
        os.remove(tarName)
        
    os.chdir(codeDir)
    cmd="ln -s %s ./" %(dirDB)
    os.system(cmd)
    print "Installation Completed! Please use rop.py"

else:
    print "No option is selected. Please use -h option to see available options"
    sys.exit(233)


