import sys
import os 
import argparse
from subprocess import call, Popen, PIPE


#codeDir
codeDir=os.path.dirname(os.path.realpath(__file__))

def download(dbName,md5sum,dirDB):

    os.chdir(dirDB)
    tarName=dirDB+"/"+"%s.tar.gz" %(dbName)
    link="https://googledrive.com/host/0Bx1fyWeQo3cOMjFNMzBrcWZfXzA/%s.tar.gz" %(dbName)
    print "-----------------------------"
    print "Downloading %s" %(link)
    call(["wget", link, "--no-check-certificate"])
    print "Checking md5sum of %s database" %(dbName)
    downloaded = Popen(["md5sum",tarName], stdout=PIPE)
    checksum_downloaded = downloaded.communicate()[0].split()[0]
    print checksum_downloaded
        
    if checksum_downloaded != md5sum:
            print "DOWNLOAD of %s failed. Please re-run the script" %(link)
            sys.exit(23)
        
        
    cmd="tar -zxvf %s" %(tarName)
    os.system(cmd)
    os.remove(tarName)


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


input_option_arguments = ap.add_argument_group('Select Database (multi-selection is possible)')
input_option_arguments.add_argument("--repeat", help="Set up database for repeat sequences ONLY (lost repeat reads) ", action="store_true")
input_option_arguments.add_argument("--immune", help="Set up database for VDJ gene segments of B cell and T cell receptors ONLY (immune reads)", action="store_true")
input_option_arguments.add_argument("--metaphlan", help = " Set up database for Metaphlan2 ONLY", action = "store_true")
input_option_arguments.add_argument("--circRNA", help = "Set up database for circular RNAs ONLY", action="store_true")
input_option_arguments.add_argument("--microbiome", help = "Set up database for microbime  ONLY", action = "store_true")


release = ap.add_argument_group('Connect database with the new release')
release.add_argument("--link2db", help="Connect the reference database with ROP", action="store_true")
release.add_argument("--f", help="Reconnect the ROP to the new database provided. Please note the existing link will removed", action="store_true")



args = ap.parse_args()

#=====================================================================================

targeted=False

# IF none of them are selected: make everything true
if not args.repeat and not args.immune and not args.circRNA and not args.microbiome and not args.metaphlan:
    args.repeat = True
    args.immune = True
    args.circRNA = True
    args.microbiome = True
    args.metaphlan = True
else:
    #It is gonna be non-reductive for now
    targeted = True




dirDB=args.dirDB+"/db/"

currentDB=codeDir+'/db'





if os.path.exists(currentDB):
    if args.f:
        cmd='rm -f %s' %(currentDB)
        os.system(cmd)
    else:
        print "Error : ROP was already linked to the reference database. To remove the current link please use --f option"
        sys.exit(23)


dict={}
dict['rRNA']='3f85a14df97ad844e5175265e54148f5'
dict['bowtie2Index']='5988f3f86b41252e005f78b6a022b7b8'
dict['repeats']='34948f44b0a393228e93ce6bd2c4d512'
dict['BWAIndex']='a48d8b93a69ac5a4301ede641e51e429'
dict['metaphlan_db']='c0e06113269d1042c84b32abef99951c'
dict['antibody']='8f2d8884677fe2e05736a3b92a57739b'
dict['virus']='82980243db39b2d96bd2ffc633e40ce2'
dict['eupathdb']='34c6ffb1f0d8664ba549ce45fefa8724'
dict['bacteria']='19f09bfe965693ff9778e177c57a21bb'


if args.link2db:
    os.chdir(codeDir)
    
    if os.path.exists(args.dirDB):
        cmd="ln -s %s db" %(args.dirDB)
        os.system(cmd)
        print "Reference databases are connected. Please use rop.py"
    else:
        print "Error : Directory with reference databases doesn't exist! Please download the reference databases first. "
else:
    if not os.path.exists(dirDB):
        os.makedirs(dirDB)
    else:
        print "Error: Directory %s already exists. Please use  --link2db to connect the ROP with the reference databases " %(dirDB)
        sys.exit(23)





if not targeted:
    
    print "The complete refrence database will be downloaded"
    




    
    
    print "Downloading the database (~65GB) takes up to 45 minute"
    print "Please wait until the installation is completed."





    for dbName in ['rRNA','bowtie2Index','repeats','BWAIndex','metaphlan_db','antibody','virus','eupathdb','bacteria']:
        print "Downloading %s database files ..." %(dbName)
        download(dbName,dict[dbName],dirDB)

        
    os.chdir(codeDir)
    cmd="ln -s %s ./" %(dirDB)
    os.system(cmd)
    print "Reference databases are ready. Please use rop.py"

else:
    
    
    if args.repeat:
        for dbName in ['rRNA','bowtie2Index','repeats']:
            print "Downloading %s database files ..." %(dbName)
            download(dbName,dict[dbName],dirDB)
        os.chdir(codeDir)
        cmd="ln -s %s ./" %(dirDB)
        os.system(cmd)
        print "Reference databases are ready. Please use rop.py"

    if args.immune:
        for dbName in ['rRNA','bowtie2Index','antibody']:
            print "Downloading %s database files ..." %(dbName)
            download(dbName,dict[dbName],dirDB)
        os.chdir(codeDir)
        cmd="ln -s %s ./" %(dirDB)
        os.system(cmd)
        print "Reference databases are ready. Please use rop.py"
    if args.metaphlan:
        for dbName in ['rRNA','bowtie2Index','metaphlan_db']:
            print "Downloading %s database files ..." %(dbName)
            download(dbName,dict[dbName],dirDB)
        os.chdir(codeDir)
        cmd="ln -s %s ./" %(dirDB)
        os.system(cmd)
        print "Reference databases are ready. Please use rop.py"
    if args.circRNA:
        for dbName in ['rRNA','bowtie2Index','BWAIndex']:
            print "Downloading %s database files ..." %(dbName)
            download(dbName,dict[dbName],dirDB)
        os.chdir(codeDir)
        cmd="ln -s %s ./" %(dirDB)
        os.system(cmd)
        print "Reference databases are ready. Please use rop.py"
    if args.microbiome:
        print "--mirobiome options was selected. Refrence database for microbiome will be downloaded"
        for dbName in ['bacteria','rRNA','bowtie2Index','virus','eupathdb']:
            print "Downloading %s database files ..." %(dbName)
            download(dbName,dict[dbName],dirDB)
        os.chdir(codeDir)
        cmd="ln -s %s ./" %(dirDB)
        os.system(cmd)
        print "Reference databases are ready. Please use rop.py"


