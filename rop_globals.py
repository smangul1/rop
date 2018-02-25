"""********************************************************************************
ROP is a computational protocol aimed to discover the source of all reads, 
originated from complex RNA molecules, recombinant antibodies and microbial 
communities. 

Written by Serghei Mangul (smangul@ucla.edu), Harry Taegyun Yang 
(harry2416@gmail.com), Kevin Hsieh (kevin.hsieh@ucla.edu), and Linus Chen 
(u6.30cl@gmail.com). University of California, Los Angeles (UCLA). (c) 2016. 

Released under the terms of the General Public License version 3.0 (GPLv3)

For more details see: https://sergheimangul.wordpress.com/rop/
ROP Tutorial: https://github.com/smangul1/rop/wiki 
********************************************************************************"""

################################################################################
### IMPORTS
################################################################################

import sys, os  # system
import argparse, subprocess  # utilities
import csv, gzip  # file I/O
import pysam


CD = os.path.dirname(os.path.abspath(__file__))
sys.path.append(CD + "/tools/biopython/biopython-1.66/")
from Bio import SeqIO  # module needed for sequence input
# for future use:
# sys.path.append(CD + "/tools/pysam-master/")


################################################################################
### ARGUMENTS 
################################################################################

ap = argparse.ArgumentParser("python rop.py")

necessary_arguments = ap.add_argument_group("Necessary Inputs")
necessary_arguments.add_argument("unmappedReads", 
	help="unmapped reads in .fastq format unless options requiring a " +\
	  "format are selected")
necessary_arguments.add_argument("dir", 
	help="directory to save results of the analysis")

job_option_arguments = ap.add_argument_group('Job Options')
job_option_arguments.add_argument("--qsub", 
	help="submit qsub jobs on Hoffman2 (UCLA) cluster. If planning to use " +\
	  "on your cluster, contact smangul@ucla.edu", 
	action="store_true")
job_option_arguments.add_argument("--qsubArray", 
	help="prepare qsub scripts to be run later using job array. Working on " +\
	  "Hoffman2 (UCLA) cluster. If planning to use on your cluster, contact " +\
	  "smangul@ucla.edu", 
	action="store_true")
job_option_arguments.add_argument("--maui", 
	help="use this option together with --qsub to submit jobs via Maui " +\
	  "scheduler. Maui is a job scheduler developped by Adaptive Computing. " +\
	  "More details are here: https://wiki.calculquebec.ca/w/Maui/en", 
	action="store_true")
job_option_arguments.add_argument("--organism", 
	help="run ROP for a specified organism (mouse) instead of human", 
	action="store", default="human")

input_option_arguments = ap.add_argument_group("Input Options")
input_option_arguments.add_argument("--b", "-b", 
	help="unmapped reads in .bam format. Do not use with --gzip" +\
	  "option in this group", 
	action="store_true")
input_option_arguments.add_argument("--gzip", "-z", 
	help="unmapped reads in .gz format. Do not use with --b", 
	action="store_true")
input_option_arguments.add_argument("--skipLowq", 
	help="skip filtering low quality reads. The input reads need to be in " +\
	  ".fasta or .fasta.gz format", 
	action="store_true")
input_option_arguments.add_argument("--skipQC", 
	help="skip entire QC step (filtering low-quality, low-complexity and " +\
	  "rRNA reads). The input reads need to be in .fasta or .fasta.gz format", 
	action="store_true")
input_option_arguments.add_argument("--skipPreliminary", "-s", 
	help="skip the preliminary steps including (1) QC and (2) remapping to " +\
	  "references (lost reads). The input reads need to be in .fasta or " +\
	  ".fasta.gz format", 
	action="store_true")
input_option_arguments.add_argument("--max", "-m",help=" account for maximum number of unmapped reads. We use liberal threshold to account for maximum number of reads",action="store_true")
input_option_arguments.add_argument("--pe", "-pe",help="report number of discordant read pairs, with the reads from the same pair classified into tdifferent classes",action="store_true")


run_only_options = ap.add_argument_group("Run Options - select analysis " +\
  "(runs all if none selected)")
run_only_options.add_argument("--repeat", 
	help="Run lost repeat profiling", 
	action="store_true")
run_only_options.add_argument("--circRNA", 
	help="Run circular RNA profiling", 
	action="store_true")
run_only_options.add_argument("--immune", 
	help="Run BCR/TCR profiling", 
	action="store_true")
run_only_options.add_argument("--microbiome", 
	help="Run microbime profiling", 
	action="store_true")
run_only_options.add_argument("--bacteria",help="Run bacteria profiling",action="store_true")
run_only_options.add_argument("--viral",help="Run viral profiling",action="store_true")
run_only_options.add_argument("--fungi",help="Run fungi profiling",action="store_true")
run_only_options.add_argument("--protozoa",help="Run protozoa profiling",action="store_true")



misc_option_arguments = ap.add_argument_group('Miscellenous Options')
misc_option_arguments.add_argument("--outGz", 
	help="(FUTURE USE) intermediate .fasta files are stored as .fasta.gz", 
	action="store_true")
misc_option_arguments.add_argument("--rezip", 
	help="(FUTURE USE) rezip fasta files after analysis", 
	action="store_true")
misc_option_arguments.add_argument("--clean", 
	help="clean more intermediate files for maximum space efficiency",
	action="store_true")
misc_option_arguments.add_argument("--quiet", 
	help="suppress progress report and warnings", 
	action="store_true")
misc_option_arguments.add_argument("--dev", 
	help="keep intermediate files", 
	action="store_true")
misc_option_arguments.add_argument("--nonReductive", 
	help="(USE WITH CAUTION) non-reductive analysis", 
	action="store_true")
misc_option_arguments.add_argument("--f", "-f", 
	help="(USE WITH CAUTION) overwrite the analysis directory (provided as " +\
	  "dir option)", 
	action="store_true")

ARGS = ap.parse_args()

# if no analysis mode is selected, select everything
if (not ARGS.repeat and not ARGS.immune and not ARGS.circRNA and not ARGS.microbiome):
    ARGS.repeat = True
    ARGS.immune = True
    ARGS.circRNA = False # starting from release v1.0.8 circRNA is no longer dafault, as this is too slow
    ARGS.metaphlan = True
    ARGS.microbiome = True
    ARGS.bacteria = False # starting from release v1.0.8 mapping to bacteria is no longer dafault, as this is too slow
    ARGS.viral = True
    ARGS.fungi = True
    ARGS.protozoa = True

# necessary for parallelization
if ARGS.qsub or ARGS.qsubArray:
	ARGS.nonReductive = True
	
# relative path to absolute path
ARGS.unmappedReads = os.path.abspath(ARGS.unmappedReads)
ARGS.dir = os.path.abspath(ARGS.dir)

# check input file type
filename, file_extension = os.path.splitext(ARGS.unmappedReads)
if file_extension == ".bam" and not ARGS.b:
	print ("WARNING: Detected bam input, but --b option is not selected. Processing as .bam")
	ARGS.b = True
elif file_extension == ".gz" and not ARGS.gzip:
	print ("WARNING: Detected gzip input, but --gzip option is not selected. Processing as .gzip")
	ARGS.gzip = True

# check if ARGS.dir exists
if os.path.exists(ARGS.dir) and not ARGS.f:
	print ("ERROR: The directory " + ARGS.dir + " exists. Please choose a different directory in which to save results of the analysis. Alternatively, use the --f option to overwrite the results into " + ARGS.dir + ".")
	sys.exit(1)
if os.path.exists(ARGS.dir) and ARGS.f:
    cmd = "rm -fr " + ARGS.dir + " &>/dev/null"
    os.system(cmd)
    print (ARGS.dir,"was deleted")


# database folder path, BASENAME (name of unmapped reads file without extension)
DB_FOLDER = "/db_" + ARGS.organism
BASENAME = os.path.basename(ARGS.unmappedReads).split('.')[0]


################################################################################
### VARIABLES
################################################################################
# Analysis directories

DIRS = {}
DIRS["QC"] = ARGS.dir + "/QC/"  # step 1
DIRS["lostReads"] = ARGS.dir + "/lost.human.reads/"  # step 2
DIRS["lostRepeat"] = ARGS.dir + "/lost.repeats/"  # step 3
DIRS["antibody"] = ARGS.dir + "/immune.repertoire/"  # step 4
DIRS["NCL"] = ARGS.dir + "/NCL/"  # step 5
DIRS["microbiome"] = ARGS.dir + "/microbiome/"  # step 6
DIRS["metaphlan"] =  DIRS["microbiome"] + "metaphlan/"
DIRS["bacteria"] = DIRS["microbiome"] + "bacteria/"
DIRS["virus"] = DIRS["microbiome"] + "viral/"
DIRS["fungi"] = DIRS["microbiome"] + "fungi/"
DIRS["protozoa"] = DIRS["microbiome"] + "protozoa/"

for dir in DIRS:
	if not os.path.exists(DIRS[dir]):
		os.makedirs(DIRS[dir])

################################################################################
# Intermediate file names

INTFNS = {}
INTFNS["unmappedFastq"] = ARGS.dir + "/unmapped_" + BASENAME + ".fastq"  # begin
INTFNS["lowQFileFasta"] = DIRS["QC"] + BASENAME + "_lowQ.fa"  # after step 1a
INTFNS["lowQCFile"] = DIRS["QC"] + BASENAME + "_lowQC.fa"  # after step 1b
INTFNS["rRNAFile"] = DIRS["QC"] + BASENAME + "_rRNA_blastFormat6.csv"
INTFNS["afterrRNAFasta"] = DIRS["QC"] + BASENAME + "_after_rRNA.fasta"  # after step 1c

INTFNS["gBamFile"] = DIRS["lostReads"] + BASENAME + "_genome.bam"
INTFNS["tBamFile"] = DIRS["lostReads"] + BASENAME + "_transcriptome.bam"

INTFNS["gSamFile"] = DIRS["lostReads"] + BASENAME + "_genome.sam"
INTFNS["tSamFile"] = DIRS["lostReads"] + BASENAME + "_transcriptome.sam"

INTFNS["afterlostReadsFasta"] = DIRS["lostReads"] + BASENAME + "_after_rRNA_lostReads.fasta"  # after step 2
INTFNS["repeatFile"] = DIRS["lostRepeat"] + BASENAME + "_lostRepeats_blastFormat6.tsv"
INTFNS["afterlostRepeatFasta"] = DIRS["lostRepeat"] + BASENAME + "_after_lostRepeat.fasta"  # after step 3
INTFNS["NCL_CIRI_file"] = DIRS["NCL"] + BASENAME +  "_NCL_CIRI_after_bwa.sam"
INTFNS["after_NCL_CIRI_file_prefix"] = BASENAME + "_circRNA.tsv"
INTFNS["afterNCLFasta"] = DIRS["NCL"] + BASENAME + "_after_NCL.fasta"  # after step 4
INTFNS["imrep"] = DIRS["antibody"] + "full_cdr3_" + BASENAME + ".txt"
#INTFNS["igkFile"] = DIRS["igk"] + BASENAME + "_IGK_igblast.tsv"
#INTFNS["iglFile"] = DIRS["igl"] + BASENAME + "_IGL_igblast.tsv"
#INTFNS["tcraFile"] = DIRS["tcra"] + BASENAME + "_TCRA_igblast.tsv"
#INTFNS["tcrbFile"] = DIRS["tcrb"] + BASENAME + "_TCRB_igblast.tsv"
#INTFNS["tcrdFile"] = DIRS["tcrd"] + BASENAME + "_TCRD_igblast.tsv"
#INTFNS["tcrgFile"] = DIRS["tcrg"] + BASENAME + "_TCRG_igblast.tsv"
#INTFNS["immuneAlignments"] = DIRS["antibody"] + BASENAME + "_immuneAlignments.vdjca"
#INTFNS["immuneAligned"] = DIRS["antibody"] + BASENAME + "_immuneAligned.txt"
#INTFNS["afterImmuneFastq"] = DIRS["antibody"] + BASENAME + "_afterImmune.fastq"
INTFNS["afterImmuneFasta"] = DIRS["antibody"] + BASENAME + "_afterImmune.fasta"  # after step 5
INTFNS["metaphlan_intermediate_map"] = DIRS["metaphlan"] + BASENAME + "_metaphlan.map"
INTFNS["metaphlan_intermediate_bowtie2out"] =  DIRS["metaphlan"] + BASENAME + "_bowtie2out.txt"
INTFNS["metaphlan_output"] =  DIRS["metaphlan"]  +  BASENAME  +  "_metaphlan_output.tsv"
INTFNS["afterBacteriaFasta"] = DIRS["bacteria"] + BASENAME + "_afterBacteria.fasta"
INTFNS["afterVirusFasta"] = DIRS["virus"] + BASENAME + "_afterVirus.fasta"
INTFNS["afterFungiFasta"] = DIRS["fungi"] + BASENAME + "_afterFungi.fasta"
INTFNS["afterProtozoaFasta"] = DIRS["protozoa"] + BASENAME + "_afterProtozoa.fasta"


INTFNS["sam_viral"] = DIRS["virus"] + BASENAME + ".viral.sam"
INTFNS["sam_viral_vipr"] = DIRS["virus"] + BASENAME + ".viral.vipr.sam"
INTFNS["sam_fungi"] = DIRS["fungi"] + BASENAME + ".fungi.sam"
INTFNS["sam_protozoa"] = DIRS["protozoa"] + BASENAME + ".protozoa.sam"

INTFNS["bam_viral"] = DIRS["virus"] + BASENAME + ".viral.bam"
INTFNS["bam_viral_vipr"] = DIRS["virus"] + BASENAME + ".viral.vipr.bam"
INTFNS["bam_fungi"] = DIRS["fungi"] + BASENAME + ".fungi.bam"
INTFNS["bam_protozoa"] = DIRS["protozoa"] + BASENAME + ".protozoa.bam"


INTFNS["unaccountedReadsFasta"] = ARGS.dir + "/" + BASENAME + "_unaccountedReads.fasta"  # after step 6

################################################################################
# Log file names

LOGFNS = {}
LOGFNS["gLogfile"] = ARGS.dir + "/" + BASENAME + ".log"
LOGFNS["cmdLogfile"] = ARGS.dir + "/" + "dev.log"
LOGFNS["toolsLogfile"] = ARGS.dir + "/"+"tools.log"
LOGFNS["logQC"] = DIRS["QC"] + BASENAME + "_QC.log"  # step 1b
LOGFNS["logrRNA"] = DIRS["QC"] + BASENAME + "_rRNA.log"  # step 1c
LOGFNS["logHuman"] = DIRS["lostReads"] + BASENAME + "_lostReads.log"  # step 2
LOGFNS["log_bowtieWG"] = DIRS["lostReads"] + BASENAME + "_bowtieWG.log"
LOGFNS["log_bowtieTR"] = DIRS["lostReads"] + BASENAME + "_bowtieTR.log"
LOGFNS["logLostRepeat"] = DIRS["lostRepeat"] + BASENAME + "_lostRepeat.log"  # step 3
LOGFNS["logNCL"] = DIRS["NCL"] + BASENAME + "_NCL.log"  # step 4
#LOGFNS["logIGH"] = DIRS["igh"] + BASENAME + "_igh.log"  # step 5
#LOGFNS["logIGL"] = DIRS["igl"] + BASENAME + "_igl.log"
#LOGFNS["logIGK"] = DIRS["igk"] + BASENAME + "_igk.log"
#LOGFNS["logTCRA"] = DIRS["tcra"] + BASENAME + "_tcra.log"
#LOGFNS["logTCRB"] = DIRS["tcrb"] + BASENAME + "_tcrb.log"
#LOGFNS["logTCRG"] = DIRS["tcrg"] + BASENAME + "_tcrg.log"
#LOGFNS["logTCRD"] = DIRS["tcrd"] + BASENAME + "_tcrd.log"
LOGFNS["logMetaphlan"] = DIRS["metaphlan"] + BASENAME + "_metaphlan.log"  # step 6
LOGFNS["logBacteria"] = DIRS["bacteria"] + BASENAME + "_bacteria.log"
LOGFNS["logVirus"] = DIRS["virus"] + BASENAME + "_virus.log"
LOGFNS["logFungi"] = DIRS["fungi"] + BASENAME + "_fungi.log"
LOGFNS["logProtozoa"] = DIRS["protozoa"] + BASENAME + "_protozoa.log"
LOGFNS["bacteriaFile"] = DIRS["bacteria"] + BASENAME + "_bacteria_blastFormat6.csv"
LOGFNS["virusFile"] = DIRS["virus"] + BASENAME + "_virus_blastFormat6.csv"
LOGFNS["bacteriaFileFiltered"] = DIRS["bacteria"] + BASENAME + "_bacteriaFiltered.csv"
LOGFNS["virusFileFiltered"] = DIRS["virus"] + BASENAME + "_virusFiltered.csv"
LOGFNS["fungiFileFiltered"] = DIRS["fungi"] + BASENAME + "_fungiFiltered.csv"
LOGFNS["protozoaFileFiltered"] = DIRS["protozoa"] + BASENAME + "_protozoaFiltered.csv"


################################################################################
# Run file names

RUNFNS = dict()
RUNFNS["runLostReadsFile"] = DIRS["lostReads"] + "/runLostReads_" + BASENAME + ".sh"
RUNFNS["runLostRepeatFile"] = DIRS["lostRepeat"] + "/runLostRepeat_" + BASENAME + ".sh"
RUNFNS["runNCL_CIRIfile"] = DIRS["NCL"] + "/run_NCL_CIRI" + BASENAME + ".sh" 
#RUNFNS["runIGHFile"] = DIRS["igh"] + "/runIGH_" + BASENAME + ".sh"
#RUNFNS["runIGKFile"] = DIRS["igk"] + "/runIGK_" + BASENAME + ".sh"
#RUNFNS["runIGLFile"] = DIRS["igl"] + "/runIGL_" + BASENAME + ".sh"
#RUNFNS["runTCRAFile"] = DIRS["tcra"] + "/runTCRA_" + BASENAME + ".sh"
#RUNFNS["runTCRBFile"] = DIRS["tcrb"] + "/runTCRB_" + BASENAME + ".sh"
#RUNFNS["runTCRDFile"] = DIRS["tcrd"] + "/runTCRD_" + BASENAME + ".sh"
#RUNFNS["runTCRGFile"] = DIRS["tcrg"] + "/runTCRG_" + BASENAME + ".sh"
RUNFNS["runAntibodyFile"] = DIRS["antibody"] + "/runAntibody_" + BASENAME + ".sh"
RUNFNS["runBacteriaFile"] = DIRS["bacteria"] +"/runBacteria_" + BASENAME + ".sh"
RUNFNS["runVirusFile"] = DIRS["virus"] + "/runVirus_" + BASENAME + ".sh"
RUNFNS["runFungiFile"] = DIRS["fungi"] + "/runFungi_" + BASENAME + ".sh"
RUNFNS["runProtozoaFile"] = DIRS["protozoa"] + "/runProtozoa_" + BASENAME + ".sh"

RUNFNS["run_metaphlan_file"] = DIRS["metaphlan"] + "/run_metaphlan_" + BASENAME + ".sh"
