Below we explain how to prepare the ROP database for a new organism using the mouse as an example. Prepared by 
Linus Chen (u6.30cl@gmail.com) and Kevin Hsieh(kevin.hsieh@ucla.edu). 

###1. Quality Control (folder: rRNA)
**Files:** `rRNA.fa, rRNA[.nhd/.nhi/.nhd/.nin/.nnd/.nni/.nog/.nsd/.nsi/.nsq/(.shd)/.00.idx]`  
The only data needed for this step is ribosomal DNA data. Such data can be downloaded from the NCBI GenBank (http://www.ncbi.nlm.nih.gov/nuccore/BK000964.3) and then constructed into a BLAST database and index using these commands:  
  
	makeblastdb -in rRNA.fa -dbtype nucl -parse_seqids -hash_index -out rRNA -title rRNA -max_file_sz 1GB
	makembindex -input rRNA -output  rRNA -iformat blastdb

###2. Lost Reads (folder: bowtie2Index)
**Files:** `genome.fa, genome[.1.bt2/.2.bt2/.3.bt2/.4.bt2/.rev.1.bt2/.rev.2.bt2]`  
Genomes and indices can be downloaded directly from Ensembl (https://ccb.jhu.edu/software/tophat/igenomes.shtml).  

###3. Lost Repeats (folder: repeats)
**Files:** `repbase.fa.raw, repbase.fa, repbase.fa[.nhr/.nin/.nog/.nsd/.nsi/.nsq/(.shd)]`  
The .fa.raw can be downloaded from RepBase (http://www.girinst.org/repbase/). The BLAST database and index are then constructed as follows:  

	cat repbase.fa.raw | sed 's/\t/-/g' | sed 's/x/N/g' >repbase.fa
	makeblastdb -parse_seqids -dbtype nucl -in repbase.fa
	makembindex -input repbase.fa -output repbase.fa -iformat blastdb

Note that we give the downloaded file the .fa.raw extension to distinguish it from the filtered file.

###4. Non-co-linear RNA (folder: BWAIndex)  
**Files:** `genome.fa (same as in db_mouse/bowtie2Index), genome.fa[.amb/.ann/.bwt/.pac/.sa]`  
Genomes and indices can be downloaded directly from Ensembl (https://ccb.jhu.edu/software/tophat/igenomes.shtml).  

###5. Lymphocytes (VDJ Recombinations) (folder: antibody)
**Directories:** `internal_data, optional_file`  
**Files:** `[IGHD/IGHJ/IGHV/IGKJ/IGKV/IGLJ/IGLV/TRAJ/TRAV/TRBD/TRBJ/TRBV/TRDJ/TRDV/TRGJ/TRGV][.fa.raw/.fa], [IGHD/IGHJ/IGHV/IGKJ/IGKV/IGLJ/IGLV/TRAJ/TRAV/TRBD/TRBJ/TRBV/TRDJ/TRDV/TRGJ/TRGV].fa[.nhr/.nin/.nog/.nsd/.nsi/.nsq/.gNames]`  
The contents of internal_data and optional_file can be downloaded from ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/ for the desired species, preserving directory structure (see the IgBlast README for more information about these directories). The .fa.raw files can be downloaded from IMGT (http://www.imgt.org/vquest/refseqh.html) for the desired species. Upon downloading, this script can be applied to generate the database files:  

	for f in *raw
	do 
		f=`printf $f | sed "s \..*  "`
		echo "--------------------------------------------------------------------------------"
		echo "> PROCESSING $f"
		echo "--------------------------------------------------------------------------------"
		cat $f.fa.raw | 
			tr ">" "\0" |  # use > as a delimiter for grep
			grep -z "^$\|Mus musculus" |  # only use Mus musculus data
			head -c -1 |  # grep appends an extra null byte to the end
			tr "\0" ">" |  # undo previous tr
			sed "s >[^|]*|\([^|]*\)|.* >\1 " |  # remove extraneous data
			tr [:lower:] [:upper:] > $f.fa  # convert to uppercase
		cat $f.fa | grep ">" | sed "s >\(.*\) \1 " > $f.gNames  # write seq_ids to file
		makeblastdb -parse_seqids -dbtype nucl -in $f.fa  # make db
	done
 
Note that we give the downloaded files the .fa.raw extension to distinguish them from filtered files. We remove entries that do not pertain to "Mus musculus", the species we are working with. Then, we remove irrelevant FASTA fields and capitalize all letters. Finally, we run `makeblastdb`.   

###6. Microbiome (folders: bacteria, eupathdb, metaphlan, virus)
The contents of these folders are static and can be copied from the existing human or mouse database.