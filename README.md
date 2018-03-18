# Read Origin Protocol

The Read Origin Protocol (ROP) is a computational protocol that aims to
discover the source of all reads, including those originating from complex RNA
molecules, recombinant antibodies, and microbial communities. 

Written by:

- Serghei Mangul (<smangul@ucla.edu>)
- Harry Taegyun Yang (<harry2416@gmail.com>)
- Kevin Hsieh (<kevin.hsieh@ucla.edu>)
- Linus Chen (<u6.30cl@gmail.com>)

at the University of California, Los Angeles (UCLA). 

Released under the terms of the General Public License, version 3.0 (GPLv3).
For more information, please visit: <https://github.com/smangul1/rop/wiki>

## Installing ROP

To install ROP, first clone this repository, then run

```
./install.sh
```

from the repository's directory. This will download dependencies and databases.
The default installation will generally suffice, but the following options are
available:

- `-c|--clean`: Just remove installed tools.
    - The installation script then must be re-run in order to use ROP again.
- `-n|--native`: Use native python.
    - MiniConda will not be downloaded.
    - May lead to dependency errors.
- `-f|--force`: Unlink databases.
    - Use with caution.
- `-l|--link LINK`: Link databases instead of downloading.
    - Useful if you previously downloaded an ROP database.
    - A symlink will be created in the current directory.
- `-d|--db-dest DB_DEST` (default: `.`): Change database download location.
    - Useful for managing space.
    - A symlink will be created in the current directory.
- `-o|--organism ORGANISM` (default: `human`): Organism to download databases
  for.
- `-s|--select-db SELECT_DB` (default: all): Database(s) to download for the
  specified organism.
    - A comma-separated list of one or more of the following: repeat, immune,
      microbiome metaphlan, viral, fungi, protozoa.
- `-h|--help`: Displays usage information.

## Using ROP

To use ROP, run

```
rop.sh unmapped_reads output_dir
```

Unless otherwise specified using an option, `unmapped_reads`
must be a .fastq/.fq file, and `output_dir` must not exist (it will be created).
Results will be written to `output_dir`, with one subdirectory for every stage
of the pipeline. The following options are available:

- `-o|--organism` (default: `human`): Run for the specified organism instead of
  human.
- `-s|--steps` (default: all except lowq and bacteria): Select the analysis modes to use.
    - A comma-separated list of one or more of the following: lowq, rdna,
      reference, repeats, circrna, immune, microbiome (which may be subdivided
      into bacteria, metaphlan, viral, fungi, protozoa).
    - `-s all` selects everything.
    - circrna and bacteria are not available in this release.
- `-m|--max`: Use a liberal threshold when remapping to reference.
    - May account for more reads.
- `-f|--force`: Overwrite the analysis destination directory.
- `-d|--dev`: Keep intermediate FASTA files.
    - Consumes extra space.
- `-z|--gzip`: gunzip the input file.
- `-b|--bam`: Input unmapped reads in .bam format instead of .fastq format.
- `-a|--fasta`: Input unmapped reads in .fasta format instead of .fastq format.
  Forcibly disables low-quality read filtering.
- `-h|--help`: Displays usage information.

A small example file is included in the repository in various formats. To try it
out, run one of the following commands from the repository directory:

```
rop.sh -b example/example.bam ropout
rop.sh example/example.fastq ropout
rop.sh -z example/example.fastq.gz ropout
rop.sh -a example/example.fasta ropout
rop.sh -az example/example.fasta.gz ropout
```

Then, browser to the `ropout` directory to see the analysis results!
