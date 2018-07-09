# Read Origin Protocol

The Read Origin Protocol (ROP) is a computational protocol that aims to
discover the source of all reads, including those originating from complex RNA
molecules, recombinant antibodies, and microbial communities. 

Written by:

- Serghei Mangul (<smangul@ucla.edu>)
- Kevin Hsieh (<kevin.hsieh@ucla.edu>)
- Linus Chen (<u6.30cl@gmail.com>)
- Harry Taegyun Yang (<harry2416@gmail.com>)

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
    - To reinstall, use `-r|--reinstall` instead.
    - Overrides conflicting options.
- `-f|--force`: Unlink databases.
    - Use with caution.
- `-n|--native`: Use native python.
    - MiniConda will not be downloaded.
    - You may use `environment.yml` to set up your python environment.
- `-l|--link LINK`: Link databases instead of downloading.
    - Useful if you previously downloaded an ROP database.
    - A symlink will be created in the current directory.
    - Overrides conflicting options.
- `-d|--db-dest DB_DEST` (default: `.`): Change database download location.
    - Useful for managing space.
    - A symlink will be created in the current directory.
- `-o|--organism ORGANISM` (default: `human`): Organism to download databases
  for.
    - Exactly one of the following: human, mouse.
- `-r|--reinstall`: Reinstall tools, even if they're already present.
- `-s|--select-db SELECT_DB` (default: `all`): Database(s) to download for the
  specified organism.
    - A comma-separated list of one or more of the following: basic, repeats,
      microbiome (which may be subdivided into metaphlan, viral, fungi,
      protozoa).
    - `-s all` selects everything.
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

- `-o|--organism` (default: `human`): Organism to run ROP for.
    - Exactly one of the following: human, mouse.
- `-s|--steps` (default: all except lowq): Select the analysis modes to use.
    - A comma-separated list of one or more of the following: lowq, rdna,
      reference, repeats, circrna, immune, microbiome (which may be subdivided
      into metaphlan, bacteria, viral, fungi, protozoa).
    - circrna, metaphlan, and bacteria are not available in the current release.
    - `-s all` selects everything.
- `-a|--fasta`: Input unmapped reads in .fasta format instead of .fastq format.
  Forcibly disables low-quality read filtering.
- `-b|--bam`: Input unmapped reads in .bam format instead of .fastq format.
- `-z|--gzip`: gunzip the input file.
- `-d|--dev`: Keep intermediate FASTA files.
    - Consumes extra space.
- `-f|--force`: Overwrite the analysis destination directory.
- `-i|--ignore-extensions`: Ignore incorrect .fastq/.fq/.fasta/.fa file
  extensions. Does not ignore incorrect .gz/.bam file extensions.
- `-m|--max`: Use a liberal threshold when remapping to reference.
    - May account for more reads.
- `-x|--commands`: Print all commands (diagnostic mode).
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

Then, browse to the `ropout` directory to see the analysis results!
