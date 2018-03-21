#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# INTRO
# ------------------------------------------------------------------------------

set -e
echo '--------------------------------------------------------------------------------'
echo 'Read Origin Protocol: Main Program'
echo '--------------------------------------------------------------------------------'
DIR=`dirname $(readlink -e "$0")`
sed '/##/ q' "$DIR/README.md" | head -n -2 | tail -n +3
echo '--------------------------------------------------------------------------------'

# Add MiniConda to PATH if it's available.
if [ -d "$DIR/tools/MiniConda/bin" ]; then
    export PATH="$DIR/tools/MiniConda/bin:$PATH"
fi

# ------------------------------------------------------------------------------
# PARSE OPTIONS
# ------------------------------------------------------------------------------

# Test for getopt availability.
set +e
getopt --test
if [ $? -ne 4 ]; then
    echo "Error: Environment doesn't support getopt." >&2
    exit 1
fi
set -e

# Call getopt.
SHORT_OPTIONS='o:s:abzdfimpqxh'
LONG_OPTIONS='organism:,steps:,fasta,bam,gzip,dev,force,ignore-extensions,max,\
pe,quiet,commands,help'
set +e
PARSED=`getopt --options="$SHORT_OPTIONS" --longoptions="$LONG_OPTIONS" \
--name "$0" -- "$@"`
if [ $? -ne 0 ]; then
    exit 1  # getopt will have printed the error message
fi
set -e
eval set -- "$PARSED"

# Set default options.
ORGANISM='human'
STEPS='rdna reference repeats immune metaphlan viral fungi protozoa'
    # Non-default: lowq (too slow).
    # Disabled: circrna bacteria (databases missing).
FASTA=false
BAM=false
GZIP=false
DEV=false
FORCE=false
IGNORE_EXTENSIONS=false
MAX=''
PE=''
QUIET=false
COMMANDS=false
UNMAPPED_READS=''
OUTPUT_DIR=''

# Review parsed options.
while true; do
    case "$1" in
        -o|--organism)
            # Run for the specified organism instead of human.
            ORGANISM="$2"
            shift 2
            ;;
        -s|--steps)
            # Select the analysis modes to use.
            STEPS=`tr ',' ' ' <<<"$2"`
            shift 2
            ;;
        -a|--fasta)
            # Input unmapped reads in .fasta format instead of .fastq format.
            # Forcibly disables low-quality read filtering.
            FASTA=true
            shift
            ;;
        -b|--bam)
            # Input unmapped reads in .bam format instead of .fastq format.
            BAM=true
            shift
            ;;
        -z|--gzip)
            # gunzip the input file.
            GZIP=true
            shift
            ;;
        -d|--dev)
            # Keep intermediate FASTA files.
            DEV=true
            shift
            ;;
        -f|--force)
            # Overwrite the analysis destination directory.
            FORCE=true
            shift
            ;;
        -i|--ignore-extensions)
            # Ignore incorrect .fastq/.fq/.fasta/.fa file extensions.
            # Does not ignore incorrect .gz/.bam file extensions.
            IGNORE_EXTENSIONS=true
            shift
            ;;
        -m|--max)
            # Use a liberal threshold when remapping to reference.
            MAX='--max'
            shift
            ;;
        -p|--pe)
            # Not implemented (usage unclear).
            # Report the number of discordant read pairs, with reads from the
            # same pair classified into different classes.
            PE='--pe'
            shift
            ;;
        -q|--quiet)
            # Not implemented (usage unclear).
            # Suppress progress report and warnings.
            QUIET=true
            shift
            ;;
        -x|--commands)
            # Print all commands.
            COMMANDS=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [-o ORGANISM] [-s STEPS] [-abz] [-dfimxh]" \
                "unmapped_reads output_dir" >&2
            exit 0
            ;;
        --)
            # Mandatory arguments.
            UNMAPPED_READS="$2"
            OUTPUT_DIR="$3"
            if [ "$UNMAPPED_READS" = '' ] || [ "$OUTPUT_DIR" = '' ]; then
                echo 'Error: Insufficient arguments.'
                exit 1
            fi
            shift 3
            break
            ;;
        *)
            echo "Error parsing options." >&2
            exit 1
            ;;
    esac
done

# Add all steps if selected.
if [ "$STEPS" = 'all' ]; then
    STEPS='lowq rdna reference repeats circrna immune microbiome'
fi

# Convert to absolute paths.
UNMAPPED_READS=`readlink -m "$UNMAPPED_READS"`
OUTPUT_DIR=`readlink -m "$OUTPUT_DIR"`

# Check if UNMAPPED_READS exists.
if [ ! -e "$UNMAPPED_READS" ]; then
    echo "Error: $UNMAPPED_READS doesn't exist." >&2
    exit 1
fi

# Check if OUTPUT_DIR exists, then make it.
if [ -d "$OUTPUT_DIR" ]; then
    if [ $FORCE = true ]; then
        rm -fr "$OUTPUT_DIR"
    else
        echo "Error: The directory $OUTPUT_DIR exists. Please choose a" \
            'different directory in which to save results of the analysis, or' \
            'use the -f option to overwrite the directory.' >&2
        exit 1
    fi
fi
mkdir -p "$OUTPUT_DIR"

# ------------------------------------------------------------------------------
# CONSTANTS
# ------------------------------------------------------------------------------

# Sample name and database location.
SAMPLE=`basename "$UNMAPPED_READS" | sed 's \([^\.]*\)\..* \1 '`
DB="$DIR/db_$ORGANISM"

# Duplicate stdout and stderr to the log file. Print commands if selected.
touch "$OUTPUT_DIR/$SAMPLE--general.log"
exec &> >(tee -i "$OUTPUT_DIR/$SAMPLE--general.log")
if [ $COMMANDS = true ]; then
    set -x
fi
echo "Input file: $UNMAPPED_READS"

# Declare output directories.
declare -A DIRS=(
    ['01_lowq']="$OUTPUT_DIR/01_lowq"
    ['02_rdna']="$OUTPUT_DIR/02_rdna"
    ['03_reference']="$OUTPUT_DIR/03_reference"
    ['04_repeats']="$OUTPUT_DIR/04_repeats"
    ['05_circrna']="$OUTPUT_DIR/05_circrna"
    ['06_immune']="$OUTPUT_DIR/06_immune"
    ['07a_metaphlan']="$OUTPUT_DIR/07a_metaphlan"
    ['07b_bacteria']="$OUTPUT_DIR/07b_bacteria"
    ['07c_viral']="$OUTPUT_DIR/07c_viral"
    ['07d_fungi']="$OUTPUT_DIR/07d_fungi"
    ['07e_protozoa']="$OUTPUT_DIR/07e_protozoa"
)

# Make output directories if necessary.
for val in "${DIRS[@]}"; do
    mkdir -p "$val"
done

# Declare intermediate files.
declare -A INTFNS=(
    # Step 1 (lowq).
    ['01_lowq_post']="${DIRS['01_lowq']}/$SAMPLE--01_lowq_post.fasta"

    # Step 2 (rdna).
    ['02_rdna_output']="${DIRS['02_rdna']}/$SAMPLE--02_rdna_output.bam"
    ['02_rdna_reads']="${DIRS['02_rdna']}/$SAMPLE--02_rdna_reads.txt"
    ['02_rdna_post']="${DIRS['02_rdna']}/$SAMPLE--02_rdna_post.fasta"

    # Step 3 (reference).
    ['03_reference_genomeoutput']="${DIRS['03_reference']}/$SAMPLE--03_reference_genomeoutput.bam"
    ['03_reference_transcriptomeoutput']="${DIRS['03_reference']}\
/$SAMPLE--03_reference_transcriptomeoutput.bam"
    ['03_reference_reads']="${DIRS['03_reference']}/$SAMPLE--03_reference_reads.txt"
    ['03_reference_post']="${DIRS['03_reference']}/$SAMPLE--03_reference_post.fasta"

    # Step 4 (repeats).
    ['04_repeats_output']="${DIRS['04_repeats']}/$SAMPLE--04_repeats_output.tsv"
    ['04_repeats_reads']="${DIRS['04_repeats']}/$SAMPLE--04_repeats_reads.txt"
    ['04_repeats_post']="${DIRS['04_repeats']}/$SAMPLE--04_repeats_post.fasta"

    # Step 5 (circrna).
    ['05_circrna_reads']="${DIRS['05_circrna']}/$SAMPLE--05_circrna_reads.txt"
    ['05_circrna_post']="${DIRS['05_circrna']}/$SAMPLE--05_circrna_post.fasta"

    # Step 6 (immune).
    ['06_immune_output']="${DIRS['06_immune']}/$SAMPLE--06_immune_output.cdr3"
    ['06_immune_clonality']="${DIRS['06_immune']}/$SAMPLE--06_immune_clonality"
    ['06_immune_reads']="${DIRS['06_immune']}/$SAMPLE--06_immune_reads.txt"
    ['06_immune_post']="${DIRS['06_immune']}/$SAMPLE--06_immune_post.fasta"

    # Step 7a (metaphlan).
    # No post file (don't reduce unmapped reads using MetaPhlAn results).
    ['07a_metaphlan_map']="${DIRS['07a_metaphlan']}/$SAMPLE--07a_metaphlan_map.map"
    ['07a_metaphlan_bowtie2out']="${DIRS['07a_metaphlan']}/$SAMPLE--07a_metaphlan_bowtie2out.txt"
    ['07a_metaphlan_output']="${DIRS['07a_metaphlan']}/$SAMPLE--07a_metaphlan_output.tsv"

    # Step 7b (bacteria).
    ['07b_bacteria_output']="${DIRS['07b_bacteria']}/$SAMPLE--07b_bacteria_output.bam"
    ['07b_bacteria_reads']="${DIRS['07b_bacteria']}/$SAMPLE--07b_bacteria_reads.txt"
    ['07b_bacteria_post']="${DIRS['07b_bacteria']}/$SAMPLE--07b_bacteria_post.fasta"

    # Step 7c (viral).
    ['07c_viral_output']="${DIRS['07c_viral']}/$SAMPLE--07c_viral_output.bam"
    ['07c_viral_viproutput']="${DIRS['07c_viral']}/$SAMPLE--07c_viral_viproutput.bam"
    ['07c_viral_reads']="${DIRS['07c_viral']}/$SAMPLE--07c_viral_reads.txt"
    ['07c_viral_post']="${DIRS['07c_viral']}/$SAMPLE--07c_viral_post.fasta"

    # Step 7d (fungi).
    ['07d_fungi_output']="${DIRS['07d_fungi']}/$SAMPLE--07d_fungi_output.bam"
    ['07d_fungi_reads']="${DIRS['07d_fungi']}/$SAMPLE--07d_fungi_reads.txt"
    ['07d_fungi_post']="${DIRS['07d_fungi']}/$SAMPLE--07d_fungi_post.fasta"

    # Step 6e (protozoa).
    ['07e_protozoa_output']="${DIRS['07e_protozoa']}/$SAMPLE--07e_protozoa_output.bam"
    ['07e_protozoa_reads']="${DIRS['07e_protozoa']}/$SAMPLE--07e_protozoa_reads.txt"
    ['07e_protozoa_post']="${DIRS['07e_protozoa']}/$SAMPLE--07e_protozoa_post.fasta"

    ['unaccounted']="$OUTPUT_DIR/$SAMPLE--unaccounted.fasta"
)

# Declare log files.
declare -A LOGFNS=(
    #['01_lowq']="${DIRS['01_lowq']}/$SAMPLE--01_lowq.log"
    ['02_rdna']="${DIRS['02_rdna']}/$SAMPLE--02_rdna.log"
    ['03_reference']="${DIRS['03_reference']}/$SAMPLE--03_reference.log"
    ['04_repeats']="${DIRS['04_repeats']}/$SAMPLE--04_repeats.log"
    ['05_circrna']="${DIRS['05_circrna']}/$SAMPLE--05_circrna.log"
    ['06_immune']="${DIRS['06_immune']}/$SAMPLE--06_immune.log"
    ['07a_metaphlan']="${DIRS['07a_metaphlan']}/$SAMPLE--07a_metaphlan.log"
    ['07b_bacteria']="${DIRS['07b_bacteria']}/$SAMPLE--07b_bacteria.log"
    ['07c_viral']="${DIRS['07c_viral']}/$SAMPLE--07c_viral.log"
    ['07d_fungi']="${DIRS['07d_fungi']}/$SAMPLE--07d_fungi.log"
    ['07e_protozoa']="${DIRS['07e_protozoa']}/$SAMPLE--07e_protozoa.log"
    ['counts']="$OUTPUT_DIR/$SAMPLE--counts.csv"
)

# ------------------------------------------------------------------------------
# UTILITY
# ------------------------------------------------------------------------------

reads_present () {
    if [ `wc -l <"$1"` -le 1 ]; then
        echo 'No more reads!'
        return 1  # false
    else
        return 0  # true
    fi
}

clean () {
    if [ $DEV = false ]; then
        rm "$1"
    fi
}

# ------------------------------------------------------------------------------
# PREPROCESS INPUT FILE
# ------------------------------------------------------------------------------

# Do preprocessing in OUTPUT_DIR.
cd "$OUTPUT_DIR"
cp "$UNMAPPED_READS" .
current=`basename "$UNMAPPED_READS"`

# Unpack the unmapped reads if -z and/or -b are selected.
if [ $GZIP == true ]; then
    echo 'Unpacking gzip...'
    post=`basename "$current" .gz`
    if [ "$post" == "$current" ]; then
        echo 'Error: input file missing .gz extension' >&2
        exit 1
    fi
    gunzip "$current"
    current="$post"
fi
if [ $BAM == true ]; then
    echo 'Unpacking bam...'
    post="$(basename "$current" .bam).fastq"
    if [ "$post" == "$current.fastq" ]; then
        echo 'Error: input file missing .bam extension' >&2
        exit 1
    fi
    samtools bam2fq "$current" >"$post"
    clean "$current"
    current="$post"
fi

# Inspect the input file, then restore current to a full path.
if [ $FASTA == true ]; then
    if [ $IGNORE_EXTENSIONS = false ] && \
        [ `basename $current .fasta` == "$current" ] && \
        [ `basename $current .fa` == "$current" ]; then
        echo 'Error: input file missing .fasta/.fa extension' >&2
        exit 1
    fi
    N=`grep -c '^>' "$current"`
    READ_LENGTH=$(($(sed -n '2 p' <"$current" | wc -m) - 1))
else
    if [ $IGNORE_EXTENSIONS = false ] && \
        [ `basename $current .fastq` == "$current" ] && \
        [ `basename $current .fq` == "$current" ]; then
        echo 'Error: input file missing .fastq/.fq extension' >&2
        exit 1
    fi
    line_count=`wc -l <"$current"`
    N=`bc <<<"$line_count/4"`
    READ_LENGTH=$(($(sed -n '2 p' <"$current" | wc -m) - 1))
fi
echo "Processing $N unmapped reads. The first unmapped read has length $READ_LENGTH."
current=`readlink -e "$current"`

# Record the number of reads accounted for in each step.
declare -A n_reads=(
    ['01_lowq']=0
    ['02_rdna']=0
    ['03_reference']=0
    ['04_repeats']=0
    ['05_circrna']=0
    ['05_immune']=0
    ['07b_bacteria']=0
    ['07c_viral']=0
    ['07d_fungi']=0
    ['07e_protozoa']=0
)

# ------------------------------------------------------------------------------
# 1. LOW QUALITY READ MARKING
# ------------------------------------------------------------------------------

echo "1. Low quality read marking (-s lowq)..."
cd "${DIRS['01_lowq']}"
post="${INTFNS['01_lowq_post']}"

if ! grep -q 'lowq' <<<"$STEPS" || [ $FASTA = true ] || \
    ! reads_present "$current"; then
    echo '--> Skipped low quality read marking.'
    
    # Must convert to fasta to continue.
    if [ $FASTA = false ]; then
        fastq_to_fasta -n <"$current" >"$post"
        clean "$current"
        current="$post"
    fi
else
    n_reads['01_lowq']=`python "$DIR/helper.py" lowq $MAX $PE \
        --pre "$current" --post "$post"`
    echo "--> Marked lowq in the names of ${n_reads['01_lowq']} low quality" \
        'reads.'
    echo '    These reads are not filtered.'
    clean "$current"
    current="$post"
fi

# ------------------------------------------------------------------------------
# 2. rDNA PROFILING
# ------------------------------------------------------------------------------

echo "2. rDNA profiling (-s rdna)..."
cd "${DIRS['02_rdna']}"
post="${INTFNS['02_rdna_post']}"

if ! grep -q 'rdna' <<<"$STEPS" || ! reads_present "$current"; then
    echo '--> Skipped rDNA profiling.'
else
    bowtie2 -f -x "$DB/ribosomal.DNA/ribosomal.DNA" --end-to-end -D 15 -R 2 \
        -L 22 -i S,1,1.15 "$current" 2>"${LOGFNS['02_rdna']}" \
        | samtools sort - >${INTFNS['02_rdna_output']}
    samtools index "${INTFNS['02_rdna_output']}"
    n_reads['02_rdna']=`python "$DIR/helper.py" rdna $MAX $PE \
        -i "${INTFNS['02_rdna_output']}" \
        -o ${INTFNS['02_rdna_reads']} \
        --pre "$current" --post "$post"`
    echo "--> Filtered ${n_reads['02_rdna']} reads from ribosomal DNA."
    clean "$current"
    current="$post"
fi

# ------------------------------------------------------------------------------
# 3. REMAPPING TO REFERENCE
# ------------------------------------------------------------------------------

echo '3. Remapping to reference (-s reference)...'
cd "${DIRS['03_reference']}"
post="${INTFNS['03_reference_post']}"

if ! grep -q 'reference' <<<"$STEPS" || ! reads_present "$current"; then
    echo '--> Skipped remapping to reference.'
else
    bwa mem "$DB/BWA.index/genome.fa" "$current" \
        2>>"${LOGFNS['03_reference']}" \
        | samtools sort - >"${INTFNS['03_reference_genomeoutput']}"
    bwa mem "$DB/BWA.index/isoforms_GRCh38_Ensembl.fasta" "$current" \
        2>>"${LOGFNS['03_reference']}" \
        | samtools sort - >"${INTFNS['03_reference_transcriptomeoutput']}"
    samtools index "${INTFNS['03_reference_genomeoutput']}"
    samtools index "${INTFNS['03_reference_transcriptomeoutput']}"
    n_reads['03_reference']=`python "$DIR/helper.py" reference $MAX $PE \
        -i "${INTFNS['03_reference_genomeoutput']},${INTFNS['03_reference_transcriptomeoutput']}" \
        -o "${INTFNS['03_reference_reads']}" \
        --pre "$current" --post "$post"`
    echo "--> Filtered ${n_reads['03_reference']} reads from reference genome" \
        'or transcriptome.'
    clean "$current"
    current="$post"
fi

# ------------------------------------------------------------------------------
# 4. REPEAT PROFILING
# ------------------------------------------------------------------------------

echo '4. Repeat profiling (-s repeats)...'
cd "${DIRS['04_repeats']}"
post="${INTFNS['04_repeats_post']}"

if ! grep -q 'repeats' <<<"$STEPS" || ! reads_present "$current"; then
    echo '--> Skipped repeat profiling.'
else
    blastn -task megablast -index_name "$DB/repeats/repbase.fa" \
        -use_index true -query "$current" -db "$DB/repeats/repbase.fa" \
        -outfmt 6 -evalue 1e-05 >"${INTFNS['04_repeats_output']}" \
        2>"${LOGFNS['04_repeats']}"
    n_reads['04_repeats']=`python "$DIR/helper.py" repeats $MAX $PE \
        -i "${INTFNS['04_repeats_output']}" \
        -o "${INTFNS['04_repeats_reads']}" \
        --pre "$current" --post "$post"`
    echo "--> Filtered ${n_reads['04_repeats']} reads from repeat sequences."
    clean "$current"
    current="$post"
fi

# ------------------------------------------------------------------------------
# 5. CIRCULAR RNA PROFILING
# ------------------------------------------------------------------------------

echo '5. Circular RNA profiling (-s circrna)...'
cd "${DIRS['05_circrna']}"
post="${INTFNS['05_circrna_post']}"

# Disabled (database missing).
#if ! grep -q 'circrna' <<<"$STEPS" || ! reads_present "$current"; then
if true; then
    echo '--> Skipped circular RNA profiling.'
else  # WARNING: This branch is untested!
    tophat2 -o . --fusion-search --keep-fasta-order --no-coverage-search \
        "$DB/Bowtie2Index/genome" "$current" 2>"${LOGFNS['05_circrna']}"
    samtools bam2fq 'accepted_hits.bam' >'accepted_hits.fastq'
    n_reads['05_circrna']=`python "$DIR/helper.py" circrna $MAX $PE \
        -i 'accepted_hits.fastq' \
        -o "${INTFNS['05_circrna_reads']}" \
        --pre "$current" --post "$post"`
    echo "--> Filtered ${n_reads['05_circrna']} reads from circular RNA."
    clean "$current"
    current="$post"
fi

# ------------------------------------------------------------------------------
# 6. IMMUNE PROFILING
# ------------------------------------------------------------------------------

echo '6. Immune profiling (-s immune)...'
cd "${DIRS['06_immune']}"
post="${INTFNS['06_immune_post']}"

if ! grep -q 'immune' <<<"$STEPS" || ! reads_present "$current"; then
    echo '--> Skipped immune profiling.'
else
    python "$DIR/tools/imrep/imrep.py" -f -1 --extendedOutput "$current" \
        "${INTFNS['06_immune_output']}" &>"${LOGFNS['06_immune']}"
    python "$DIR/tools/imrep/clonality.py" \
        "${INTFNS['06_immune_output']}" \
        "${INTFNS['06_immune_clonality']}" &>>"${LOGFNS['06_immune']}"
    n_reads['06_immune']=`python "$DIR/helper.py" immune $MAX $PE \
        -i $(ls full_cdr3_*) \
        -o "${INTFNS['06_immune_reads']}" \
        --pre "$current" --post "$post"`
    echo "--> Filtered ${n_reads['06_immune']} reads from T and B cell"\
        'repetoires.'
    clean "$current"
    current="$post"
fi

# ------------------------------------------------------------------------------
# 7. MICROBIOME PROFILING
# ------------------------------------------------------------------------------

echo '7a. MetaPhlAn profiling (-s metaphlan)...'
cd "${DIRS['07a_metaphlan']}"
# No post file (don't reduce unmapped reads using MetaPhlAn results).

if ! grep -qE 'metaphlan|microbiome' <<<"$STEPS" || ! reads_present "$current"; then
    echo '--> Skipped MetaPhlAn profiling.'
else
    python "$DIR/tools/metaphlan2/metaphlan2.py" "$current" \
        --input_type multifasta --nproc 8 \
        --bowtie2out "${INTFNS['07a_metaphlan_bowtie2out']}" \
        >"${INTFNS['07a_metaphlan_output']}" 2>"${LOGFNS['07a_metaphlan']}"
    n_reads_07a_metaphlan=`wc -l <"${INTFNS['07a_metaphlan_output']}"`
    echo "--> Identified $n_reads_07a_metaphlan reads using MetaPhlAn."
    echo '    These reads are neither filtered nor included in the total.'
    # Don't clean or change $current.
fi

echo '7b. Bacterial profiling (-s bacteria)...'
cd "${DIRS['07b_bacteria']}"
post="${INTFNS['07b_bacteria_post']}"

# Disabled (database missing).
#if ! grep -qE 'bacteria|microbiome' <<<"$STEPS" || ! reads_present "$current"; then
if true; then
    echo '--> Skipped bacterial profiling.'
else
    bwa mem "$DB/bacteria/bacteria.ncbi.february.3.2018.fasta" "$current" \
        | samtools sort - >"${INTFNS['07b_bacteria_output']}"
    n_reads['07b_bacteria']=`python "$DIR"/helper.py microbiome $MAX $PE \
        -i "${INTFNS['07b_bacteria_output']}" \
        -o "${INTFNS['07b_bacteria_reads']}" \
        --pre "$current" --post "$post"`
    echo "--> Filtered ${n_reads['07b_bacteria']} reads from bacterial genomes."
    clean "$current"
    current="$post"
fi

echo '7c. Viral profiling (-s viral)...'
cd "${DIRS['07c_viral']}"
post="${INTFNS['07c_viral_post']}"

if ! grep -qE 'viral|microbiome' <<<"$STEPS" || ! reads_present "$current"; then
    echo '--> Skipped viral profiling.'
else
    bwa mem -a "$DB/viral/viral.ncbi.february.3.2018.fasta" "$current" \
        2>"${LOGFNS['07c_viral']}" \
        | samtools sort - >"${INTFNS['07c_viral_output']}"
    bwa mem -a "$DB/viral.vipr/NONFLU_All.fastq" "$current" \
        2>>"${LOGFNS['07c_viral']}" \
        | samtools sort - >"${INTFNS['07c_viral_viproutput']}"
    samtools index "${INTFNS['07c_viral_output']}"
    samtools index "${INTFNS['07c_viral_viproutput']}"
    n_reads['07c_viral']=`python "$DIR/helper.py" microbiome $MAX $PE \
        -i "${INTFNS['07c_viral_output']},${INTFNS['07c_viral_viproutput']}" \
        -o "${INTFNS['07c_viral_reads']}" \
        --pre "$current" --post "$post"`
    echo "--> Filtered ${n_reads['07c_viral']} reads from viral genomes."
    clean "$current"
    current="$post"
fi

echo '7d. Fungal profiling (-s fungi)...'
cd "${DIRS['07d_fungi']}"
post="${INTFNS['07d_fungi_post']}"

if ! grep -qE 'fungi|microbiome' <<<"$STEPS" || ! reads_present "$current"; then
    echo '--> Skipped fungal profiling.'
else
    bwa mem -a "$DB/fungi/fungi.ncbi.february.3.2018.fasta" "$current" \
        2>"${LOGFNS['07d_fungi']}" \
        | samtools sort - >"${INTFNS['07d_fungi_output']}"
    samtools index "${INTFNS['07d_fungi_output']}"
    n_reads['07d_fungi']=`python "$DIR/helper.py" microbiome $MAX $PE \
        -i "${INTFNS['07d_fungi_output']}" \
        -o "${INTFNS['07d_fungi_reads']}" \
        --pre "$current" --post "$post"`
    echo "--> Filtered ${n_reads['07d_fungi']} reads from fungal genomes."
    clean "$current"
    current="$post"
fi

echo '7e. Protozoan profiling (-s protozoa)...'
cd "${DIRS['07e_protozoa']}"
post="${INTFNS['07e_protozoa_post']}"

if ! grep -qE 'protozoa|microbiome' <<<"$STEPS" || ! reads_present "$current"; then
    echo '--> Skipped protozoan profiling.'
else
    bwa mem -a "$DB/protozoa/protozoa.ncbi.february.3.2018.fasta" "$current" \
        2>"${LOGFNS['07e_protozoa']}" \
        | samtools sort - >"${INTFNS['07e_protozoa_output']}"
    samtools index "${INTFNS['07e_protozoa_output']}"
    n_reads['07e_protozoa']=`python "$DIR/helper.py" microbiome $MAX $PE \
        -i "${INTFNS['07e_protozoa_output']}" \
        -o "${INTFNS['07e_protozoa_reads']}" \
        --pre "$current" --post "$post"`
    echo "--> Filtered ${n_reads['07e_protozoa']} reads from protozoan genomes."
    clean "$current"
    current="$post"
fi

# ------------------------------------------------------------------------------
# CLEANUP
# ------------------------------------------------------------------------------

# Revise low quality read count.
if grep -q 'lowq' <<<"$STEPS"; then
    n_reads['01_lowq']=`grep -c '^>lowq_' "$current"`
fi

# Sum accounted reads and write to file.
steps=''
counts=''
sum=0
for key in "${!n_reads[@]}"; do
    steps+="$key,"
    counts+="${n_reads[$key]},"
    ((sum += ${n_reads[$key]})) || true
done
steps=`sed 's .$  ' <<<"$steps"`
counts=`sed 's .$  ' <<<"$counts"`
pct=`bc -l <<<"scale=2; 100*$sum/$N"`
echo "$steps" >"${LOGFNS['counts']}"
echo "$counts" >>"${LOGFNS['counts']}"

# Unaccounted reads.
cp "$current" "${INTFNS['unaccounted']}"
clean "$current"

# Final message.
echo "Summary: ROP has accounted for $sum reads ($pct% of unmapped reads)."
if grep -q 'lowq' <<<"$STEPS"; then
    echo "Of those, ${n_reads['01_lowq']} were considered accounted due to low quality."
fi
