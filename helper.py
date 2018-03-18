#!/u/project/zarlab/kahsieh/rop-wip/tools/MiniConda/bin/python

"""
--------------------------------------------------------------------------------
Read Origin Protocol: Helper Script
For more information, please visit: <https://github.com/smangul1/rop/wiki>
--------------------------------------------------------------------------------
"""

import argparse
import csv
import pysam
from Bio import SeqIO
import sys

# ------------------------------------------------------------------------------
# PARSE ARGUMENTS
# ------------------------------------------------------------------------------

ap = argparse.ArgumentParser('Read Origin Protocol: Helper Script')

req_args = ap.add_argument_group('Required Arguments')
req_args.add_argument('step')

opt_args = ap.add_argument_group('Optional Arguments')
req_args.add_argument('--input_file', '-i', action='store')
opt_args.add_argument('--output_file', '-o', action='store')
opt_args.add_argument('--pre', action='store')
opt_args.add_argument('--post', action='store')
opt_args.add_argument('--max', action='store_true')
opt_args.add_argument('--pe', action='store_true')

ARGS = ap.parse_args()

# ------------------------------------------------------------------------------
# UTILITY
# ------------------------------------------------------------------------------

# Excludes reads from the pre fasta and outputs the result to the post fasta.
def exclude(pre_file, reads, post_file):
    with open(pre_file) as pre, open(post_file, 'w') as post:
        for seq in SeqIO.parse(pre, 'fasta'):
            if seq.name not in reads:
                SeqIO.write([seq], post, 'fasta')

def writeset(set_o, output_file):
    with open(output_file, 'w') as output_fo:
        for el in sorted(set_o):
            output_fo.write(el + '\n')

# ------------------------------------------------------------------------------
# 1. LOW QUALITY READ MARKING
# ------------------------------------------------------------------------------

if ARGS.step == 'lowq':
    count = 0
    with open(ARGS.post, 'w') as post:
        for record in SeqIO.parse(ARGS.pre, 'fastq'):
            j = record.letter_annotations['phred_quality']
            prc = sum(1 for i in j if i >= 20) / float(len(j))
            read_length = len(str(record.seq))
            if prc > 0.75 and all(i in 'ACTGN' for i in record.seq):
                post.write(str('>' +\
                    record.name.replace('/', '---') + '_length_' +\
                    str(read_length)) + '\n')
                post.write(str(record.seq) + '\n')
            else:
                post.write(str('>lowq_' +\
                    record.name.replace('/', '---') + '_length_' +\
                    str(read_length)) + '\n')
                post.write(str(record.seq) + '\n')
                count += 1
    print(count)

# ------------------------------------------------------------------------------
# 2. rDNA PROFILING
# ------------------------------------------------------------------------------

elif ARGS.step == 'rdna':
    reads = set()
    with pysam.AlignmentFile(ARGS.input_file, 'rb', check_sq=False) as input_fo:
        for read in input_fo.fetch():
            number_mismatches = int(read.get_tag('NM'))
            read_length = int(read.infer_read_length())
            identity = 1 - number_mismatches/read_length
            if identity > 0.94:
                reads.add(read.query_name)
    exclude(ARGS.pre, reads, ARGS.post)
    writeset(reads, ARGS.output_file)
    print(len(reads))

# ------------------------------------------------------------------------------
# 3. REMAPPING TO REFERENCE
# ------------------------------------------------------------------------------

elif ARGS.step == 'reference':
    ed_human = 10 if ARGS.max else 6
    reads = set()
    for input_file in ARGS.input_file.split(','):
        with pysam.AlignmentFile(input_file, 'rb', check_sq=False) as input_fo:
            for read in input_fo.fetch():
                number_mismatches = int(read.get_tag('NM'))
                read_length = int(read.infer_read_length())
                alignment_length = int(read.query_alignment_length)
                soft = read_length - alignment_length 
                number_mismatches += soft
                if number_mismatches <= ed_human:
                    reads.add(read.query_name)
    exclude(ARGS.pre, reads, ARGS.post)
    writeset(reads, ARGS.output_file)
    print(len(reads))

# ------------------------------------------------------------------------------
# 4. REPEAT PROFILING
# ------------------------------------------------------------------------------

elif ARGS.step == 'repeats':
    reads = set()
    with open(ARGS.input_file, 'r') as input_fo:
        for line in csv.reader(input_fo, delimiter='\t'):
            element, identity = line[0], float(line[2])
            alignment_len = float(line[3])
            start, end = int(line[6]), int(line[7])
            read_len = end - start + 1
            e = float(line[10])
            if e < 1e-05 and alignment_len >= 0.8 * read_len and \
                identity >= 0.9 * read_len:
                reads.add(element)
    exclude(ARGS.pre, reads, ARGS.post)
    writeset(reads, ARGS.output_file)
    print(len(reads))

# ------------------------------------------------------------------------------
# 5. CIRCULAR RNA PROFILING
# ------------------------------------------------------------------------------

elif ARGS.step == 'circrna':
    reads = set()
    for seq in SeqIO.parse(ARGS.input_file, 'fastq'):
        reads.add(seq.id)
    exclude(ARGS.pre, reads, ARGS.post)
    writeset(reads, ARGS.output_file)
    print(len(reads))

# ------------------------------------------------------------------------------
# 6. IMMUNE PROFILING
# ------------------------------------------------------------------------------

elif ARGS.step == 'immune':
    reads = set()
    for input_file in ARGS.input_file.split(','):
        with open(input_file, 'r') as input_fo:
            reader = csv.reader(input_fo)
            reader.next()
            for line in reader:
                reads.add(line[0])
    exclude(ARGS.pre, reads, ARGS.post)
    writeset(reads, ARGS.output_file)
    print(len(reads))

# ------------------------------------------------------------------------------
# 7. MICROBIOME
# ------------------------------------------------------------------------------

elif ARGS.step == 'microbiome':
    reads = set()
    for input_file in ARGS.input_file.split(','):
        with pysam.AlignmentFile(input_file, 'rb', check_sq=False) as input_fo:
            for read in input_fo.fetch():
                number_mismatches = int(read.get_tag('NM'))
                read_length = int(read.infer_read_length())
                alignment_length = int(read.query_alignment_length)
                identity = 1 - number_mismatches/float(alignment_length)
                if alignment_length >= 0.8 * read_length and identity >= 0.9:
                    reads.add(read.query_name)
    exclude(ARGS.pre, reads, ARGS.post)
    writeset(reads, ARGS.output_file)
    print(len(reads))
