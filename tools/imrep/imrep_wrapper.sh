#!/bin/bash
# usage: imrep_wrapper.sh input_fasta

python `dirname $0`/imrep.py $1 output.txt --extendedOutput

echo "Compiling read names."
cat full_cdr3.txt partial_cdr3.txt | sed "s \t.*  " >immune_read_names.txt
