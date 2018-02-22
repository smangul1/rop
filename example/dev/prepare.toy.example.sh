#obtained HMP sample prepared by https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2
wget https://bitbucket.org/biobakery/biobakery/raw/tip/demos/biobakery_demos/data/metaphlan2/input/SRS014476-Supragingival_plaque.fasta.gz


wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/fasta-to-fastq/fasta_to_fastq.pl

perl fasta_to_fastq.pl SRS014476-Supragingival_plaque.fasta >SRS014476-Supragingival_plaque.fastq

