#!/usr/bin/env ruby

# 13/03/2016
# version: ruby 2.0.0
# 
# Wrapper script for circRNA detection.
#
#
# Options (update!):
# u - Fastq-file with unmapped reads from tophat in bam-format.
# p - Prefix used for all output files.
# x - Path to bowtie2-index directory and base name <index_directory>/<bt2-index>.
# a - Length of anchors, 20 bp is recommended.
# l - Length of read.
# f - Path to directory with one fasta files for each chromosome.
# s - File with chromosomes that should excluded.
#			One chromosome per line.
# singletones - File with paired reads for which only one mate mapped.
# sequencing-type - SE or PE for single-end or paired-end.
#
# Remarks:
# Implement additional option to delete intermediate files.


require 'optparse'
require 'open3'
require_relative "../lib/alignments.rb"
require_relative "../lib/analysis.rb"
require_relative "../lib/transSplicing.rb"

# options
##########################################################################################

options = {}
optparse = OptionParser.new do |opts|
  opts.banner = "ncSplice version 0.0.0 by Franziska Gruhl (franziska.gruhl@isb-sib.ch)\n
  Usage: ncSplice.rb -u <unmapped.fastq> -p <prefix> -x <index-directory>/<bt2-index> -a <anchor-length> -l <read-length> -c <chromsomes>/*.fa -s <exclude.txt> [options]"

  opts.on('-h', '--help', 'Display help screen.') do
    puts opts
    exit
  end
	
	opts.on('-v', '--version', 'Print ncSplice version and dependencies.') do
		puts '# ncSplice'
		puts 'ncSplice v0.1.0'
		puts "\n# Dependencies"
		puts ['ruby >= 2.0.0', 'samtools >= 1.0.0', 'bowtie2 >= 2.1.0'].join("\n")
    exit
  end

  options[:phred] = 25
  options[:anchor] = 20
	
 	opts.on('-u', '--unmapped <filename>', String, 'fastq file with unmapped reads') {|u| options[:unmapped] = u}
 	opts.on('-q', '--quality <integer>', Integer, 'Minimal phred quality unmapped reads need to have for further analysis.') {|q| options[:phred] = q}
  opts.on('-p', '--prefix <string>', String, 'Prefix for all files.') {|b| options[:prefix] = b}
  opts.on('-x', '--bowtie-index <directory>', String, 'Bowtie-index diretory and base name: <index-directory>/<bt2-index>.') {|x| options[:bowtie] = x}
  opts.on('-a', '--anchor-length <integer>', Integer, 'Length of the read anchor for remapping, default is 20 bp, shorter anchors will decrease the mapping precision and longer anchors will cause a reduction in candidates.') {|a| options[:anchor] = a}
  opts.on('-l', '--read-length <integer>', Integer, 'Length of the sequencing read.') {|l| options[:readlength] = l}  
  opts.on('-c', '--chromosome-files <directory>', 'Directory with chromosome fasta files, one fasta-file per chromosome.') {|c| options[:chromosomes] = c}
  opts.on('-s', '--skip-chr <filename>', String, 'Text file with chromosomes to exclude, such as the mitochondrial chromosome (recommended), chromosomes need to be listed in a separate text-file with one chromosome per line.') {|s| options[:skip] = s}
end

optparse.parse!(ARGV)
optparse.parse!(ARGV << '-h') if options.length <= 7

# local
unmapped_reads = options[:unmapped]
phred_quality = options[:phred]
prefix = options[:prefix]
bowtie_index = options[:bowtie]
anchor_length = options[:anchor]
read_length = options[:readlength]
options[:chromosomes][-1] == '/' ? chromosome_files = options[:chromosomes] : chromosome_files = "#{options[:chromosomes]}/"
skip = options[:skip]

# global
$logfile = File.open("#{prefix}_logfile.log", 'w')


# run
##########################################################################################

begin	
	Analysis.prepare_anchorpairs(unmapped_reads, anchor_length, "#{prefix}_anchors.fastq")
 	Analysis.bowtie_map(bowtie_index, "#{prefix}_anchors.fastq", "#{prefix}.bam")
 
 	anchor_pairs = Open3.popen3("samtools view #{prefix}.bam") do |stdin, stdout, stderr, t|
 		TransSplice.process_bam(stdout, chromosome_files, skip)
 	end

  TransSplice.seed_extension(anchor_pairs, anchor_length, read_length, chromosome_files, "#{prefix}_candidateReads.txt")	
  TransSplice.collaps_qnames("#{prefix}_candidateReads.txt", "#{prefix}_candidates.txt")
  TransSplice.candidates2fa("#{prefix}_candidates.txt", chromosome_files, read_length, "#{prefix}_faIndex.fa")
  Analysis.bowtie_build("#{prefix}_faIndex.fa", "#{prefix}")
  	
  Open3.popen3("mkdir #{prefix}_index") if !Dir.exists?("#{prefix}_index")
  Open3.popen3("mv #{prefix}_faIndex.fa #{prefix}_index/; mv *bt2 #{prefix}_index/")
  	
  Analysis.bowtie_map2("#{prefix}_index/#{prefix}", unmapped_reads, "#{prefix}_remapping.bam")
  
  Open3.popen3("samtools view #{prefix}_remapping.bam") do |stdin, stdout, stderr, t|
  	TransSplice.remapped_reads(stdout, "#{prefix}_remapped.txt", read_length)
  end
  	
  TransSplice.final_candidates("#{prefix}_candidates.txt", "#{prefix}_remapped.txt", "#{prefix}_final.txt")
	
rescue StandardError => err
	$logfile.puts "#{Time.new}: Error in ncSplice_transSplice.rb"
	$logfile.puts err.message
	err.backtrace.each {|line| $logfile.puts line}
	exit
end
