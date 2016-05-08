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
require_relative "../lib/bamClass.rb"

# options
##########################################################################################

options = {}
optparse = OptionParser.new do |opts|
  opts.banner = "ncSplice version 0.0.0 by Franziska Gruhl (franziska.gruhl@isb-sib.ch)\n
  Formatting of input file\n
  Usage: readPreperation.rb -p <prefix> -i <input> [options]"

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
	
	opts.on('-u', '--input <filename>', String, 'Input file with unmapped reads') {|u| options[:unmapped] = u}
	opts.on('-f', '--format <string>', String, 'Input format of unmapped reads, can be \'fastq\', \'bam\' or \'fasta\'') {|f| options[:format] = f}
	opts.on('-q', '--quality <integer>', Integer, 'Minimal phred quality unmapped reads need to have for further analysis') {|q| options[:phred] = q}
  opts.on('-s', '--sequencing-type <string>', String, 'Sequencing type, \'se\' for single-end and \'pe\' for paired-end sequencing, default is \'se\'') {|seq| options[:sequencing] = seq}
  opts.on('-p', '--prefix <string>', String, 'Prefix for all files.') {|p| options[:prefix] = p}
  opts.on('-a', '--anchor-length <integer>', Integer, 'Length of the read anchor for remapping, default is 20 bp, shorter anchors will decrease the mapping precision and longer anchors will cause a reduction in candidates.') {|a| options[:anchor] = a}
  opts.on('-x', '--bowtie-index <directory>', String, 'Bowtie-index diretory and base name: <index-directory>/<bt2-index>.') {|x| options[:bowtie] = x}
  opts.on('-m', '--mapped-reads <filename>', String, 'Mapped reads for filtering of singletons') {|m| options[:mapped] = m}
end

optparse.parse!(ARGV)
optparse.parse!(ARGV << '-h') if options.length < 7

# local
unmapped_reads = options[:unmapped]
input_format = options[:format]
phred_quality = options[:phred]
sequencing_type = options[:sequencing]
prefix = options[:prefix]
anchor_length = options[:anchor]
bowtie_index = options[:bowtie]
mapped_reads = options[:mapped]

# global
$logfile = File.open("#{prefix}_logfile.log", 'w')


# run
##########################################################################################

begin

	# if fasta, prepare reads directly in fasta-format
	if input_format == 'fasta'
			Analysis.fasta2anchors(unmapped_reads, anchor_length, sequencing_type, "#{prefix}_anchors.fasta")
			Analysis.bowtie_map(bowtie_index, "#{prefix}_anchors.fasta", input_format, "#{prefix}.bam")	
	
	else
		# if bam, transform unmapped reads into fastq-format before preparing anchors
		if input_format == 'bam'
			Open3.popen3("samtools view #{unmapped_reads}") do |stdin, stdout, stderr, t|
				Analysis.bam2fastq(stdout, "#{prefix}_unmapped.fastq", phred_quality)
			end
			Analysis.prepare_anchorpairs("#{prefix}_unmapped.fastq", anchor_length, sequencing_type, "#{prefix}_anchors.fastq")	
		
		# if fastq, prepare anchors from fastq-format
		else
			Analysis.prepare_anchorpairs(unmapped_reads, anchor_length, sequencing_type, "#{prefix}_anchors.fastq")
		end	

		# map in fastq format
		Analysis.bowtie_map1(bowtie_index, "#{prefix}_anchors.fastq", input_format, "#{prefix}.bam")	
	end
		
	if sequencing_type == 'pe'
		Open3.popen3("samtools view -bF 1 #{mapped_reads} > singletons.bam")
	end
	
rescue StandardError => err
	$logfile.puts "#{Time.new}: Error in ncSplice_prepareUnmapped.rb"
	$logfile.puts err.message
	err.backtrace.each {|line| $logfile.puts line}
	exit
end
