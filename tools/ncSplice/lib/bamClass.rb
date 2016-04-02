#!/usr/bin/env ruby

class ReadBam

	# Holds methods for backsplice, intra- and inter-chromosomal junction detection. 
	# Usable for unstranded, paired-end sequencing data.
	# Takes sam/bam line to create object.
	#
	#
	# Create new SplicedRead object:
	#
	# id      - Read id of the read anchor.
	# strand  - Strand on which the anchor maps. 
	# chr     - Chromosome on which the anchor maps.
	# start   - Start position of anchor mapping.
	# cigar   - Cigar string of anchor mapping.
	attr_accessor :id, :strand, :chr, :start, :cigar
	
	def initialize(line)
		@id, @strand, @chr = line[0..2]
		@start = line[3].to_i
		@cigar = line.find {|l| l.match(/MD:Z:[[:digit:]]+/) }
	end
	
	def start
		@start.to_i
	end
	
	def strand
		@strand.to_i & 0x10 > 0 ? @strand = -1 : @strand = 1
	end
	
	
	# Examines whether read pairs fulfills circular RNA conditions.
	# Those are:
	# 1. Anchors have to be on the same strand.
	# 2. Anchors need to be on the same chromosome.
	# 3. Anchors need to match completely.
	# 4. Chromosome should not be on the "to exclude"-list.
	# 5. Anchors should map a max. of 100 kb apart from each other.
	# 6. Anchors have to be in head-to-tail configuration.
	#
	# anchor		- Other half of anchor pair, given as object created by SplicedRead class.
	# distance	- Integer of the max. distance anchors should be apart form each other.
	# exclude		- Array with chromosome names to ignore.
	#
	# Returns boolean.
	def circular?(anchor, distance, exclude)
		conditions = []
		
		conditions << (strand == anchor.strand) # 1.
		conditions << (chr == anchor.chr) # 2.
		conditions << [cigar, anchor.cigar].all? {|c| c == 'MD:Z:20'} # 3.
		conditions << exclude.all? {|c| ![chr, anchor.chr].include?(c)}	# 4.
		conditions << ((start - anchor.start).abs <= distance) # 5.
		if strand == 1 && anchor.strand == 1 # 6.
			conditions << (start > anchor.start)
		elsif strand == -1 && anchor.strand == -1
			conditions << (start < anchor.start)
		end
	
		conditions.all? {|a| a == true}
	end
	
	
	# Examines whether read pairs fulfills intra-chromosomal transcript conditions.
	# Those are:
	# 1. Anchors need to be on the same chromosome.
	# 2. Anchors need to match completely.
	# 3. Chromosome should not be on the "to exclude"-list.
	# 4. Anchors should map a min. of 1 mb apart from each other.
	#
	# anchor		- Other half of anchor pair given as object created by SplicedRead class.
	# distance	- Integer of the max. distance anchors should be apart form each other.
	# exclude		- Array with chromosome names to ignore
	#
	# Returns boolean.
	def intraChimeric?(anchor, distance, exclude)
		conditions = []
		
		conditions << (chr == anchor.chr) # 1.
		conditions << [cigar, anchor.cigar].all? {|c| c == 'MD:Z:20'} # 2.
		conditions << exclude.all? {|c| ![chr, anchor.chr].include?(c)}	# 3.
		conditions << ((start - anchor.start).abs >= distance) # 4.
	
		conditions.all? {|a| a == true }
	end


	# Examines whether read pairs fulfills intra-chromosomal transcript conditions.
	# Those are:
	# 1. Anchors need to be on different chromosomes.
	# 2. Anchors need to match completely.
	# 3. Chromosome should not be on the "to exclude"-list.
	#
	#	anchor	- Other half of anchor pair given as object created by SplicedRead class.
	# exclude	- Array with chromosome names to ignore.
	#
	# Returns boolean.
	def interChimeric?(anchor, exclude)
		conditions = []
		
		conditions << (chr != anchor.chr) # 1.
		conditions << [cigar, anchor.cigar].all? {|c| c == 'MD:Z:20'}	# 2.	
		conditions << exclude.all? {|c| ![chr, anchor.chr].include?(c)} # 3.

		conditions.all? {|a| a == true}
	end
end