#!/usr/bin/env ruby

module Alignment

	extend self
	
	# Trim end of alignment if mismatch is in the last 2 bp, likely to be a mapping error.
	#
	#
	# array	- Array of Integers from mapping with breakpoint_downstream or
	#					breakpoint_upstream methods.
	# mm		- Number of allowed mismatches.
	#
	# Return integer of trimmed alignment array.
	def trim(array, mm)
		array = array[0..-3] if array[-2..-1] == [0,1] && mm == 1
		array = array[0..-2] if array[-1] == 0 && mm == 1
		array.length
	end


	# Extend 3' anchor to upstream genomic region.
	# Use the trim-method to remove mismatches from the alignment end.
	#
	#
	# ref - Reference sequence to which DNA sequence will be compared.
	# mm  - Number of allowed mismatches.
	#
	# Return integer of alignment length.
	def upstream(seq, reference, mm)
		a = []
		mismatch = 0

		(seq.length - 1).downto(0).each do |i|
			if seq[i] == reference[i]
				a << 1
			else
				mismatch += 1
				return trim(a, mm) if mismatch > mm
				a << 0
			end
		end
	end


	# Extend 5' anchor to downstream genomic region.
	# Use the trim-method to remove mismatches from the alignment end.
	#
	#
	# ref - Reference sequence to which DNA sequence will be compared.
	# mm  - Number of allowed mismatches.
	#
	# Return integer of alignment length.
	def downstream(seq, reference, mm)
		a = []
		mismatch = 0
	
		0.upto(seq.length - 1).each do |i|
			if seq[i] == reference[i]
				a << 1
			else
				mismatch += 1
				return trim(a, mm) if mismatch > mm
				a << 0
			end
		end
	end


	# Create reverse complement of DNA.
	#
	#
	# dna - DNA sequence
	#
	# Return reverse complement.
	def reverse_complement(dna)
		complement = []
		dna.each_char do |s|
			if s == 'A'
				complement << 'T'
			elsif s == 'T'
				complement << 'A'
			elsif s == 'G'
				complement << 'C'
			elsif s == 'C'
				complement << 'G'
			elsif s == 'N'
				complement << 'N'
			end	
		end
		complement.join('').reverse
	end
	
	
	# Check read quality based on phred value.
	#
	#
	# quality	- Obserevd phred value of read.
	# phred		- Minimal phred value allowed.
	#
	# Return boolean.
	def quality_ok?(quality, phred)
  	quality.each_char.all? {|x| (x.ord - 33) >= phred}
	end


	# Compare number of mismatches reported by MD:Z tag to allowed number of mismatches.
	#
	#
	# mdz - String of MD:Z tag from bowtie2 alignment.
	# mm  - Integer of max. number of allowed mismatches.
	#
	# Return boolean.
	def mismatches?(mdz, mm)
		mdz = mdz.split(':').last
		counter = 0
		%w[A T G C].each {|x| counter += mdz.count(x) if !mdz.nil?}
		counter > mm
	end
	
	
	# Report genompic distance over which read mapped.
	#
	#
	# cigar				- Cigar string from bowtie2 alignment.
	# read_length	- Length of sequecning reas.
	#
	# Return false if insertions etc present, otherwise genomic length as interger.
	def genomic_mappinglength(cigar, read_length)
		if cigar == '{read_length}M'
			read_length
		elsif cigar.include?('N') && %w(I D S H P).all? {|x| !cigar.include?(x)}
			cigar = cigar.split(/N|M/).collect {|x| x.to_i}
			cigar.inject {|sum, x| sum + x}
		else
			false
		end
	end
	
	
	# Check whether spliced read has mapped mate in the right position and orientation.
	#
	#
	# mate1	- Mapping information for mate 1.
	# mate2	- Mapping information for mate 2.
	#
	# Return false if insertions etc present, otherwise genomic length as interger.
	def paired?(mate1, mate2)
		mate1_chr, mate1_start, mate1_stop, mate1_strand = mate1
		mate2_chr, mate2_start, mate2_stop, mate2_strand = mate2
		conditions = []
		
		conditions.push(mate1_chr == mate2_chr)
		conditions.push(mate1_start <= mate2_start && mate1_stop >= mate2_stop)
		conditions.push(mate1_strand != mate2_strand)
		conditions.all? { |con| con == true }
	end
end
