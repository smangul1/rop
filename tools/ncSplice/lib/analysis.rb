#!/usr/bin/env ruby

require_relative "alignments.rb"
require_relative "bamClass.rb"

module Analysis

	extend self
	
	# Catch exit status of system process and report it to logfile.
	# Exit program if subprocess was not finished successfully.
	#
	#
	# t       - Open3 object for thread.
	# stderr  - Open3 object for STDERR.
	# name    - Name for logfile.
	# logfile - Name of logfile to write to.
	#
	# Return string that is written to logfile. 
	def system_exitcode(t, stderr, name)	
		if t.value.success?
			$logfile.puts "#{Time.new.strftime("%c")}: #{name} succeded."
			if stderr.any?
				$logfile.puts "#{name} output:"
				stderr.readlines.each {|line| $logfile.puts line}
			end
		else
			$logfile.puts "#{Time.new.strftime("%c")}: Error in #{name}:"
			stderr.readlines.each {|line| $logfile.puts line}
			exit
		end
	end
 
 
  # Convert unmapped.bam to fastq-format
 	#
 	#
 	# input_file		- Unmapped reads from tophat in bam format.
 	# output_file		-	Name of output file (in fastq).
 	# phred_quality	- Minimal phredquality for read to be kept.
 	#
 	# Unmapped reads in fastq format.
 	def bam2fastq(input_file, output_file, phred_quality)
 		File.open(output_file, 'w') do |output|
			input_file.each do |line|
  			line = line.strip.split(/\s+/)
  
  			flag = line[1].to_i
  			flag & 0x40 > 0 ? mate = '1' : mate = '2'
  			
  			qname, sequence, quality = line[0], line[9], line[10] 
  			output.puts "@#{qname}/#{mate}", sequence, '+', quality if quality_ok?(quality, phred_quality)
  		end
  	end
  	$logfile.puts "#{Time.new.strftime("%c")}: Converted unmapped.bam into fastq-format."	
	end
	
	
	# Prepare fastq-file with anchors from unmapped.fastq.
	#
	#
	# input_file    - unmapped.fastq.
	# anchor_length - Length of anchor, default is 20 bp.
	# output_file   - Name of output file.
	#
	# Return fastq-file with anchor pairs.
	def prepare_anchorpairs(input_file, anchor_length, output_file)	
		name, mate, seq, quality = nil
		counter = -1

		File.open(output_file, 'w') do |output| 
			File.open(input_file, 'r').readlines.each do |line|
				counter += 1
				line = line.strip
			
				if counter % 4 == 0 
					if line.match('/')
						name, mate = line.split('/')
					else
						name, mate = line, '1'
					end
				elsif counter % 4 == 1
					seq = line
				
				elsif counter % 4 == 3
					quality = line
			
					name_A = "#{name}_#{mate}_#{seq}_A"
					name_B = "#{name}_#{mate}_#{seq}_B"
			
					seq_A = seq[0..anchor_length - 1]
					seq_B = seq[-anchor_length..-1]
	
					quality_A = quality[0..anchor_length - 1]
					quality_B = quality[-anchor_length..-1]
			
					output.puts [name_A, seq_A, '+', quality_A, name_B, seq_B, '+', quality_B].join("\n")
				
					name, mate, seq, quality = nil
					counter = -1
				end 
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Anchor preparation succeded."	
	end


	# Mapping of anchors via bowtie2 system call.
	#
	#
	# bowtie_index  - Path to bowtie index and index base name
	# fastq_file    - Fastq-file with anchors
	# output_file   - Name of output file
	# logfile       - Name of logfile to write to.
	#
	# Return mapped anchors in bam-format.
	def bowtie_map(bowtie_index, fastq_file, output_file)
		stdin, stdout, stderr, t = Open3.popen3("bowtie2 -x #{bowtie_index} -q -U #{fastq_file} | samtools view -bS - > #{output_file}")
		system_exitcode(t, stderr, 'Bowtie2')
	end


	# Read bam-file from IO and process lines pair-wise to filter circRNA candidates.
	#
	#
	# input_file - Input file from IO
	# fasta      - Path to chromsomal fasta-files.
	# skip       - Text-file with chromsomes to skip, one/line.
	#
	# Return hash with valid anchor pairs.
	def process_bam(input_file, fasta, skip, distance=100000)

		# general settings
		exclude = []
		File.open(skip, 'r').readlines.each {|line| exclude << line.strip}
		firstline = TRUE 
		anchor_left = nil
		anchor_right = nil
		input_hash = {}

		# Initiate chromosome hash
		Dir.foreach(fasta) do |item|
			chr = item.sub('.fa', '')
			next if item == '.' || item == '..' || exclude.include?(chr)  
			input_hash[chr] = []
		end

		# read bam file
		input_file.each do |line|
			line = line.strip.split(/\s+/)
		
			if firstline 
				anchor_left = ReadBam.new(line)
				firstline = FALSE
			else
				anchor_right = ReadBam.new(line)
				chr = anchor_right.chr

				if anchor_left.circular?(anchor_right, distance, exclude)
		
					# store coordinate sorted
					if anchor_left.strand == 1
						input_hash[chr] << [anchor_right, anchor_left] 
					else
						input_hash[chr] << [anchor_left, anchor_right]
					end
				end
	
				anchor_left, anchor_right = nil
				firstline = TRUE
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Found anchor pairs."		
		input_hash
	end


	# Seed extension for anchor pairs. They are kept as candidate if extension is succesful.
	#
	#
	# input_hash    - Input hash with anchor pairs.
	# anchor_length - Length of anchors
	# read_length   - Length of sequencing read.
	# fasta         - Path to chromsomal fasta-files.
	#	output_file   - Name of output_file.
	#
	# Return tab-delimited file with initial candidates.
	def seed_extension(input_hash, anchor_length, read_length, fasta, output_file, mm = 1, max_overhang = read_length + 8)

		output_hash = {}
	
		input_hash.each do |chr, anchorpairs|
			chr = chr.to_s
	
			# Load reference
			fasta_file = File.open("#{fasta}#{chr}.fa", 'r')
			header = fasta_file.gets.strip
			dna = fasta_file.read.gsub(/\n/, '')
	
			# Loop through hash to extend seeds for each pair
			anchorpairs.each do |pair|
				upstream, downstream = pair
				qname, mate, read = upstream.id.split('_')[0..2]

				strand, chr, upstream_start, downstream_start = upstream.strand, upstream.chr, upstream.start - 1, downstream.start - 1

				read = Alignment.reverse_complement(read) if strand == -1
				upstream_dna = dna[upstream_start - read_length + anchor_length..upstream_start + anchor_length - 1].upcase
				downstream_dna = dna[downstream_start..downstream_start + read_length - 1].upcase				
				upstream_alignmentlength = Alignment.upstream(read, upstream_dna, mm)
				downstream_alignmentlength = Alignment.downstream(read, downstream_dna, mm)
				total_alignmentlength = upstream_alignmentlength + downstream_alignmentlength

				if total_alignmentlength >= read_length && total_alignmentlength <= max_overhang
					upstream_breakpoint = upstream_start - upstream_alignmentlength + anchor_length	
					downstream_breakpoint = downstream_start + downstream_alignmentlength - 1
					overhang = total_alignmentlength - read_length
	
					qname = qname.to_sym
					summary = [chr, upstream_breakpoint, downstream_breakpoint, strand, total_alignmentlength, mate] 
					
					# Candidates for which both, R1 and R2, are present are deleted
					# One read can neither fall on two different non-canonical nor the same junction
					if !output_hash.has_key?(qname)
						output_hash[qname] = summary
					else
						output_hash.delete(qname)
					end
				end
			end
		end
		
		File.open(output_file, 'w') do |output|
			output_hash.each do |qname, v| 
				output.puts ["#{qname.to_s}/#{v[-1]}", v[0..-2]].join("\t") if v[2] - v[1] >= read_length
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Seed extension succeded."
	end
	
	
	# Collaps candidates reads onto common loci.
	#
	#
	# input_file	- Input file.
	# output_file	- Name of output_file.
	#
	# Return tab-delimited file with candidate loci.
	def collaps_qnames(input_file, output_file)
	
		loci = {}
	
		# Read candidate loci and count reads/locus
		File.open(input_file, 'r').readlines.each do |line|
			line = line.strip.split("\t")
			qname = line[0]
			base = qname.gsub(/\/[1,2]/, '')
			pos = line[1..3].join(':')
			alignment_length = line[5]
	
			if !loci.has_key?(pos)
				loci[pos] = {:count => 1, :qnames => [qname], :l => alignment_length}
			else 
				loci[pos][:qnames] << qname
				loci[pos][:count] += 1
			end
		end

		# Output
		File.open(output_file, 'w') do |output|
			loci.each do |pos, v| 
				output.puts [pos.split(':'), v[:count], v[:l], v[:qnames].join(';')].join("\t") if v[:count] > 0
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Collapsed anchor pairs to single loci."
	end
	
	
	# Transform candidate loci into fasta-format.
	#
	#
	# input_file  - Input file.
	# fasta       - Path to chromsomal fasta-files.
	# read_length - Length of sequencing read.
	# output_file	- Name of output_file.
	#
	# Return tab-delimited file with candidate loci.
	def candidates2fa(input_file, fasta, read_length, output_file, exoncov=8)
		chromosomes = {}
		positions = []
		
		# Input into hash sorted by chromosomes
		File.open(input_file, 'r').readlines.each do |line|
			line = line.strip.split("\t")[0..-2]
			pos = line[0..2].join(':')
			chr = line[0]
	
			if !chromosomes.has_key?(chr)
				chromosomes[chr] = [line]
		
			# 2nd elsif to exclude reads that map on same junction but opposite ends		
			elsif chromosomes.has_key?(chr) && !positions.include?(pos)
				chromosomes[chr].push(line)
				positions << pos
			end
		end

		# Output
		output = File.open(output_file, 'w') do |output|
			chromosomes.each do |chr, values|
				fasta_file = File.open("#{fasta}#{chr}.fa", 'r')
				header = fasta_file.gets.strip
				dna = fasta_file.read.gsub(/\n/, '')
		
				values.each do |v|
					bp_a, bp_b = v[1..2].collect {|x| x.to_i}
					overlap = v[-1].to_i - read_length
					l = read_length - exoncov 
			
					upstream = dna[bp_a..bp_a + overlap + l - 1].upcase	
					downstream = dna[bp_b - l - overlap + 1..bp_b - overlap].upcase
			
					output.puts [">#{v[0..2].join(':')}", downstream + upstream].join("\n")
				end
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Wrote loci to fasta-file."
	end
	
	
	# Built bowtie-index from candidate loci.
	#
	#
	# input_file - Input file.
	# logfile    - Name of logfile to write to.
	#
	# Return tab-delimited file with candidate loci.
	def bowtie_build(input_file, prefix)
		stdin, stdout, stderr, t = Open3.popen3("bowtie2-build -q -f #{input_file} #{prefix}")
	system_exitcode(t, stderr, 'Building bowtie2 index')
	end
	
	
	# Filter remapped reads, keep those that haven't been used yet and fullfill criteria
	#
	#
	# input_file  - Input file.
	# output_file - Name of output_file.
	# read_length - Length of sequencing read.
	# mm          - Number of mismatchs, default 2.	
	#
	# Return tab-delimited file with candidate loci.
	def remapped_reads(input_file, output_file, read_length, mm=2)
		remapped = {}
		
		# Filter remapped reads
		input_file.each do |line|
			mdz = line.match(/MD:Z:[[:digit:]]+/).to_s
			line = line.strip.split(/\s+/)
			qname, mate = line[0].split('/')
			chr, start, stop = line[2].split(':')
			cigar = line[5]
	
			if !remapped.has_key?(qname) && !Alignment.mismatches?(mdz, mm) && cigar == "#{read_length}M"
				remapped[qname] = [chr, start, stop, mate]
			else	
				remapped.delete(qname)
			end
		end

		# Output
		File.open(output_file, 'w') do |output|
			remapped.each {|k, v| output.puts ["#{k}/#{v[-1]}", v[0..2]].join("\t")}
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Found remapped reads."
	end
	
	
	# Compare first mapping and remapping to get final candidates.
	#
	#
	# before      - Input file.
	# after       - Name of logfile to write to.
	# output_file - Name of output_file.
	#
	# Return tab-delimited file with candidate loci.
	def final_candidates(before, after, output_file)
		circles = {}
		all_ids = {}

		# Read circular candidates into hash
		File.open(before, 'r').readlines.each do |line|
			line = line.strip.split("\t")
	
			pos = line[0..2].join(':')
			read_count = line[3].to_i
			qname = line[-1].split(';')
	
			# Create qname index to make search faster
			# Remark 2
			qname.each do |q|
				k1, k2 = q.split(':')[3..4]
		
				all_ids[k1] = {} if !all_ids.has_key?(k1)
		
				if !all_ids[k1].has_key?(k2)
					all_ids[k1][k2] = [q]
				else
					all_ids[k1][k2] << q
				end
			end

			circles[pos] = {:counts => read_count, :qnames => qname}
		end

		# Read remapped readpairs and compare them to initial candidates
		File.open(after, 'r').readlines.each do |line|
			line = line.strip.split("\t")
	
			qname = line[0]
			pos = line[1..3].join(':')
			k1, k2 = qname.split(':')[3..4]

			read_unused = (!all_ids.has_key?(k1) || !all_ids[k1].has_key?(k2) || !all_ids[k1][k2].include?(qname)) 
						
			# Add read if read is not already used (condition 2)
			if circles.has_key?(pos) && read_unused 
				circles[pos][:counts] += 1
				circles[pos][:qnames] << qname
			end
		end

		# Output
		File.open(output_file, 'w') do |output|
			output.puts %w(chr pos_a pos_b readCounts qnames).join("\t")
	
			circles.each do |pos, v| 
				output.puts [pos.split(':'), v[:counts], v[:qnames].join(';')].join("\t")
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Final candidate list finished."
	end
end