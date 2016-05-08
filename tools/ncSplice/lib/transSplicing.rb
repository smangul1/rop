#!/usr/bin/env ruby

require_relative "./alignments.rb"
require_relative "./bamClass.rb"

module TransSplice

	extend self

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
	def interChimeric?(anchor_left, anchor_right, exclude)
		conditions = []
		
		conditions << (anchor_left.chr != anchor_right.chr) # 1.
		conditions << [anchor_left.cigar, anchor_right.cigar].all? {|c| c == 'MD:Z:20'}	# 2.	
		conditions << exclude.all? {|c| ![anchor_left.chr, anchor_right.chr].include?(c)} # 3.

		conditions.all? {|a| a == true}
	end
	
	
	# Read bam-file from IO and process lines pair-wise to filter circRNA candidates.
	#
	#
	# input_file - Input file from IO
	# fasta      - Path to chromsomal fasta-files.
	# skip       - Text-file with chromsomes to skip, one/line.
	#
	# Return hash with valid anchor pairs.
	def process_bam(input_file, fasta, skip)

		# general settings
		exclude = []
		File.open(skip, 'r').readlines.each {|line| exclude << line.strip}
		firstline = TRUE 
		anchor_left = nil
		anchor_right = nil
		chr_a = nil
		chr_b = nil
		input_hash = {}

		# Initiate chromosome hash
		Dir.foreach(fasta) do |item|
			chr = item.sub('.fa', '')
			next if item == '.' || item == '..' || exclude.include?(chr)  
			input_hash[chr] = {}
		end

		input_hash.each_key do |chr_a|
			input_hash.keys.each {|chr_b| input_hash[chr_a][chr_b] = []}
		end

		# read bam file
		input_file.each do |line|
			line = line.strip.split(/\s+/)
		
			if firstline 
				anchor_left = ReadBam.new(line)
				firstline = FALSE
				chr_a = anchor_left.chr
			else
				anchor_right = ReadBam.new(line)
				chr_b = anchor_right.chr
				
				if input_hash.has_key?(chr_a) && interChimeric?(anchor_left, anchor_right, exclude)
					
					if anchor_left.strand == 1 && anchor_right.strand == 1
						input_hash[chr_b][chr_a] << [anchor_right, anchor_left] 
					elsif anchor_left.strand == -1 && anchor_right.strand == -1
						input_hash[chr_a][chr_b]  << [anchor_left, anchor_right] 
					else
						input_hash[chr_b][chr_a]  << [anchor_right, anchor_left] 
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
	
		input_hash.each do |chr_a, chromosomes|
			# Load reference
			fasta_file = File.open("#{fasta}#{chr_a}.fa", 'r')
			header = fasta_file.gets.strip
			dna_a = fasta_file.read.gsub(/\n/, '')

			chromosomes.each do |chr_b, anchorpairs|
			  fasta_file = File.open("#{fasta}#{chr_b}.fa", 'r')
  			header = fasta_file.gets.strip
  			dna_b = fasta_file.read.gsub(/\n/, '')

				# Loop through hash to extend seeds for each pair
				anchorpairs.each do |pair|
					upstream, downstream = pair
					qname, mate, read = upstream.id.split('_')[0..2]

					upstream.strand == 1 ? upstream_read = read : upstream_read = Alignment.reverse_complement(read)
					downstream.strand == 1 ? downstream_read = read : downstream_read = Alignment.reverse_complement(read)
					
					up = dna_a[upstream.start - read_length + anchor_length..upstream.start + anchor_length - 1].upcase
					down = dna_b[downstream.start..downstream.start + read_length - 1].upcase	
				
					if upstream.strand == downstream.strand
						upstream_alignmentlength = Alignment.upstream(upstream_read, up, mm)
						downstream_alignmentlength = Alignment.downstream(downstream_read, down, mm)
						upstream_breakpoint = upstream.start - upstream_alignmentlength + anchor_length	
						downstream_breakpoint = downstream.start + downstream_alignmentlength - 1

					elsif upstream.strand == 1 && downstream.strand == -1
						down = dna_b[downstream.start - read_length + anchor_length..downstream.start + anchor_length - 1].upcase
						upstream_alignmentlength = Alignment.upstream(upstream_read, up, mm)
						downstream_alignmentlength = Alignment.upstream(downstream_read, down, mm)
						upstream_breakpoint = upstream.start - upstream_alignmentlength + anchor_length	
						downstream_breakpoint = downstream.start - downstream_alignmentlength + anchor_length	
				
					else
						up = dna_a[upstream.start..upstream.start + read_length - 1].upcase	
						upstream_alignmentlength = Alignment.downstream(upstream_read, up, mm)
						downstream_alignmentlength = Alignment.downstream(downstream_read, down, mm)
						upstream_breakpoint = upstream.start + upstream_alignmentlength - 1
						downstream_breakpoint = downstream.start + downstream_alignmentlength - 1
					end

					total_alignmentlength = upstream_alignmentlength + downstream_alignmentlength

					if total_alignmentlength >= read_length && total_alignmentlength <= max_overhang
						overhang = total_alignmentlength - read_length
	
						qname = qname.to_sym
						summary = [chr_a, upstream_breakpoint, upstream.strand, chr_b, downstream_breakpoint, downstream.strand, total_alignmentlength, mate] 
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
		end

		File.open(output_file, 'w') do |output|
			output_hash.each do |qname, v| 
				output.puts ["#{qname.to_s}/#{v[-1]}", v[0..-2]].join("\t") if (v[2] - v[1]).abs >= read_length
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
			pos_a = line[1..3].join(':')
			pos_b = line[4..6].join(':')
			pos = [pos_a, pos_b].join(':')

			alignment_length = line[-1]
	
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
			chr_a, pos_a, strand_a, chr_b, pos_b, strand_b = line[0..5]
			pos = [chr_a, pos_a, chr_b, pos_b].join(':')
	
			chromosomes[chr_a] = {} if !chromosomes.has_key?(chr_a)
			
			if !chromosomes.has_key?(chr_b)
				chromosomes[chr_a][chr_b] = [line]
		
			# 2nd elsif to exclude reads that map on same junction but opposite ends		
			elsif chromosomes[chr_a].has_key?(chr_b) && !positions.include?(pos)
				chromosomes[chr_a][chr_b].push(line)
				positions << pos
			end
		end

		# Output
		output = File.open(output_file, 'w') do |output|
			chromosomes.each do |chr_a, values|
				fasta_file = File.open("#{fasta}#{chr_a}.fa", 'r')
				header = fasta_file.gets.strip
				dna_a = fasta_file.read.gsub(/\n/, '')
				
				values.each do |chr_b, values|
				  fasta_file = File.open("#{fasta}#{chr_b}.fa", 'r')
					header = fasta_file.gets.strip
					dna_b = fasta_file.read.gsub(/\n/, '')

					values.each do |v|
						bp_a, bp_b = v[1].to_i, v[4].to_i
						strand_a, strand_b = v[2], v[5]
						overlap = v[-1].to_i - read_length
						l = read_length - exoncov 
										
						upstream = dna_a[bp_a..bp_a + overlap + l].upcase	
						downstream = dna_b[bp_b - l - overlap + 1..bp_b - overlap].upcase
					
						if strand_a == '1' && strand_b == '-1'
							downstream = Alignment.reverse_complement(dna_b[bp_b..bp_b + l].upcase)
						elsif strand_a == '-1' && strand_b == '1'
							upstream = Alignment.reverse_complement(dna_a[bp_a - l + 1..bp_a].upcase)
						end
		
						id = [v[0..1], v[3..4]].join(':')
						output.puts [">#{id}", downstream + upstream].join("\n")
					end
				end
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Wrote loci to fasta-file."
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
			mdz = line.match(/MD:Z:\S*/).to_s
			line = line.strip.split(/\s+/)
			qname, mate = line[0].split('/')
			pos = line[2].split(':')
			cigar = line[5]
	
			if !remapped.has_key?(qname) && Alignment.max_mismatches?(mdz, mm) && cigar == "#{read_length}M"
				remapped[qname] = [pos, mate]
			else	
				remapped.delete(qname)
			end
		end

		# Output
		File.open(output_file, 'w') do |output|
			remapped.each {|k, v| output.puts ["#{k}/#{v[-1]}", v[0]].join("\t")}
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
		candidates = {}
		all_ids = {}

		# Read circular candidates into hash
		File.open(before, 'r').readlines.each do |line|
			line = line.strip.split("\t")
	
			pos = [line[0..1], line[3..4]].join(':')
			read_count = line[6].to_i
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

			candidates[pos] = {:counts => read_count, :qnames => qname}
		end

		# Read remapped readpairs and compare them to initial candidates
		File.open(after, 'r').readlines.each do |line|
			line = line.strip.split("\t")
	
			qname = line[0]
			pos = [line[1..2], line[3..4]].join(':')
			k1, k2 = qname.split(':')[3..4]

			read_unused = (!all_ids.has_key?(k1) || !all_ids[k1].has_key?(k2) || !all_ids[k1][k2].include?(qname)) 
						
			# Add read if read is not already used (condition 2)
			if candidates.has_key?(pos) && read_unused 
				candidates[pos][:counts] += 1
				candidates[pos][:qnames] << qname
			end
		end

		# Output
		File.open(output_file, 'w') do |output|
			output.puts %w(chr pos_a chr_b pos_b readCounts qnames).join("\t")
	
			candidates.each do |pos, v| 
				output.puts [pos.split(':'), v[:counts], v[:qnames].join(';')].join("\t")
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Final candidate list finished."
	end
end
