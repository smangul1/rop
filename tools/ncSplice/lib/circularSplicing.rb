#!/usr/bin/env ruby

require_relative "./bamClass.rb"
require_relative "./alignments.rb"

module CircRNA

	extend self
	
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
	def circular?(anchor_left, anchor_right, distance, exclude)
		conditions = []
		
		conditions << (anchor_left.strand == anchor_right.strand) # 1.
		conditions << (anchor_left.chr == anchor_right.chr) # 2.
		conditions << [anchor_left.cigar, anchor_right.cigar].all? {|c| c == 'MD:Z:20'} # 3.
		conditions << exclude.all? {|c| ![anchor_left.chr, anchor_right.chr].include?(c)}	# 4.
		conditions << ((anchor_left.start - anchor_right.start).abs <= distance) # 5.
		if anchor_left.strand == 1 && anchor_right.strand == 1 # 6.
			conditions << (anchor_left.start > anchor_right.start)
		elsif anchor_left.strand == -1 && anchor_right.strand == -1
			conditions << (anchor_left.start < anchor_right.start)
		end
	
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

				if circular?(anchor_left, anchor_right, distance, exclude)
		
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
	def seed_extension(input_hash, anchor_length, read_length, fasta, sequencing_type, output_file, mm = 1, max_overhang = read_length + 8)

		candidates = {}
	
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
				read = Alignment.reverse_complement(read) if upstream.strand == -1

				upstream_dna = dna[upstream.start - read_length + anchor_length..upstream.start + anchor_length - 1].upcase
				downstream_dna = dna[downstream.start..downstream.start + read_length - 1].upcase	
							
				upstream_alignmentlength = Alignment.upstream(read, upstream_dna, mm)
				downstream_alignmentlength = Alignment.downstream(read, downstream_dna, mm)
				total_alignmentlength = upstream_alignmentlength + downstream_alignmentlength
				
				upstream_breakpoint = upstream.start - upstream_alignmentlength + anchor_length	
				downstream_breakpoint = downstream.start + downstream_alignmentlength - 1

				mapping_conditions = [total_alignmentlength >= read_length, total_alignmentlength <= max_overhang, (downstream_breakpoint - upstream_breakpoint).abs >= read_length]
				
				if mapping_conditions.all? {|condition| condition == true}
					overhang = total_alignmentlength - read_length
					summary = [upstream.chr, upstream_breakpoint, downstream_breakpoint, upstream.strand, total_alignmentlength, 'se', mate] 
					
					if !candidates.has_key?(qname)
						candidates[qname] = summary 
					
					# R1 and R2 are deleted if both are present, but fall on different junctions
					else
						same_pos = candidates[qname][0..2] == summary[0..2]
						diff_strand = candidates[qname][3] != summary[3]
						diff_mate = candidates[qname][-1] != summary[-1]
						#$stdout.puts [same_pos, diff_mate].join("\t")
						if same_pos && diff_strand && diff_mate
							candidates[qname][-2] = 'pe'
						else 
							candidates.delete(qname)
						end
					end
				end
			end
		end
		
		File.open(output_file, 'w') do |output|
			candidates.each do |qname, v| 
				output.puts ["#{qname}/#{v[-1]}", v[0..-2]].join("\t") 
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
			support = line[-1]
			
			if !loci.has_key?(pos)
				loci[pos] = {:count => 1, :qnames => [qname], :l => alignment_length, :type => support}
			else 
				loci[pos][:qnames] << qname
				loci[pos][:count] += 1
				loci[pos][:type] = support if support == 'pe'
			end
		end

		# Output
		File.open(output_file, 'w') do |output|
			loci.each do |pos, v| 
				output.puts [pos.split(':'), v[:count], v[:l], v[:type], v[:qnames].join(';')].join("\t") if v[:count] > 0
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
			line = line.strip.split("\t")[0..-3]
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
			
					upstream = dna[bp_a..bp_a + overlap + l].upcase	
					downstream = dna[bp_b - l - overlap + 1..bp_b - overlap].upcase

					output.puts [">#{v[0..2].join(':')}", downstream + upstream].join("\n")
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
	def remapped_reads(input_file, output_file, read_length, sequencing_type, mm=2)
		remapped = {}

		# Filter remapped reads
		input_file.each do |line|
			mdz = line.match(/MD:Z:\S*/).to_s			
			line = line.strip.split(/\s+/)
			qname, mate = line[0].split('/')[0..1].join('/'), line[0].split('/')[-1]
			flag = line[1].to_i
			
			flag & 0x10 > 0 ? strand = -1 : strand = 1
			pos = line[2].split(':')
			cigar = line[5]
			
			if Alignment.max_mismatches?(mdz, mm) && cigar == "#{read_length}M"
				if !remapped.has_key?(qname)
					remapped[qname] = [pos, strand, 'se', mate]
				else
					
					# if second mate is unmapped, then chr/start/stop need to be the same and strand/mate need to be opposite
					if remapped[qname][0] == pos && remapped[qname][1] != strand && remapped[qname][-1] != mate
						remapped[qname][2] = 'pe'
					else
						remapped.delete(qname)
					end
				end
			end
		end

		if sequencing_type == 'pe'
			Open3.popen3("samtools view singletons.bam") do |stdin, stdout, stderr, t|
				stdout.each do |line|
  		
  				line = line.strip.split(/\s+/)
  				qname, flag, chr, start = line[0..3]  
  			
  				if remapped.has_key?(qname)
  				
  					flag.to_i & 0x10 > 0 ? strand = -1 : strand = 1
  					cigar = line[5]
						distance = Alignment.genomic_mappinglength(cigar, read_length)
			
						if distance != false
							strand == 1 ? stop = start.to_i + distance : stop = start.to_i - distance
							mate_1 = [remapped[qname][0], remapped[qname][1]].flatten
							mate_2 = [chr, start.to_i, stop, strand]
							remapped[qname][2] = 'pe' if Alignment.paired?(mate_1, mate_2)
						end
					end
				end
			end
		end
		# Output
		File.open(output_file, 'w') do |output|
			remapped.each {|k, v| output.puts ["#{k}/#{v[-1]}", v[0..-2]].join("\t")}
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
	def final_candidates(input_file, output_file)

		loci = {}
		
		File.open(input_file, 'r').readlines.each do |line|
			line = line.strip.split("\t")
			qname = line[0]
			base = qname.gsub(/\/[1,2]/, '')
			pos = line[1..3].join(':')
			support = line[-1]
			qname = ["#{base}/1", "#{base}/2"] if support == 'pe'

			if !loci.has_key?(pos)
				loci[pos] = {:count => 1, :qnames => [qname], :type => support}
			else 
				loci[pos][:qnames] << qname
				loci[pos][:count] += 1
				loci[pos][:type] = support if support == 'pe'
			end
		end
		
		# Output
		File.open(output_file, 'w') do |output|
			loci.each do |pos, v| 
				output.puts [pos.split(':'), v[:count], v[:type], v[:qnames].join(';')].join("\t")
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Final candidate list finished."
	end
end