#!/usr/bin/env ruby

require_relative "alignments.rb"

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
 
 
  # Convert bam- into fastq-format
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
  			output.puts "@#{qname}/#{mate}", sequence, '+', quality if Alignment.quality_ok?(quality, phred_quality)
  		end
  	end
  	$logfile.puts "#{Time.new.strftime("%c")}: Converted unmapped.bam into fastq-format."	
	end
	
	
	# Prepare fasta-file with anchors from unmapped reads in fasta-format.
	#
	#
	# input_file    	- Unmapped reads in fasta-format.
	# anchor_length 	- Length of anchor, default is 20 bp.
	# sequencing_type	- Sequencing type, either se for single-end, or pe for paired-end.
	# output_file   	- Name of output file.
	#
	# Return fastq-file with anchor pairs.
	def fasta2anchors(input_file, anchor_length, sequencing_type, output_file)
		counter = -1
		name, mate, seq = nil, nil, nil
		
		File.open(output_file, 'w') do |output|	
			File.open(input_file, 'r').each do |line|
				counter += 1
				
				if counter % 2 == 0
					if sequencing_type == 'se'
						name, mate = line.strip.match(/(?<=\>)(\S*)/).to_s, 1
					else
						name, mate = line.strip.match(/(?<=\>)(\S*)/).to_s.split('/')
					end
				elsif counter % 2 == 1
					seq = line.strip	

					name_A = "#{name}_#{mate}_#{seq}_A"
					name_B = "#{name}_#{mate}_#{seq}_B"
					
					seq_A = seq[0..anchor_length - 1]
					seq_B = seq[-anchor_length..-1]

					output.puts [">#{name_A}", seq_A, ">#{name_B}", seq_B].join("\n")
					name, mate, seq = nil, nil, nil
					counter = -1
				end
			end
		end
	end
	
	
	# Prepare fastq-file with anchors from unmapped reads in fastq-format.
	#
	#
	# input_file    	- Unmapped reads in fastq-format.
	# anchor_length 	- Length of anchor, default is 20 bp.
	# sequencing_type	- Sequencing type, either se for single-end, or pe for paired-end.
	# output_file   	- Name of output file.
	#
	# Return fastq-file with anchor pairs.
	def prepare_anchorpairs(input_file, anchor_length, sequencing_type, output_file)	
		name, mate, seq, quality = nil, nil, nil
		counter = -1

		File.open(output_file, 'w') do |output| 
			File.open(input_file, 'r').readlines.each do |line|
				counter += 1
				line = line.strip
			
				if counter % 4 == 0 
					if sequencing_type == 'se'
						name, mate = line.strip, 1
					else
						name, mate = line.strip.split('/')
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
				
					name, mate, seq, quality = nil, nil, nil
					counter = -1
				end 
			end
		end
		$logfile.puts "#{Time.new.strftime("%c")}: Anchor preparation succeded."	
	end


	# 1st mapping-round (anchors) via bowtie2 system call.
	#
	#
	# bowtie_index  - Path to bowtie index and index base name.
	# input_file    - Input file with anchors in fasta- or fastq-format.
	# input_format	- Input format of reads, can be fasta or fastq.
	# output_file   - Name of output file.
	#
	# Return mapped anchors in bam-format.
	def bowtie_map1(bowtie_index, input_file, input_format, output_file)
		if input_format == 'fasta'
			stdin, stdout, stderr, t = Open3.popen3("bowtie2 -x #{bowtie_index} -f -U #{input_file} | samtools view -bS - > #{output_file}")
		else
			stdin, stdout, stderr, t = Open3.popen3("bowtie2 -x #{bowtie_index} -q -U #{input_file} | samtools view -bS - > #{output_file}")
		end
		system_exitcode(t, stderr, 'Bowtie2')
	end
	
	
	# 2nd mapping-round (unmapped reads) via bowtie2 system call.
	#
	#
	# bowtie_index  - Path to bowtie index and index base name
	# input_file    - Input file with unmapped reads in fasta- or fastq-format.
	# input_format	- Input format of reads, can be fasta or fastq.
	# output_file   - Name of output file
	#
	# Return mapped anchors in bam-format.
	def bowtie_map2(bowtie_index, input_file, input_format, output_file)
		if input_format == 'fasta'
			stdin, stdout, stderr, t = Open3.popen3("bowtie2 -x #{bowtie_index} -f -U #{input_file} --no-unal | samtools view -bS - > #{output_file}")
		else
			stdin, stdout, stderr, t = Open3.popen3("bowtie2 -x #{bowtie_index} -q -U #{input_file} --no-unal | samtools view -bS - > #{output_file}")
		end	
		system_exitcode(t, stderr, 'Bowtie2')
	end

	
	# Built bowtie-index from candidate loci.
	#
	#
	# input_file - Input file.
	# logfile    - Name of logfile to write to.
	#
	# Return tab-delimited file with candidate loci.
	def bowtie_build(input_file, prefix)
		stdin, stdout, stderr, t = Open3.popen3("bowtie2-build -f #{input_file} #{prefix}")
		system_exitcode(t, stderr, 'Building bowtie2 index')
	end
	
	
	# Read singletons.
 	# Do directly from accepted_hits.bam or from separate file??
 	#
 	#
 	# singletons	- Text-file with singletons (keep as bam-file?)
 	# read_length	- Length of sequencing read.
 	#
 	# Return hash with all singletons, key = qname, values = mapping information.
	def read_singletons(singletons, read_length)
		single_reads = {}
		
		File.open(singletons, 'r').readlines.each do |line|
  
  		line = line.strip.split(/\s+/)
  		qname, flag, chr, start = line[0..3]  	
  		flag.to_i & 0x10 > 0 ? strand = -1 : strand = 1
  		cigar = line[5]
			distance = genomic_mappinglength(cigar, read_length)
			
			if distance != false
				strand == 1 ? stop = start + distance : stop = start - distance
				single_reads[qname] = [chr, start, stop, strand]
			end
		end
		single_reads
	end

end