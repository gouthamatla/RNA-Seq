#This script takes a query name sorted bam file and fixes the CIGAR string of STAR.
# Eg:
# HISEQ:479:C78UFACXX:1:1101:10371:52836  99      chr1    121133905       255     97M3S   =       121133911
# HISEQ:479:C78UFACXX:1:1101:10371:52836  147     chr1    121133911       255     1S99M   =       121133905

#Will become
# HISEQ:479:C78UFACXX:1:1101:10371:52836  99      chr1    121133905       27      100M    =       121133910
# HISEQ:479:C78UFACXX:1:1101:10371:52836  147     chr1    121133910       27      100M    =       121133905


import pysam,sys,re
BAM = pysam.Samfile(sys.argv[1])

header=BAM.header
outfile = pysam.AlignmentFile(sys.argv[2], "wh", header=header)


def fix_CIGAR(read):
	
	if "S" in read.cigarstring:
		
		#If the read is clipped at the begening and at the end
		if  re.search("^([0-9]+S).*([0-9]+)S$",read.cigarstring):
			soft_clip= re.search("(^[0-9]+)S([0-9]+)(\S)(.*)",read.cigarstring)
			softclipped_bases=soft_clip.group(1)
			following_matches=soft_clip.group(2)
			rest_ofCIGAR=soft_clip.group(4)
			new_startString=int(softclipped_bases)+int(following_matches)
			new_startString=str(new_startString)+soft_clip.group(3)+rest_ofCIGAR
						
			cigar_end=re.search("([0-9]+)(\S)([0-9]+)S$",new_startString)
			softclipped_bases=cigar_end.group(3)
			preceeding_matches=cigar_end.group(1)
			preceeding_letter=cigar_end.group(2)
			new_CIGAR= int(softclipped_bases)+int(preceeding_matches)
			new_CIGAR=str(new_CIGAR)+preceeding_letter
			read.cigarstring=new_CIGAR			
			read.reference_start=read.reference_start-int(softclipped_bases)
		
		#if the read is softclipped only at the begining
		elif  re.search("^([0-9]+S)",read.cigarstring):
			cigar_start= re.search("(^[0-9]+)S([0-9]+)(\S)(.*)",read.cigarstring)
			softclipped_bases=cigar_start.group(1)
			following_matches=cigar_start.group(2)
			rest_ofCIGAR=cigar_start.group(4)
			new_startString=int(softclipped_bases)+int(following_matches)
			new_CIGAR=str(new_startString)+cigar_start.group(3)+rest_ofCIGAR
			read.cigarstring=new_CIGAR
			read.reference_start=read.reference_start-int(softclipped_bases)
		
		elif  re.search("([0-9]+S)$",read.cigarstring):
			cigar_end=re.search("([0-9]+)(\S)([0-9]+)S$",read.cigarstring)
			softclipped_bases=cigar_end.group(3)
			preceeding_matches=cigar_end.group(1)
			preceeding_letter=cigar_end.group(2)
			new_CIGAR= int(softclipped_bases)+int(preceeding_matches)
			new_CIGAR=str(new_CIGAR)+preceeding_letter
			read.cigarstring=new_CIGAR
	
	else:
		return read
	
	return read

for read in BAM:
	
	if read.is_paired and not read.mate_is_unmapped:
		next_read=BAM.next()
		print read.qname, next_read.qname
		read=fix_CIGAR(read)
		next_read=fix_CIGAR(next_read)
		
		#Change the coordinates of next read accordingly
		next_read.next_reference_start=read.reference_start
		read.next_reference_start=next_read.reference_start
		
		outfile.write(read)
		outfile.write(next_read)
		
	elif  read.is_paired and read.mate_is_unmapped or not read.is_paired:
		print "Orphan", read.qname
		read=fix_CIGAR(read)
		read.next_reference_start=read.reference_start
		outfile.write(read)
	
