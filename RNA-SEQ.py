# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  # #
#
# A Python wrapper to run tophat2 and cuffdiff on paired-end illumina files
#  
#  1.Align the R1 and R2 to given bowtie2 index
#  2.Run Samtools to keep reads with uniqe mapping, sort and index
#  3.Run CuffDiff on sorted files. Experimental grouping information will be taken from groupingInfo file
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


import os,re

#Please give full path
Input_FilePath=""
Path_to_Genome=""
Path_to_GenomeIndex=""
Path_to_GTF=""
No_of_processors=""
groupingInfo=""

###{{{ Check the inputs

if not os.path.isdir(Input_FilePath):
	print "Path to fastq(gz) does not exists"
	exit()
if not os.path.exists(Path_to_Genome):
	print "Genome fasta file does not exists"
	exit()
#if not os.path.exists(Path_to_GenomeIndex):
#	print "Path to genome build does not exists"
#	exit()

if not os.path.exists(Path_to_GTF):
	print "GTF file does not exists in the given path"
	exit()
if not os.path.exists(groupingInfo):
	print "Experimental grouping file does not exists in the given path"
	exit()

###}}}


def parseGroupInfo(groupingInfo):
	Group1=[]
	Group2=[]
	for line in open(groupingInfo,'r'):
		if line.strip().split("\t")[1] == "1":
			Group1.append("sort."+(line.strip().split("\t")[0]).split("_R1_")[0]+".bam")
		elif line.strip().split("\t")[1] == "2":
			Group2.append("sort."+str(line.strip().split("\t")[0]).split("_R1_")[0]+".bam")
		else:
			print "Only two groups are supported in experimental grouping file"
			exit()
			
	return ",".join(Group1),",".join(Group2)

group1,group2=parseGroupInfo(groupingInfo)


listOf_R1_File= sorted([fastq for fastq in os.listdir(Input_FilePath) if re.search("R1",fastq)])

listOf_R2_File= sorted([fastq for fastq in os.listdir(Input_FilePath) if re.search("R2",fastq)])


Paired_Files = zip(listOf_R1_File,listOf_R2_File)


if len(listOf_R1_File) != len(listOf_R2_File):
	print "error: The R1 and R2 files are not in equal number. Program is exiting"
	print "number of R1 files: "+str(len(listOf_R1_File))
	print "number of R2 files: "+str(len(listOf_R2_File))
	exit()


def CheckPairs(pair):
	if pair[0]==pair[1]:
		print "error"
		exit()
	elif str(pair[0])==str(pair[1]).replace("_R2_", "_R1_", 1):
		return True

def output(pair):
	if str(pair[0])==str(pair[1]).replace("_R2_", "_R1_", 1):
		outdir=str(pair[0]).split("_R1_")
		return outdir[0]


tophat_Commands=[]
tophat_Outdirds=[]

for pair in Paired_Files:
	if len(pair) != 2:
		print "error: The R1 and R2 files are not in equal number"
		exit()
	elif CheckPairs(pair):
		Command = "tophat2 -o "+output(pair)+" -p "+No_of_processors+" -G "+Path_to_GTF+" "+Path_to_GenomeIndex+" "+Input_FilePath+"/"+pair[0]+" "+Input_FilePath+"/"+pair[1]+" &> "+output(pair)+".log"
		print Command
		os.system(Command)
		tophat_Commands.append(Command)
		tophat_Outdirds.append(output(pair))
		
'''
for command in tophat_Commands:
	print command

exit()
'''

for dir in tophat_Outdirds:
	filter="samtools view -b -q 50 " + dir+"/accepted_hits.bam"+" -o "+dir+".bam"
	sort="samtools sort "+dir+".bam"+" sort."+dir
	index="samtools index "+" sort."+dir+".bam"
	print filter
	os.system(filter)
	print sort
	os.system(sort)
	print index
	os.system(index)

cuffdiff_cmd="cuffdiff -o CuffDiff_Output -p "+ No_of_processors+" -u -b " +Path_to_Genome+" "+Path_to_GTF+" -L 1,2 "+group1+" "+group2+" &> CuffDiff_Output.log" 
print cuffdiff_cmd

os.system(cuffdiff_cmd)
