import gzip,numpy
from Bio import SeqIO

read_count_matrix = open("test_matrix.txt","r")

number_of_samples=6

transcript_dict={}

for seq_record in SeqIO.parse(gzip.open("101_transcriptome_chr1.fa.gz"), "fasta"):
    transcript_dict[seq_record.id]=seq_record.seq

def single_num(mean,sd):
    num = None
    while True:
        loc = mean
        scale = sd
        num = numpy.random.normal(loc=loc,scale=scale)

        if True:
            break

    return num

def generate_fasta(sample_name,count,transcript_id,read1,read2):
    print sample_name
    sample_R1=open("sample_"+sample_name+"_R1.fastq","a")
    sample_R1.write(">"+"sample_"+sample_name+"_"+transcript_id+"_"+str(count)+"/1")
    sample_R1.write("\n")
    sample_R1.write(read1)
    sample_R1.write("\n")
    sample_R2=open("sample_"+sample_name+"_R2.fastq","a")
    sample_R2.write(">"+"sample_"+sample_name+"_"+transcript_id+"_"+str(count)+"/1")
    sample_R2.write("\n")
    sample_R2.write(read1)
    sample_R2.write("\n")


frag_len=open("frag_len.txt","w")


def get_NormallyDistributedInsertLengths(sample_name,transcript_name,copy_number):
    
    transcript_length=len(transcript_dict[transcript_name])
    transcript_seq=transcript_dict[transcript_name]
    
    if transcript_length>202:
        
        insert_length=[int(x) for x in numpy.random.normal(loc=50,scale=20,size=copy_number)]        
        
        count=0
        for insert_len in insert_length:
            break_point=numpy.random.randint(0, high=(transcript_length-202))
            fragment=(transcript_seq[break_point:break_point+insert_len+202])
            frag_len.write(str(len(fragment))+"\n")
            read1=fragment[0:101]
            read2=fragment[-101:]
            generate_fasta(str(sample_name),count,transcript_name, str(read1), str(read2))
            count+=1
        
    elif transcript_length<=202:
        count=0
        for insert_len in range(copy_number):
            break_point=numpy.random.randint(0, high=transcript_length-101)
            fragment=(transcript_seq[break_point:])
            frag_len.write(str(len(fragment))+"\n")
            read1=fragment[0:101]
            read2=fragment[-101:]
            generate_fasta(str(sample_name),count, transcript_name, str(read1),str(read2))
            count+=1
            
for sample in range(1,number_of_samples+1):
    sample_name="Sample_"+str(sample)    
    
    for transcript in read_count_matrix:
        transcript=transcript.strip()
        transcripts_matrix=transcript.split("\t")    
        transcript_name=transcripts_matrix[0]
        transcript_count=transcripts_matrix[int(sample)]
        #print sample_name,transcript_name,len(transcript_dict[transcript_name]),transcript_count
        get_NormallyDistributedInsertLengths(sample,transcript_name,int(transcript_count))
    read_count_matrix.seek(0)
