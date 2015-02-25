import gzip
import numpy,random

#10000 should not have any foldchange
#3000 should have 10 fold change in condition 1
#2000 should have 10 fold change in condition 2
#200 should be expressed only in condition 1
#200 should express only in condition 2

foldchange_array=['10000,1','3000,10','2000,-10','200,+','200,-']
true_fold_change=open("true_fc.txt","w")

num_reads=100
num_replicates=[3,3]
dispersion=0.3


transcriptome=gzip.open("101_transcriptome_chr1.fa.gz")
transcript_list=[]
transcript_dict={}
fc_dict={}

for line in transcriptome:
    line=line.strip()
    if line.startswith(">"):
        transcript_list.append(line.replace(">",""))

number_of_transcript = sum(int(item.split(",")[0]) for item in foldchange_array)
numpy.random.shuffle(transcript_list)
select_randome_transcripts=transcript_list[:number_of_transcript]


line_count=0
start=0
end=0

for item in range(len(foldchange_array)):
    transcript_num=int(foldchange_array[item].split(",")[0])
    numpy.random.shuffle(transcript_list)
    fold_change=foldchange_array[item].split(",")[1]
    end+=transcript_num
    
    for transcript in select_randome_transcripts[start:end]:
        fc_dict[transcript]=fold_change
        true_fold_change.write(str(transcript+"\t"+fold_change))
        true_fold_change.write("\n")
    
        line_count+=1
    start += end-start


def get_NB(transcript,fc):
    if fc=="+":
        #The transcript is expressed only in first condition
        condition1=numpy.random.negative_binomial(num_reads, dispersion,num_replicates[0])
        condition2=[0 for x in range(num_replicates[1])]
        counts=list(condition1)+condition2
        print transcript+"\t"+"\t".join(str(v) for v in counts)
    elif fc=="-":
        #The transcript is expressed only in first condition
        condition1=numpy.random.negative_binomial(num_reads, dispersion,num_replicates[0])
        condition2=["0" for x in range(num_replicates[1])]
        counts=list(condition1)+condition2
        print transcript+"\t"+"\t".join(str(v) for v in counts)
    elif int(fc)==1:
        condition1=numpy.random.negative_binomial(num_reads, dispersion,num_replicates[0])
        condition2=numpy.random.negative_binomial(num_reads, dispersion,num_replicates[1])
        counts=list(condition1)+list(condition2)
        print transcript+"\t"+"\t".join(str(v) for v in counts)
    elif int(fc)>1:
        condition1=numpy.random.negative_binomial(num_reads*int(fc), dispersion,num_replicates[0])
        condition2=numpy.random.negative_binomial(num_reads, dispersion,num_replicates[1]) 
        counts=list(condition1)+list(condition2)
        print transcript+"\t"+"\t".join(str(v) for v in counts)
    elif int(fc)<0:
        condition1=numpy.random.negative_binomial(num_reads, 0.3,num_replicates[0])
        condition2=numpy.random.negative_binomial(num_reads*abs(int(fc)), 0.3,num_replicates[1]) 
        counts=list(condition1)+list(condition2)
        print transcript+"\t"+"\t".join(str(v) for v in counts)        



for key,value in fc_dict.items():
    if isinstance(value,int):
        get_NB(key,int(value))
    else:
        get_NB(key,str(value))
