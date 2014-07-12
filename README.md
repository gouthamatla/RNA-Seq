RNA-Seq
=======

Automated TopHat-cuffdiff pipeline for illumina paired-end data.

This script performs:

1. Alignment to the genome with tophat2
2. Then Removes the reads with mapping quality less than 50 ( unique mapping )
3. Sort and index the bam files.
4. Run cuffdiff based on experimental grouping.

This script takes 

1. Path to set of input fastq files ( Illumina paired-end files with illumina naming convention).
2. Path to genome fasta file.
3. Path and base name of Bowtie2 index files
4. Path to GTF file
5. Number of processors
6. Experimental grouping file

This script assumes:

1. cufflinks package v2 installed and added to Path
2. samtools is installed and added to path

Experimental grouping file should have the R1 files and the group number seperated by tab...an example file is given below:

c1r1_R1_001.fastq.gz    1
c1r2_R1_001.fastq.gz    1
c1r3_R1_001.fastq.gz    1
c2r1_R1_001.fastq.gz    2
c2r2_R1_001.fastq.gz    2
c2r3_R1_001.fastq.gz    2


Usage:

1.Open the script in any text editor and give full path to all the reguired files/folders.
2. Then run the script using the following command:
    
    python RNA-SEQ.py

