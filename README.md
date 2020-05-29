# direct RNA sequencing with ONT for gene expression profiling


This repository contains scripts to produce results and figures presented in the paper.

Please cite:
Paper.

Data access:

Instructions and descriptions of each script:

dRNA.sh
Usage:
./dRNA.sh sample_genome.fasta sample_transcriptome.fasta sample_annotation.gtf sequin_genome.fasta sequin_transcriptonme.fasta sequin_annotation.gtf pass.fastq sample_name
Description: 
Takes FASTQ file from a sample that includes 5Y and sequin sequence. Maps to both human and sequin genome and transcriptome. Performs gene and transcript quantification. 

wrapper-dRNA.pbs
Usage:
x
Description:
The pbs job script used to run dRNA.sh on all samples.

extractDataBam.R
Usage: 
Rscript extractDataBam.R yourfile.bam gencode.gtf outprefix
Description:
Pulls data from a bam file. Calculates various useful metrics such as % full-length transcripts and filters all primary and secondary alignments to retain only the 'true' alignment per read. (in progress)





