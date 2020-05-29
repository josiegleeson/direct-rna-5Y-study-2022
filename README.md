# direct RNA sequencing with ONT for gene expression profiling


This repository contains scripts to produce results and figures presented in the paper.

Please cite: x

Data access: x

<b>Instructions and descriptions of each script:</b></br>
<b>dRNA.sh</b></br>
Usage:</br>
./dRNA.sh sample_genome.fasta sample_transcriptome.fasta sample_annotation.gtf sequin_genome.fasta sequin_transcriptonme.fasta sequin_annotation.gtf pass.fastq sample_name</br>
Description:</br> 
Takes FASTQ file from a sample that includes 5Y and sequin sequence. Maps to both human and sequin genome and transcriptome. Performs gene and transcript quantification. 

<b>dRNA-wrapper.pbs</b></br>
Usage:</br>
qsub dRNA-wrapper.pbs</br>
Description:</br>
The pbs job script used to run dRNA.sh on all samples.

<b>extractDataBam.R</b></br>
Usage:</br> 
Rscript extractDataFromBam.R yourfile.bam gencode.gtf outprefix</br>
Description:</br>
Pulls data from a bam file. Calculates various useful metrics such as % full-length transcripts and filters all primary and secondary alignments to retain only the 'true' alignment per read. (in progress)

<b>filterBamTrueTranscript</b></br>
Usage:</br> 
x</br>
Description:</br>
Filters a bam file of transcriptome primary and secondary alignments (up to 10 secondary alignments with scores equal to primary alignemnt) to keep the best alignment per read. This file can be subsequently used for quantification with Salmon to provide a more accurate representation of transcripts.







