#!/bin/bash

#http://homer.salk.edu/homer/basicTutorial/mapping.html

#~/.local/bin/cutadapt --help

STAREX=/home/escott/bin/STAR

#MUSREF=resources/Mus_musculus.GRCm38.dna.primary_assembly.fa
#GENOMEIDX=resources/star-mus_musculus

HUMREF=/home/escott/resources/reference/human_g1k_v37.fasta

#STAR --runMode genomeGenerate \
    #--genomeFastaFiles /home/escott/resources/reference/human_g1k_v37.fasta \
    #--genomeDir /home/escott/resources/reference/star-human_b37 \
    #--runThreadN 4

GENOMEIDX=/home/escott/resources/reference/star-human_b37

fadir=/home/escott/workspace/toe1/raw/new
outdir=/home/escott/workspace/toe1/raw/aln

#while read -r line;
  #do name=$line;
  ## Rum alignment using STAR, this gives SAM file
  #echo "fastq file ${fadir}/${name}_L008_R1_001.fastq.gz"
  #echo "Out prefix: ${outdir}/${name}"
  #$STAREX --genomeDir $GENOMEIDX \
    #--readFilesIn ${fadir}/${name}_L008_R1_001.fastq.gz \
                  #${fadir}/${name}_L008_R2_001.fastq.gz \
    #--runThreadN 16 \
    #--outFileNamePrefix ${outdir}/${name} \
    #--readFilesCommand zcat
#done < '/home/escott/workspace/toe1/raw/SampleID.txt'

while read -r line;
   do name=$line;
   echo "Sam file: ${outdir}/${name}.out.sam"
   samtools view -bS "${outdir}/${name}Aligned.out.sam" -o "${outdir}/${name}.bam"
done < '/home/escott/workspace/toe1/raw/SampleID.txt'


#FASTA=/projects/ps-gleesonlab/mboat7_rnaseq/raw/GleesonKO_1_S5_L006_R1_001.fastq.gz
#FASTA=/projects/ps-gleesonlab/mboat7_rnaseq/raw/GleesonKO_25_S7_L006_R1_001.fastq.gz
#FASTA=/projects/ps-gleesonlab/mboat7_rnaseq/raw/GleesonKO_17_S6_L006_R1_001.fastq.gz
#FASTA=/projects/ps-gleesonlab/mboat7_rnaseq/raw/GleesonKO_29_S8_L006_R1_001.fastq.gz
#FASTA=/projects/ps-gleesonlab/mboat7_rnaseq/raw/GleesonKO_1_S5_L006_R1_001.fastq.gz
#FASTA=/projects/ps-gleesonlab/mboat7_rnaseq/raw/GleesonWT_10_S2_L006_R1_001.fastq.gz
#FASTA=raw/GleesonWT_9_S1_L006_R1_001.fastq.gz
#ALNDIR=aligned/GleesonWT_9_S1_L006_R1_001

#$STAREX --genomeDir $GENOMEIDX \
    #--readFilesIn $FASTA \
    #--runThreadN 16 \
    #--outFileNamePrefix $ALNDIR \
    #--readFilesCommand zcat


#/software/STAR/STAR_2.3.1z12/STAR --genomeDir ~/Homo_sapiens_assembly19_STAR/
#--readFilesIn $name_1.fastq.gz $name_2.fastq.gz \
#--readFilesCommand zcat \
#--runThreadN 20 \
#--outFileNamePrefix $name

#TFILE=raw/GleesonKO_17_S6_L006_R1_001.fastq.gz

#STAR \
 #--runThreadN 6 \
 #--runMode alignReads \
 #--genomeDir /home/escott/workspace/mboat7/index \
 #--readFilesIn $TFILE \
 #--readFilesCommand zcat \
 #--sjdbGTFfile resources/gencode.vM7.annotation.gtf.gz
 ##--genomeFastaFiles resources/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
 ##--sjdbGTFfile resources/gencode.vM7.annotation.gtf.gz

 ##--sjdbOverhang

