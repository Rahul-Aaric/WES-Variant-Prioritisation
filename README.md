# WES-Variant-Prioritisation
This repository contains shell scripts and command-line workflows for Whole Exome Sequencing (WES) analysis from raw FASTQ reads to variant calling, annotation, filtering, and variant prioritisation.

## Workflow Overview
Raw FASTQ quality control --> Read trimming using fastp/Trimmomatic --> Reference genome indexing --> Alignment using BWA-MEM --> SAM to BAM conversion and sorting --> Read group addition --> Duplicate marking --> Base Quality Score Recalibration --> Germline variant calling using GATK HaplotypeCaller --> Somatic variant calling using GATK Mutect2 --> VCF merging --> Annotation using SnpEff, SnpSift, and ANNOVAR --> Population frequency and functional filtering --> Variant prioritisation

## Tools Used
fastp
FastQC
Trimmomatic
BWA
SAMtools
GATK
bcftools
SnpEff
SnpSift
ANNOVAR
bedtools

## Reference Genome
Human GRCh38 / hg38 reference genome was used.

Points to keep in mind:- 
a) Modify all paths according to local installation.
b) GRCh38/hg38 reference genome was used.
c) Pipeline supports both germline and somatic WES analysis.
d) Recommended to use compressed FASTQ files (*.fastq.gz).
