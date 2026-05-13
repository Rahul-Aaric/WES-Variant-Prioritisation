#!/bin/bash

set -euo pipefail

#############################################
# WES Variant Prioritisation Pipeline
# From FASTQ to variant annotation/filtering
# Author: Rahul Aaric
#############################################

THREADS=20
JAVA_MEM="24G"

# --------- USER PATHS: EDIT THESE ----------
RAW_DIR="/path/to/raw_fastq"
TRIM_DIR="/path/to/trimmed"
QC_DIR="/path/to/qc"
ALIGN_DIR="/path/to/alignment"
VCF_DIR="/path/to/vcf"
ANNOT_DIR="/path/to/annotation"
INDEX_DIR="/path/to/index"

REFERENCE="${INDEX_DIR}/human_ref38_patch13.fna"
TARGET_BED="/path/to/Twist_ComprehensiveExome_targets_hg38.bed"

GATK="/path/to/gatk"
TRIMMOMATIC="/usr/share/java/trimmomatic.jar"
SNPEFF="/path/to/snpEff.jar"
SNPSIFT="/path/to/SnpSift.jar"
ANNOVAR="/path/to/annovar/table_annovar.pl"
ANNOVAR_DB="/path/to/annovar/humandb"

DBSNP="/path/to/dbsnp151.vcf"
GNOMAD="/path/to/af-only-gnomad.hg38.vcf.gz"
PON="/path/to/pon.vcf.gz"

mkdir -p "$TRIM_DIR" "$QC_DIR" "$ALIGN_DIR" "$VCF_DIR" "$ANNOT_DIR" "$INDEX_DIR"

#############################################
# 1. Reference genome download and indexing
#############################################

prepare_reference() {
    cd "$INDEX_DIR"

    if [ ! -f "$REFERENCE" ]; then
        wget https://ftp-trace.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz

        gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz
        mv GCF_000001405.39_GRCh38.p13_genomic.fna human_ref38_patch13.fna
    fi

    bwa index "$REFERENCE"
    samtools faidx "$REFERENCE"

    "$GATK" CreateSequenceDictionary \
        -R "$REFERENCE" \
        -O "${REFERENCE%.fna}.dict"
}

#############################################
# 2. FASTQ quality control
#############################################

run_fastqc_raw() {
    fastqc "$RAW_DIR"/*.fastq.gz -t "$THREADS" -o "$QC_DIR"
}

#############################################
# 3. Read trimming using fastp
#############################################

run_fastp_paired() {
    for R1 in "$RAW_DIR"/*_R1*.fastq.gz
    do
        SAMPLE=$(basename "$R1" | sed 's/_R1.*.fastq.gz//')
        R2=$(echo "$R1" | sed 's/_R1/_R2/')

        echo "Running fastp for $SAMPLE"

        fastp \
            -i "$R1" \
            -I "$R2" \
            -o "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
            -O "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz" \
            -w "$THREADS" \
            -h "$QC_DIR/${SAMPLE}_fastp.html" \
            -j "$QC_DIR/${SAMPLE}_fastp.json"
    done
}

#############################################
# 4. Alternative trimming using Trimmomatic
#############################################

run_trimmomatic_paired() {
    for R1 in "$RAW_DIR"/*_R1*.fastq.gz
    do
        SAMPLE=$(basename "$R1" | sed 's/_R1.*.fastq.gz//')
        R2=$(echo "$R1" | sed 's/_R1/_R2/')

        echo "Running Trimmomatic for $SAMPLE"

        java -jar "$TRIMMOMATIC" PE -phred33 \
            -threads "$THREADS" \
            "$R1" "$R2" \
            "$TRIM_DIR/${SAMPLE}_R1_paired.fastq.gz" \
            "$TRIM_DIR/${SAMPLE}_R1_unpaired.fastq.gz" \
            "$TRIM_DIR/${SAMPLE}_R2_paired.fastq.gz" \
            "$TRIM_DIR/${SAMPLE}_R2_unpaired.fastq.gz" \
            ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 \
            LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35
    done
}

#############################################
# 5. FASTQC after trimming
#############################################

run_fastqc_trimmed() {
    fastqc "$TRIM_DIR"/*.fastq.gz -t "$THREADS" -o "$QC_DIR"
}

#############################################
# 6. Alignment using BWA-MEM
#############################################

run_alignment() {
    for R1 in "$TRIM_DIR"/*_R1_trimmed.fastq.gz
    do
        SAMPLE=$(basename "$R1" | sed 's/_R1_trimmed.fastq.gz//')
        R2="$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz"

        echo "Aligning $SAMPLE"

        bwa mem -t "$THREADS" "$REFERENCE" "$R1" "$R2" \
            > "$ALIGN_DIR/${SAMPLE}.sam"
    done
}

#############################################
# 7. SAM to sorted BAM
#############################################

sort_bam() {
    for SAM in "$ALIGN_DIR"/*.sam
    do
        SAMPLE=$(basename "$SAM" .sam)

        echo "Sorting $SAMPLE"

        "$GATK" --java-options "-Xmx${JAVA_MEM}" SortSam \
            --INPUT "$SAM" \
            --OUTPUT "$ALIGN_DIR/${SAMPLE}_sorted.bam" \
            --SORT_ORDER coordinate
    done
}

#############################################
# 8. Add read groups
#############################################

add_read_groups() {
    for BAM in "$ALIGN_DIR"/*_sorted.bam
    do
        SAMPLE=$(basename "$BAM" _sorted.bam)

        echo "Adding read groups for $SAMPLE"

        "$GATK" --java-options "-Xmx${JAVA_MEM}" AddOrReplaceReadGroups \
            --INPUT "$BAM" \
            --OUTPUT "$ALIGN_DIR/${SAMPLE}_readgroup.bam" \
            --RGLB lib1 \
            --RGPL illumina \
            --RGPU unit1 \
            --RGSM "$SAMPLE"
    done
}

#############################################
# 9. Mark duplicates
#############################################

mark_duplicates() {
    for BAM in "$ALIGN_DIR"/*_readgroup.bam
    do
        SAMPLE=$(basename "$BAM" _readgroup.bam)

        echo "Marking duplicates for $SAMPLE"

        "$GATK" --java-options "-Xmx${JAVA_MEM}" MarkDuplicates \
            --INPUT "$BAM" \
            --OUTPUT "$ALIGN_DIR/${SAMPLE}_dedup.bam" \
            --METRICS_FILE "$ALIGN_DIR/${SAMPLE}_duplicate_metrics.txt" \
            --VALIDATION_STRINGENCY LENIENT
    done
}

#############################################
# 10. Build BAM index
#############################################

index_bam() {
    for BAM in "$ALIGN_DIR"/*_dedup.bam
    do
        "$GATK" BuildBamIndex \
            --INPUT "$BAM"
    done
}

#############################################
# 11. Alignment and insert size metrics
#############################################

collect_metrics() {
    for BAM in "$ALIGN_DIR"/*_dedup.bam
    do
        SAMPLE=$(basename "$BAM" _dedup.bam)

        "$GATK" --java-options "-Xmx${JAVA_MEM}" CollectAlignmentSummaryMetrics \
            --REFERENCE_SEQUENCE "$REFERENCE" \
            --INPUT "$BAM" \
            --OUTPUT "$QC_DIR/${SAMPLE}_alignment_metrics.txt" \
            --VALIDATION_STRINGENCY LENIENT

        "$GATK" --java-options "-Xmx${JAVA_MEM}" CollectInsertSizeMetrics \
            --INPUT "$BAM" \
            --OUTPUT "$QC_DIR/${SAMPLE}_insert_metrics.txt" \
            --Histogram_FILE "$QC_DIR/${SAMPLE}_insert_size_histogram.pdf"
    done
}

#############################################
# 12. Depth and coverage
#############################################

calculate_depth_coverage() {
    for BAM in "$ALIGN_DIR"/*_dedup.bam
    do
        SAMPLE=$(basename "$BAM" _dedup.bam)

        samtools depth -b "$TARGET_BED" "$BAM" | \
        awk '{sum+=$3} END {print "Average =",sum/NR}' \
        > "$QC_DIR/${SAMPLE}_average_depth.txt"

        samtools depth -b "$TARGET_BED" "$BAM" | \
        awk '{c++; if($3>0) total+=1} END {print (total/c)*100}' \
        > "$QC_DIR/${SAMPLE}_genome_coverage.txt"
    done
}

#############################################
# 13. Base Quality Score Recalibration
#############################################

run_bqsr() {
    for BAM in "$ALIGN_DIR"/*_dedup.bam
    do
        SAMPLE=$(basename "$BAM" _dedup.bam)

        echo "Running BQSR for $SAMPLE"

        "$GATK" --java-options "-Xmx${JAVA_MEM}" BaseRecalibrator \
            -R "$REFERENCE" \
            -I "$BAM" \
            --known-sites "$DBSNP" \
            -O "$ALIGN_DIR/${SAMPLE}_recal_data.table"

        "$GATK" --java-options "-Xmx${JAVA_MEM}" ApplyBQSR \
            -R "$REFERENCE" \
            -I "$BAM" \
            --bqsr-recal-file "$ALIGN_DIR/${SAMPLE}_recal_data.table" \
            -O "$ALIGN_DIR/${SAMPLE}_recal.bam"

        "$GATK" --java-options "-Xmx${JAVA_MEM}" BaseRecalibrator \
            -R "$REFERENCE" \
            -I "$ALIGN_DIR/${SAMPLE}_recal.bam" \
            --known-sites "$DBSNP" \
            -O "$ALIGN_DIR/${SAMPLE}_post_recal_data.table"

        "$GATK" --java-options "-Xmx${JAVA_MEM}" AnalyzeCovariates \
            -before "$ALIGN_DIR/${SAMPLE}_recal_data.table" \
            -after "$ALIGN_DIR/${SAMPLE}_post_recal_data.table" \
            -plots "$QC_DIR/${SAMPLE}_recalibration_plots.pdf"
    done
}

#############################################
# 14. Germline variant calling
#############################################

run_haplotypecaller() {
    for BAM in "$ALIGN_DIR"/*_recal.bam
    do
        SAMPLE=$(basename "$BAM" _recal.bam)

        echo "Calling germline variants for $SAMPLE"

        "$GATK" --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
            -R "$REFERENCE" \
            -I "$BAM" \
            -O "$VCF_DIR/${SAMPLE}_germline_raw.vcf"
    done
}

#############################################
# 15. Compress and index VCF files
#############################################

compress_index_vcf() {
    for VCF in "$VCF_DIR"/*_germline_raw.vcf
    do
        SAMPLE=$(basename "$VCF" _germline_raw.vcf)

        bcftools view "$VCF" -Oz -o "$VCF_DIR/${SAMPLE}_variants.vcf.gz"
        bcftools index "$VCF_DIR/${SAMPLE}_variants.vcf.gz"
    done
}

#############################################
# 16. Merge germline VCF files
#############################################

merge_vcfs() {
    bcftools merge \
        --force-samples \
        --missing-to-ref \
        "$VCF_DIR"/*_variants.vcf.gz \
        -o "$VCF_DIR/all_samples_merged.vcf"
}

#############################################
# 17. Germline annotation
#############################################

annotate_germline() {
    java -Xmx64G -jar "$SNPEFF" \
        -v GRCh38.99 \
        "$VCF_DIR/all_samples_merged.vcf" \
        > "$ANNOT_DIR/variants_snpeff.vcf"

    java -Xmx64G -jar "$SNPSIFT" annotate \
        "$DBSNP" \
        "$ANNOT_DIR/variants_snpeff.vcf" \
        > "$ANNOT_DIR/variants_snpeff_dbsnp.vcf"

    perl "$ANNOVAR" \
        "$ANNOT_DIR/variants_snpeff_dbsnp.vcf" \
        "$ANNOVAR_DB" \
        -buildver hg38 \
        -out "$ANNOT_DIR/germline_annovar" \
        -remove \
        -protocol refGene,cytoBand,1000g2015aug_all,gnomad211_exome,clinvar_20220320,cosmic70 \
        -operation g,r,f,f,f,f \
        -nastring . \
        -vcfinput \
        -polish
}

#############################################
# 18. Prepare BAM files for somatic analysis
#############################################

prepare_somatic_readgroups() {
    NORMAL_LIST="normal_sample_list.txt"
    TUMOR_LIST="tumor_sample_list.txt"

    while read SAMPLE
    do
        "$GATK" --java-options "-Xmx${JAVA_MEM}" AddOrReplaceReadGroups \
            --INPUT "$ALIGN_DIR/${SAMPLE}_recal.bam" \
            --OUTPUT "$ALIGN_DIR/${SAMPLE}_recal_SM.bam" \
            --RGLB lib1 \
            --RGPL illumina \
            --RGPU unit1 \
            --RGSM "NORMAL_${SAMPLE}"

        "$GATK" BuildBamIndex \
            --INPUT "$ALIGN_DIR/${SAMPLE}_recal_SM.bam"
    done < "$NORMAL_LIST"

    while read SAMPLE
    do
        "$GATK" --java-options "-Xmx${JAVA_MEM}" AddOrReplaceReadGroups \
            --INPUT "$ALIGN_DIR/${SAMPLE}_recal.bam" \
            --OUTPUT "$ALIGN_DIR/${SAMPLE}_recal_SM.bam" \
            --RGLB lib1 \
            --RGPL illumina \
            --RGPU unit1 \
            --RGSM "TUMOR_${SAMPLE}"

        "$GATK" BuildBamIndex \
            --INPUT "$ALIGN_DIR/${SAMPLE}_recal_SM.bam"
    done < "$TUMOR_LIST"
}

#############################################
# 19. Somatic variant calling using Mutect2
#############################################

run_mutect2() {
    PAIR_LIST="tumor_normal_pairs.txt"

    while read TUMOR NORMAL
    do
        echo "Running Mutect2 for tumor $TUMOR and normal $NORMAL"

        "$GATK" --java-options "-Xmx${JAVA_MEM}" Mutect2 \
            --native-pair-hmm-threads "$THREADS" \
            -R "$REFERENCE" \
            -I "$ALIGN_DIR/${TUMOR}_recal_SM.bam" \
            -I "$ALIGN_DIR/${NORMAL}_recal_SM.bam" \
            -tumor "TUMOR_${TUMOR}" \
            -normal "NORMAL_${NORMAL}" \
            --germline-resource "$GNOMAD" \
            --panel-of-normals "$PON" \
            --genotype-pon-sites true \
            --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
            -O "$VCF_DIR/${TUMOR}_${NORMAL}_mutect2_raw.vcf.gz"

        "$GATK" --java-options "-Xmx${JAVA_MEM}" FilterMutectCalls \
            -R "$REFERENCE" \
            -V "$VCF_DIR/${TUMOR}_${NORMAL}_mutect2_raw.vcf.gz" \
            -O "$VCF_DIR/${TUMOR}_${NORMAL}_mutect2_filtered.vcf.gz"

    done < "$PAIR_LIST"
}

#############################################
# 20. Extract PASS somatic variants
#############################################

extract_pass_somatic() {
    for VCF in "$VCF_DIR"/*_mutect2_filtered.vcf.gz
    do
        SAMPLE=$(basename "$VCF" _mutect2_filtered.vcf.gz)

        bcftools view -f PASS "$VCF" \
            -o "$VCF_DIR/${SAMPLE}_PASS.vcf"
    done
}

#############################################
# 21. Somatic annotation and filtering
#############################################

annotate_filter_somatic() {
    for VCF in "$VCF_DIR"/*_PASS.vcf
    do
        SAMPLE=$(basename "$VCF" _PASS.vcf)

        java -Xmx24G -jar "$SNPSIFT" annotate \
            "$DBSNP" \
            "$VCF" \
            > "$ANNOT_DIR/${SAMPLE}_dbsnp.vcf"

        java -Xmx4G -jar "$SNPSIFT" filter \
            "((! exists MAF) | (MAF <= 0.01))" \
            "$ANNOT_DIR/${SAMPLE}_dbsnp.vcf" \
            > "$ANNOT_DIR/${SAMPLE}_MAF_0.01.vcf"

        java -Xmx4G -jar "$SNPSIFT" filter \
            "((! exists SAS_AF) | (SAS_AF <= 0.01))" \
            "$ANNOT_DIR/${SAMPLE}_MAF_0.01.vcf" \
            > "$ANNOT_DIR/${SAMPLE}_MAF_0.01_SAS_filtered.vcf"

        perl "$ANNOVAR" \
            "$ANNOT_DIR/${SAMPLE}_MAF_0.01_SAS_filtered.vcf" \
            "$ANNOVAR_DB" \
            -buildver hg38 \
            -out "$ANNOT_DIR/${SAMPLE}_annovar" \
            -remove \
            -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp41a,clinvar_20190305,cosmic90 \
            -operation g,r,f,f,f,f,f \
            -nastring . \
            -vcfinput \
            -polish
    done
}

#############################################
# MAIN PIPELINE
#############################################

echo "Starting WES Variant Prioritisation Pipeline"

# Uncomment the steps you want to run.

# prepare_reference
# run_fastqc_raw
# run_fastp_paired
# run_fastqc_trimmed
# run_alignment
# sort_bam
# add_read_groups
# mark_duplicates
# index_bam
# collect_metrics
# calculate_depth_coverage
# run_bqsr
# run_haplotypecaller
# compress_index_vcf
# merge_vcfs
# annotate_germline

# For somatic analysis:
# prepare_somatic_readgroups
# run_mutect2
# extract_pass_somatic
# annotate_filter_somatic

echo "Pipeline completed"
