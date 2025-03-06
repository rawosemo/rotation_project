#!/bin/bash
#============================================================================================================================================
# Comprehensive bioinformatics pipeline to infer parchment source and elucidate the microbiome that colonize the surface of the parchments
#
# This pipeline performs the following steps:
#   1. Quality control using FastQC and MultiQC.
#   2. Adapter trimming with Trimmomatic.
#   3. Alignment of reads to the human genome (hg38 and CHM13).
#   4. Alignment of unmapped reads to animal genomes.
#   5. Extraction of animal reference IDs from genome files (to produce reference_id_mapping_file).
#   6. Counting reads mapped to each animal genome (the most being the source of the parchment).
#   7. Taxonomic classification using Kraken2 and krakenUniq.
#   8. Taxonomic classification using MetaPhlAn.
#   9. Extraction of read IDs for a specific microbe (e.g., Cutibacterium).
#  10. Extraction, conversion, and sampling of reads for BLAST confirmation.
#  11. Alignment of extracted reads to a microbe genome (e.g., C. acnes).
#  12. Final processing: sorting, indexing, and computing coverage/depth.
#
# Each section prints one start message (with a timestamp) and one completion  message (with the total elapsed time in seconds).
#
# Author: Rasaq Awosemo
# Date: 2025-03-04
#==========================================================================================================================================

set -euo pipefail

###############################################################################
# Section 1: FastQC Analysis
# Purpose: Assess the quality of raw FASTQ files and generate summary reports.
###############################################################################
INPUT_DIR="/home5/rbawosem/SharedDrive/sub_large/trimmed"
OUTPUT_DIR="/home5/rbawosem/SharedDrive/sub_large/trimmed/fastQC_results"
LOG_FILE="/home5/rbawosem/SharedDrive/sub_large/trimmed/fastqc_run.log"

mkdir -p "$OUTPUT_DIR"

section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FastQC analysis started." >> "$LOG_FILE"

fastqc -t 8 -o "$OUTPUT_DIR" "$INPUT_DIR"/*_paired.fastq.gz 2>&1 | tee -a "$LOG_FILE"
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR" 2>&1 | tee -a "$LOG_FILE"

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FastQC analysis completed in ${elapsed} seconds." >> "$LOG_FILE"

###############################################################################
# Section 2: Trimmomatic Adapter Trimming
# Purpose: Remove adapter sequences and low-quality bases from raw reads.
###############################################################################
TRIMMOMATIC_JAR="/home5/rbawosem/apps/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS_FILE="/home5/rbawosem/apps/Trimmomatic-0.39/adapters/custom_adapters.fa"
INPUT_DIR="/home5/rbawosem/SharedDrive/sub_large/"
OUTPUT_DIR="/home5/rbawosem/SharedDrive/sub_large/trimmed"
LOG_FILE="${OUTPUT_DIR}/trimming_run.log"

mkdir -p "$OUTPUT_DIR"

section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Trimmomatic adapter trimming started." >> "${LOG_FILE}"

for FORWARD_READ in "$INPUT_DIR"/*_R1_*.fastq.gz; do
  BASENAME=$(basename "$FORWARD_READ" | sed 's/_R1_001.fastq.gz//')
  REVERSE_READ="${INPUT_DIR}/${BASENAME}_R2_001.fastq.gz"
  FORWARD_PAIRED="${OUTPUT_DIR}/${BASENAME}_R1_paired.fastq.gz"
  FORWARD_UNPAIRED="${OUTPUT_DIR}/${BASENAME}_R1_unpaired.fastq.gz"
  REVERSE_PAIRED="${OUTPUT_DIR}/${BASENAME}_R2_paired.fastq.gz"
  REVERSE_UNPAIRED="${OUTPUT_DIR}/${BASENAME}_R2_unpaired.fastq.gz"
  SAMPLE_LOG="${OUTPUT_DIR}/${BASENAME}_trimmomatic.log"

  java -jar "$TRIMMOMATIC_JAR" PE -threads 8 -phred33 \
    "$FORWARD_READ" "$REVERSE_READ" \
    "$FORWARD_PAIRED" "$FORWARD_UNPAIRED" \
    "$REVERSE_PAIRED" "$REVERSE_UNPAIRED" \
    ILLUMINACLIP:"$ADAPTERS_FILE":2:30:10:2 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    > "$SAMPLE_LOG" 2>&1
done

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Trimmomatic adapter trimming completed in ${elapsed} seconds." >> "${LOG_FILE}"
echo "Trimming completed. Trimmed files and logs are in $OUTPUT_DIR"

###############################################################################
# Section 3: Human Genome Alignment
# Purpose: Align reads to the human genome, filter mapped/unmapped reads,
#          convert to FASTQ, compress outputs, and record summary statistics.
###############################################################################
GENOME_DIR="/home5/rbawosem/genomes/hg38"
INPUT_DIR="/home5/rbawosem/pilot_data/trimmed/align_batch/unmapped/realignment/unmapped/"
OUTPUT_DIR="/home5/rbawosem/pilot_data/trimmed/align_batch/unmapped/realignment/unmapped/realignment2"
THREADS=8

MAPPED_DIR="$OUTPUT_DIR/mapped"
UNMAPPED_DIR="$OUTPUT_DIR/unmapped"
FASTQ_DIR="$OUTPUT_DIR/fastq"
LOG_DIR="$OUTPUT_DIR/logs"

mkdir -p "$MAPPED_DIR" "$UNMAPPED_DIR" "$FASTQ_DIR" "$LOG_DIR"

SUMMARY_FILE="$OUTPUT_DIR/alignment_summary.tsv"
echo -e "Sample\tInput Read Pairs\tMapped Reads\tUnmapped Reads" > "$SUMMARY_FILE"

section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Human genome alignment started." >> "$LOG_DIR/alignment_master.log"

for SAMPLE in $(ls "$INPUT_DIR"/*_R1.fastq.gz | sed 's/_R1.fastq.gz//' | xargs -n 1 basename); do
  echo "Processing sample: $SAMPLE"
  R1_IN="${INPUT_DIR}/${SAMPLE}_R1.fastq.gz"
  R2_IN="${INPUT_DIR}/${SAMPLE}_R2.fastq.gz"
  SAM_FILE="${LOG_DIR}/${SAMPLE}_human.sam"
  BAM_FILE="${LOG_DIR}/${SAMPLE}_human.bam"
  MAPPED_BAM="${MAPPED_DIR}/${SAMPLE}_mapped.bam"
  UNMAPPED_BAM="${UNMAPPED_DIR}/${SAMPLE}_unmapped.bam"
  R1_UNMAPPED="${FASTQ_DIR}/${SAMPLE}_unmapped_R1.fastq"
  R2_UNMAPPED="${FASTQ_DIR}/${SAMPLE}_unmapped_R2.fastq"
  R1_MAPPED="${FASTQ_DIR}/${SAMPLE}_mapped_R1.fastq"
  R2_MAPPED="${FASTQ_DIR}/${SAMPLE}_mapped_R2.fastq"
  SAMPLE_LOG="${LOG_DIR}/${SAMPLE}_alignment.log"

  bowtie2 -x "${GENOME_DIR}/CHM13_index" -1 "$R1_IN" -2 "$R2_IN" \
    -S "$SAM_FILE" --very-sensitive --threads "$THREADS" 2>&1 | tee "$SAMPLE_LOG"
  samtools view -bS "$SAM_FILE" > "$BAM_FILE"
  samtools view -b -F 4 "$BAM_FILE" > "$MAPPED_BAM"
  MAPPED_READS=$(samtools view -c "$MAPPED_BAM")
  samtools view -b -f 12 -F 256 "$BAM_FILE" > "$UNMAPPED_BAM"
  samtools fastq "$UNMAPPED_BAM" -1 "$R1_UNMAPPED" -2 "$R2_UNMAPPED"
  samtools fastq "$MAPPED_BAM" -1 "$R1_MAPPED" -2 "$R2_MAPPED"
  gzip "$R1_UNMAPPED" "$R2_UNMAPPED" "$R1_MAPPED" "$R2_MAPPED"
  INPUT_READS=$(grep "reads; of these:" "$SAMPLE_LOG" | awk '{print $1}')
  UNMAPPED_READS=$(zcat "${R1_UNMAPPED}.gz" | wc -l | awk '{print $1 / 4}')
  echo -e "$SAMPLE\t$INPUT_READS\t$MAPPED_READS\t$UNMAPPED_READS" >> "$SUMMARY_FILE"
  mv "${FASTQ_DIR}/${SAMPLE}_unmapped_R1.fastq.gz" "$UNMAPPED_DIR"
  mv "${FASTQ_DIR}/${SAMPLE}_unmapped_R2.fastq.gz" "$UNMAPPED_DIR"
done

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Human genome alignment completed in ${elapsed} seconds." >> "$LOG_DIR/alignment_master.log"

###Note: the unmapped reads from this step were subjected to another human genome alignment
###using CHM13 (all you need to do is change the GENOME_DIR to where CHM13 was indexed)
###############################################################################
# Section 4: Alignment Against Animal Genomes
# Purpose: Align unmapped reads to a combined animal genome database.
###############################################################################
INPUT_DIR="/home5/rbawosem/pilot_data/trimmed/align_batch/unmapped/realignment/unmapped/"
OUTPUT_DIR="$INPUT_DIR/all_genomes"
COMBINED_GENOME="/home5/rbawosem/genomes_combined/combined_genomes"
THREADS=8

mkdir -p "$OUTPUT_DIR/mapped" "$OUTPUT_DIR/unmapped" "$OUTPUT_DIR/logs"

section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Animal genome alignment started." >> "$OUTPUT_DIR/logs/animal_alignment.log"

for R1_FASTQ in "$INPUT_DIR"/*_unmapped_R1.fastq.gz; do
    SAMPLE=$(basename "$R1_FASTQ" _unmapped_R1.fastq.gz)
    R2_FASTQ="$INPUT_DIR/${SAMPLE}_unmapped_R2.fastq.gz"
    ALIGNED_SAM="$OUTPUT_DIR/${SAMPLE}.sam"
    MAPPED_BAM="$OUTPUT_DIR/mapped/${SAMPLE}_mapped.bam"
    UNMAPPED_BAM="$OUTPUT_DIR/unmapped/${SAMPLE}_unmapped.bam"
    MAPPED_R1="$OUTPUT_DIR/mapped/${SAMPLE}_R1.fastq.gz"
    MAPPED_R2="$OUTPUT_DIR/mapped/${SAMPLE}_R2.fastq.gz"
    UNMAPPED_R1="$OUTPUT_DIR/unmapped/${SAMPLE}_unmapped_R1.fastq.gz"
    UNMAPPED_R2="$OUTPUT_DIR/unmapped/${SAMPLE}_unmapped_R2.fastq.gz"
    LOG_FILE="$OUTPUT_DIR/logs/${SAMPLE}_alignment.log"

    bowtie2 -x "$COMBINED_GENOME" -1 "$R1_FASTQ" -2 "$R2_FASTQ" --very-sensitive -S "$ALIGNED_SAM" 2>&1 | tee "$LOG_FILE"
    samtools view -bS "$ALIGNED_SAM" > "$OUTPUT_DIR/${SAMPLE}_aligned.bam"
    samtools view -b -f 2 -q 30 "$OUTPUT_DIR/${SAMPLE}_aligned.bam" > "$MAPPED_BAM"
    samtools fastq "$MAPPED_BAM" -1 "$MAPPED_R1" -2 "$MAPPED_R2"
    samtools view -b -f 12 -F 256 "$OUTPUT_DIR/${SAMPLE}_aligned.bam" > "$UNMAPPED_BAM"
    samtools fastq "$UNMAPPED_BAM" -1 "$UNMAPPED_R1" -2 "$UNMAPPED_R2"
    rm -f "$ALIGNED_SAM" "$OUTPUT_DIR/${SAMPLE}_aligned.bam"
done

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Animal genome alignment completed in ${elapsed} seconds." >> "$OUTPUT_DIR/logs/animal_alignment.log"

###############################################################################
# Section 5: Define Animal reference IDs
# Purpose: Extract reference IDs and species names from genome files.
###############################################################################
MAPPING_FILE="/home5/rbawosem/genomes/species_mapping.txt"
# It contains two column: the first column is the full path to the genome file, and the second column is the species name.
#You can manually create it
OUTPUT_MAPPING="reference_id_mapping_new.tsv"
section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Defining animal reference IDs started." 

while read -r FILE SPECIES; do
    if [[ -f "$FILE" ]]; then
        grep ">" "$FILE" | awk -v sp="$SPECIES" '{print $1"\t"sp}' | sed 's/>//'
    else
        echo "ERROR: File not found - $FILE" >&2
    fi
done < "$MAPPING_FILE" > "$OUTPUT_MAPPING"

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Animal reference IDs defined in ${elapsed} seconds." 
#The file reference_id_mapping_new.tsv then  maps reference IDs (e.g., NW_017219385.1, NC_010443.5) to their corresponding species names (e.g., Capra hircus, Sus scrofa).
###############################################################################
# Section 6: Count Reads Mapped to Each Animal Genome
# Purpose: Summarize the number of primary aligned reads per species.
###############################################################################
bam_dir="/home5/rbawosem/pilot_data/trimmed/align_batch/unmapped/realignment/unmapped/all_genomes/mapped/"
output_file="$bam_dir/primary_alignment_counts_per_species.tsv"
mapping_file="/home5/rbawosem/genomes/reference_id_mapping_new.tsv"

echo -e "Sample\tSpecies\tPrimary_Aligned_Reads" > "$output_file"

section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Counting primary aligned reads per species started." 

for bam_file in "$bam_dir"/*_sorted.bam; do
    sample=$(basename "$bam_file" | sed 's/_sorted.bam//')
    temp_counts=$(mktemp)
    while read -r ref_id species; do
        count=$(samtools view -c -F 260 "$bam_file" "$ref_id")
        echo -e "$species\t$count" >> "$temp_counts"
    done < "$mapping_file"
    awk -F'\t' '{species[$1]+=$2} END {for (s in species) print s, species[s]}' "$temp_counts" | \
    while read -r species total_count; do
        echo -e "$sample\t$species\t$total_count" >> "$output_file"
    done
    rm "$temp_counts"
done

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Read counting completed in ${elapsed} seconds." 

###############################################################################
# Section 7: Kraken2 Classification
# Purpose: Classify reads taxonomically using the Kraken2 database.
###############################################################################
kraken_db="/home5/rbawosem/kraken2_db"
input_dir="/home5/rbawosem/pilot_data/trimmed/align_batch/unmapped/realignment/unmapped/combined_align/unmapped/"
output_dir="/home5/rbawosem/pilot_data/trimmed/align_batch/unmapped/realignment/unmapped/combined_align/unmapped/kraken2_results"
threads=24

mkdir -p "$output_dir"

section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Kraken2 classification started." 

for file_R1 in "$input_dir"/*_unmapped_R1.fastq.gz; do
    sample_id=$(basename "$file_R1" | sed 's/_unmapped_R1.fastq.gz//')
    file_R2="${input_dir}/${sample_id}_unmapped_R2.fastq.gz"
    if [[ ! -f "$file_R2" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: Reverse read file $file_R2 not found. Skipping..."
        continue
    fi
    kraken2 --db "$kraken_db" --paired "$file_R1" "$file_R2" --threads "$threads" \
        --output "$output_dir/${sample_id}_output.txt" \
        --report "$output_dir/${sample_id}_report.txt"
done

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Kraken2 classification completed in ${elapsed} seconds." 
#Do the same for krakenUniq by chnaging the "kraken_db" directory to that of krakenUniq.
###############################################################################
# Section 8: MetaPhlAn Classification
# Purpose: Classify reads taxonomically using the MetaPhlAn database.
###############################################################################
input_dir="/home5/rbawosem/pilot_data/trimmed/align_batch/unmapped/realignment/unmapped/combined_align/unmapped/"
output_dir="/home5/rbawosem/pilot_data/trimmed/align_batch/unmapped/realignment/unmapped/combined_align/unmapped/metaphlan_results/"
metaphlan_db="/home5/rbawosem/metaphlan_db"
threads=8

mkdir -p "$output_dir"

section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] MetaPhlAn classification started." 

for file_R1 in "$input_dir"/*_unmapped_R1.fastq.gz; do
    sample_id=$(basename "$file_R1" | sed 's/_unmapped_R1.fastq.gz//')
    file_R2="${input_dir}/${sample_id}_unmapped_R2.fastq.gz"
    metaphlan "$file_R1,$file_R2" --input_type fastq --bowtie2db "$metaphlan_db" \
        --bowtie2out "${output_dir}/${sample_id}_bowtie2.bz2" --nproc "$threads" \
        --output_file "${output_dir}/${sample_id}_metaphlan.txt"
done

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] MetaPhlAn classification completed in ${elapsed} seconds." 

###############################################################################
# Section 9: Validate Microbe Presence (Extract Read IDs)
# Purpose: Extract read IDs for a specific microbe (e.g., taxID 1747 for Cutibacterium).
###############################################################################
section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting read IDs for validation started." 

for SAMPLE in SG-087_S4_L001 SG-226_S9_L001 SG-112_S5_L001 SG-224_S8_L001 SG-193_S6_L001 SG-200_S7_L001; do
    grep "1747" "$output_dir/${SAMPLE}_output.txt" | cut -f2 > "${SAMPLE}_Cutibacterium_readIDs.txt"
done

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Read ID extraction completed in ${elapsed} seconds." 

###############################################################################
# Section 10: Extract Reads for BLAST Confirmation
# Purpose: From the read IDs, extract corresponding reads, convert them to FASTA,
#          combine paired FASTA files, and sample a subset for BLAST.
###############################################################################
RAW_READS_DIR="/home5/rbawosem/pilot_data/trimmed/align_batch/unmapped/realignment/unmapped/combined_align/unmapped"
READ_ID_DIR="$output_dir"
EXTRACT_OUTPUT="${READ_ID_DIR}/extracted_reads"
FASTA_DIR="${EXTRACT_OUTPUT}/fasta"
COMBINED_DIR="${FASTA_DIR}/combined"
EXTRACTED_DIR="${FASTA_DIR}/extracted"

mkdir -p "$EXTRACT_OUTPUT" "$FASTA_DIR" "$COMBINED_DIR" "$EXTRACTED_DIR"

section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting reads for BLAST confirmation started." 

for READ_ID_FILE in "$READ_ID_DIR"/*Cutibacterium_readIDs.txt; do
    SAMPLE_ID=$(basename "$READ_ID_FILE" | sed 's/_Cutibacterium_readIDs.txt//')
    R1_FASTQ="${RAW_READS_DIR}/${SAMPLE_ID}_unmapped_R1.fastq.gz"
    R2_FASTQ="${RAW_READS_DIR}/${SAMPLE_ID}_unmapped_R2.fastq.gz"
    OUT_R1="${EXTRACT_OUTPUT}/${SAMPLE_ID}_C_acnes_R1.fastq.gz"
    OUT_R2="${EXTRACT_OUTPUT}/${SAMPLE_ID}_C_acnes_R2.fastq.gz"
    if [[ -f "$R1_FASTQ" && -f "$R2_FASTQ" ]]; then
        seqtk subseq "$R1_FASTQ" "$READ_ID_FILE" | gzip > "$OUT_R1" &
        seqtk subseq "$R2_FASTQ" "$READ_ID_FILE" | gzip > "$OUT_R2" &
        wait
        seqtk seq -A "$OUT_R1" > "${FASTA_DIR}/${SAMPLE_ID}_C_acnes_R1.fa" &
        seqtk seq -A "$OUT_R2" > "${FASTA_DIR}/${SAMPLE_ID}_C_acnes_R2.fa" &
        wait
        cat "${FASTA_DIR}/${SAMPLE_ID}_C_acnes_R1.fa" "${FASTA_DIR}/${SAMPLE_ID}_C_acnes_R2.fa" \
            > "${COMBINED_DIR}/${SAMPLE_ID}_C_acnes_combined.fa"
        seqtk sample "${COMBINED_DIR}/${SAMPLE_ID}_C_acnes_combined.fa" 0.2 \
            > "${EXTRACTED_DIR}/${SAMPLE_ID}_C_acnes_sub.fa"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: FASTQ files missing for $SAMPLE_ID, skipping..."
    fi
done

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extraction, conversion, and sampling completed in ${elapsed} seconds." 

###############################################################################
# Section 11: Align Extracted Reads to Microbe Genomes (C. acnes Example)
# Purpose: Align the extracted reads to a specific microbial genome (e.g., C. acnes)
#          and compute coverage/depth statistics.
###############################################################################
REF_DIR="/home5/rbawosem/genomes2/C_acnes_genome/C_acnes"
READS_DIR="${EXTRACT_OUTPUT}"
CA_OUTPUT="${READS_DIR}/C_acnes_analysis"
CA_MAPPED="${CA_OUTPUT}/mapped_reads"
CA_UNMAPPED="${CA_OUTPUT}/Unmapped"

mkdir -p "$CA_MAPPED" "$CA_UNMAPPED"

section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Aligning extracted reads to C. acnes genome started." 

for SAMPLE in $(ls "$READS_DIR" | grep _L001_C_acnes_R1.fastq.gz | sed 's/_L001_C_acnes_R1.fastq.gz//'); do
    R1="${READS_DIR}/${SAMPLE}_L001_C_acnes_R1.fastq.gz"
    R2="${READS_DIR}/${SAMPLE}_L001_C_acnes_R2.fastq.gz"
    SAM_OUT="${CA_OUTPUT}/${SAMPLE}_C_acnes.sam"
    SORTED_BAM="${CA_OUTPUT}/${SAMPLE}_C_acnes_sorted.bam"
    DEPTH_OUT="${CA_OUTPUT}/${SAMPLE}_C_acnes_depth.txt"
    COVERAGE_OUT="${CA_OUTPUT}/${SAMPLE}_C_acnes_coverage.txt"
    MAPPED_BAM="${CA_MAPPED}/${SAMPLE}_C_acnes_mapped.bam"
    UNMAPPED_BAM="${CA_UNMAPPED}/${SAMPLE}_C_acnes_unmapped.bam"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping $SAMPLE: FASTQ files not found!" 
        continue
    fi

    bowtie2 -x "$REF_DIR" -1 "$R1" -2 "$R2" -S "$SAM_OUT"
    samtools view -bS "$SAM_OUT" | samtools sort -o "$SORTED_BAM"
    samtools index "$SORTED_BAM"
    samtools view -b -F 4 "$SORTED_BAM" > "$MAPPED_BAM"
    samtools view -b -f 12 "$SORTED_BAM" > "$UNMAPPED_BAM"
    rm "$UNMAPPED_BAM"
    samtools depth "$SORTED_BAM" > "$DEPTH_OUT"
    samtools coverage "$SORTED_BAM" > "$COVERAGE_OUT"
    AVG_DEPTH=$(awk '{sum+=$3} END {print sum/NR}' "$DEPTH_OUT")
    echo -e "Sample: $SAMPLE\tAverage Depth: $AVG_DEPTH" >> "${CA_OUTPUT}/C_acnes_summary_depth.txt"
    rm "$SAM_OUT"
done

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Alignment to C. acnes genome completed in ${elapsed} seconds." 

###############################################################################
# Section 12: Final Processing: Sorting, Indexing, and Coverage Computation
# Purpose: For each sorted BAM file, compute per-base depth and coverage,
#          saving the outputs with a modified filename that excludes the "_sorted.bam" suffix.
###############################################################################
section_start=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Final processing started." 

for FILE in /path/to/sorted_bam/files/*_sorted.bam; do
    samtools depth "$FILE" > "${FILE%_sorted.bam}_depth.txt"
    samtools coverage "$FILE" > "${FILE%_sorted.bam}_coverage.txt"
done

section_end=$(date +%s)
elapsed=$((section_end - section_start))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Final processing completed in ${elapsed} seconds." 

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All pipeline processes completed successfully."

