#!/bin/bash

# Purpose:
#   To align R1 & R2 Fastqs as built and corrected by divmux.codon

# Return:
#   Aligned BAM

# Usage:
# ./02_divalign.sh <exp_type> <modality> <sample_name> <Path_R1_correct> <Path_R2_correct> <threads> <path_bwa> <path_bwarefDB> <PathGATK> <RemDups> <PathSamtools> <PathOutputBam> <PathOutputPicardDupStats> <sam_header> <min_good_reads_in_cells>
# e.g.:
# ./02_divalign.sh nanoCNT modA ScKDMA_S1 R1_correct.fq R2_correct.fq 8 /home/ahrmad/bwa-mem2-2.2.1_x64-linux/bwa-mem2 /home/ahrmad/refBWAmem2/hg19.fa /home/ahrmad/gatk-4.5.0.0/gatk true /home/ahrmad/micromamba/envs/ali/bin/samtools /home/ahrmad/testing/TEST.bam /home/ahrmad/testing/TEST_DupMetrics.txt /home/ahrmad/testing/sam_header.txt 200

# Install bwa-mem2 and index the ref
#   curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
#   ./bwa-mem2-2.2.1_x64-linux/bwa-mem2 index hg19.fa
# Install samtools (if not already installed)
# Install gatk (https://github.com/broadinstitute/gatk/releases/latest)
set -e

# Variables
exp_type="${1}";
modality="${2}";
sample_name="${3}";
R1="${4}";
R2="${5}";
threads="${6}";
path_bwa="${7}";
path_bwarefDB="${8}";
PathGATK="${9}";
RemDups="${10}";
PathSamtools="${11}";
PathOutputBam="${12}";
PathOutputPicardDupStats="${13}";
PathSam_header="${14}";
min_good_reads_in_cells="${15}";

if [ ${#@} -lt 14 ]; then
    printf '\nUsage:\n';
    printf '    02_divalign.sh \\\n';
    printf '        exp_type \\\n';
    printf '        modality \\\n';
    printf '        sample_name \\\n';
    printf '        R1 \\\n';
    printf '        R2 \\\n';
    printf '        threads \\\n';
    printf '        path_bwa \\\n';
    printf '        path_bwarefDB \\\n';
    printf '        PathGATK \\\n';
    printf '        RemDups \\\n';
    printf '        PathSamtools \\\n';
    printf '        PathOutputBam \\\n';
    printf '        PathOutputPicardDupStats\n';
    printf '        PathSam_header\n';
    printf '        min_good_reads_in_cells\n';
    printf 'Parameters:\n';
    printf '  - exp_type: Experience type (ie. nanoCT).\n';
    printf '  - modality: Modality.\n';
    printf '  - sample_name: Sample name.\n';
    printf '  - R1:   Path to corrected FASTQ R1 (uncompressed).\n';
    printf '  - R2:   Path to corrected FASTQ R2 (uncompressed).\n';
    printf '  - threads: Number of threads to align with.\n';
    printf '  - path_bwa: Path to BWA binary.\n';
    printf '  - path_bwarefDB: Path to BWA DB w/ extension (ex:hg19.fa).\n';
    printf '  - PathGATK: Path to gatk bin.\n';
    printf '  - RemDups: Remove Duplicates (true) or only mark them (false).\n';
    printf '  - PathSamtools: Path to samtools binary.\n';
    printf '  - PathOutputBam: Path to output BAM file.\n';
    printf '  - PathOutputPicardDupStats: Path to output Picard Duplication Stats file.\n';
    printf '  - PathSam_header: Path to SAM header file as produced by demux.codon.\n';
    printf '  - min_good_reads_in_cells: Minimum number of Properly paired and mapped reads in cells.\n';
    printf 'Purpose: Align R1 & R2 as built by divmux codon script\n';
    printf '         Output is a bam file with marked/removed duplicates and associated duplication metrics \n\n';
    exit 1
fi

RGID="${sample_name}-${modality}-${exp_type}" #sample_name.modality.exp_type
library="${sample_name}.${modality}" #sample_name.modality
platform="ILLUMINA" #technology

#threads split for samtools view and sort
if [ "$threads" -lt 4 ]; then
    view_threads=1
    sort_threads=1
else
    view_threads=$((threads / 4))
    sort_threads=$((threads - view_threads))
fi

# Alignement
echo "Aligning ${R1} and ${R2} with BWA-MEM2"
${path_bwa} mem -t ${threads} -C -o ${RGID}_TEMP.sam ${path_bwarefDB} ${R1} ${R2}

echo "Changing header and sorting alignment..."
# Prepare the header
${PathSamtools} view -H ${RGID}_TEMP.sam > ${RGID}_TEMPHEADER.sam
sed -e "s/.*/&\tSM:$sample_name\tLB:$library\tPL:$platform/" ${PathSam_header} > ${RGID}_RG_HEADER.sam
# Find the line number of the last occurrence of "@SQ"
last_sq_line=$(grep -n "@SQ" "${RGID}_TEMPHEADER.sam" | tail -n1 | cut -d':' -f1)
# Split the header based on the line number of the last "@SQ" occurrence
head -n "$last_sq_line" "${RGID}_TEMPHEADER.sam" > "${RGID}_TEMPHEADER1.sam"
echo "" > "${RGID}_TEMPHEADER2.sam" 
tail -n +"$((last_sq_line + 1))" "${RGID}_TEMPHEADER.sam" >> "${RGID}_TEMPHEADER2.sam"

# Append the new header to the RG-replaced SAM file and convert to BAM
{ cat ${RGID}_TEMPHEADER1.sam ${RGID}_RG_HEADER.sam ${RGID}_TEMPHEADER2.sam; ${PathSamtools} view --threads ${view_threads} ${RGID}_TEMP.sam; } | ${PathSamtools} sort --threads ${sort_threads} -m 2G -n -o ${RGID}_TEMP.bam -

# Mark duplicates
echo "Marking duplicates..."
${PathGATK} MarkDuplicates -I ${RGID}_TEMP.bam -O ${RGID}_TEMP_NoDup.bam -M ${PathOutputPicardDupStats} --REMOVE_DUPLICATES ${RemDups} --BARCODE_TAG CB --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000

echo "Get numbers of properly paired and mapped reads per barcode..."
# Get the number of Properly paired and mapped reads per barcode
${PathSamtools} view --require-flags 0x2 ${RGID}_TEMP_NoDup.bam | awk '/RG:Z:/{match($0, /RG:Z:[ACGT]+/); cnts[substr($0, RSTART+5, RLENGTH-5)]++} END {for (bc in cnts) print cnts[bc] "\t" bc}' > ${RGID}_ProperPairedMapped_reads_per_barcode.tsv

echo "Filtering and sorting..."
# Get the barcodes with more than min_good_reads_in_cells reads
awk -v threshold="$min_good_reads_in_cells" -F'\t' '$1 >= threshold {print $2}' ${RGID}_ProperPairedMapped_reads_per_barcode.tsv > ${RGID}_GoodBarcodes.txt
# Keep only the Properly paired and mapped reads from the barcodes with less than min_good_reads_in_cells reads
${PathSamtools} view --threads ${view_threads} --with-header --require-flags 0x2 --read-group-file ${RGID}_GoodBarcodes.txt ${RGID}_TEMP_NoDup.bam | ${PathSamtools} sort -@ ${sort_threads} -m 2G -o ${PathOutputBam} -
# Index the BAM
${PathSamtools} index --threads ${threads} --bai --output ${PathOutputBam}.bai ${PathOutputBam}

#Remove temp files
echo "Removing temp files..."
rm ${RGID}_TEMP.sam ${RGID}_TEMPHEADER.sam ${RGID}_TEMPHEADER1.sam ${RGID}_TEMPHEADER2.sam ${RGID}_RG_HEADER.sam ${RGID}_TEMP.bam ${RGID}_TEMP_NoDup.bam* ${RGID}_GoodBarcodes.txt

# UNUSED CODE
# Get alignement stats
#samtools flagstat -@ ${threads} ${PathOutputBam}
