#!/bin/bash

# Purpose:
#   To align R1 & R3 Fastqs as built and corrected by divmux.codon

# Return:
#   Aligned BAM

# Usage:
# ./02_divalign.sh <exp_type> <modality> <sample_name> <Path_R1_correct> <Path_R3_correct> <threads> <path_bwa> <path_bwarefDB>
# e.g.:
# ./02_divalign.sh nanoCNT modA ScKDMA_S1 ScKDMA_S1_R1_001_correct.fastq ScKDMA_S1_R3_001_correct.fastq 20 /home/bwa-mem2-2.2.1_x64-linux/bwa-mem2 /home/refBWAmem2/hg19.fa

# Install bwa-mem2 and index the ref
#   curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
#   ./bwa-mem2-2.2.1_x64-linux/bwa-mem2 index hg19.fa
# Install samtools (if not already installed)
set -e

# Variables
exp_type="${1}";
modality="${2}";
sample_name="${3}";
R1="${4}";
R3="${5}";
threads="${6}";
path_bwa="${7}";
path_bwarefDB="${8}";

if [ ${#@} -lt 8 ] ; then
    printf '\nUsage:\n';
    printf '    02_divalign.sh \\\n';
    printf '        exp_type \\\n';
    printf '        modality \\\n';
    printf '        sample_name \\\n';
    printf '        R1 \\\n';
    printf '        R3 \\\n';
    printf '        threads \\\n';
    printf '        path_bwa \\\n';
    printf '        path_bwarefDB \\\n';
    printf 'Parameters:\n';
    printf '  - exp_type: Experience type (ie. nanoCT).\n';
    printf '  - modality: Modality.\n';
    printf '  - sample_name: Sample name.\n';
    printf '  - R1:   Path to corrected FASTQ R1 (uncompressed).\n';
    printf '  - R3:   Path to corrected FASTQ R3 (uncompressed).\n';
    printf '  - threads: Number of threads to align with.\n';
    printf 'Purpose: Align R1 & R3 as built by Divmux\n';
    printf '         Output is a bam file of name <sample_name.modality.exp_type.bam> \n\n';
    exit 1
fi

RGID="${sample_name}.${modality}.${exp_type}" #sample_name.modality.exp_type
library="${sample_name}.${modality}" #sample_name.modality
platform="ILLUMINA" #technology

# Alignement
${path_bwa} mem ${path_bwarefDB} \
-t ${threads} \
-R "@RG\tID:${RGID}\tSM:${sample_name}\tLB:${library}\tPL:${platform}" \
-C \
${R1} ${R3} > TEMP.sam

awk '
BEGIN {
    FS="\t";    # Set field separator as tab
    OFS="\t";   # Set output field separator as tab
}
{
    # Skip lines starting with "@"
    if (substr($0, 1, 1) == "@") {
        print;  # Print the line as it is
        next;  # Skip to the next line
    }

    # Use the match function to extract the string following RG:Z: and CB:Z:
    match($0, /RG:Z:([^[:space:]]+)/, rg);
    match($0, /CB:Z:([^[:space:]]+)/, cb);

    # Replace the string following RG:Z: with the string following CB:Z:
    sub("RG:Z:" rg[1], "RG:Z:" cb[1]);
    # Print the modified line
    print
}' TEMP.sam > TEMPRG.sam

awk '
BEGIN {
    FS="\t";    # Set field separator as tab
    OFS="\t";   # Set output field separator as tab

    # Declare an associative array to store unique Barcodes
    # This will act as a set to keep track of seen Barcodes
    delete barcodes;
}
{
    # Skip lines starting with "@" (header lines)
    if (substr($0, 1, 1) == "@") {
        if ($1 == "@RG") {
            match($0, /SM:([^[:space:]]+)/, all_sm);
            match($0, /LB:([^[:space:]]+)/, all_lb);
            match($0, /PL:([^[:space:]]+)/, all_pl);
            next;
        } else if ($1 == "@PG") {
            last_PG_line = $0;
            next;
        } else {
            # Print other header lines as they are
            print;
            next;
        }
    }

    # extract the Barcode from CB:Z:
    match($0, /CB:Z:([^[:space:]]+)/, cb);
    if (cb[1]) {
        # Store the extracted Barcode in the associative array
        barcodes[cb[1]] = 1;
    }
}
END {
    # Output the new @RG header lines with unique Barcodes
    for (barcode in barcodes) {
        # Construct the @RG line with the unique Barcode
        rg_line = "@RG\tID:" barcode "\tSM:" all_sm[1] "\tLB:" all_lb[1] "\tPL:" all_pl[1];
        # Print the new @RG line
        print rg_line;
    }
    # Print the last @PG line
    print last_PG_line;
}' TEMP.sam > TEMPHEADER.sam

# Remove lines starting with '@' in A.txt, concatenate B.txt at the header, and output to C.txt
{ cat TEMPHEADER.sam; grep -v '^@' TEMPRG.sam; } | samtools view -bS --threads ${threads} -o ${RGID}.bam -
rm TEMP.sam TEMPRG.sam TEMPHEADER.sam

# Convert Peak BED file to SAF format
awk 'OFS="\t" {print $1"."$2"."$3"."$4, $1, $2+1, $3, "."}' ${inputBEDFile} > ${outputSAFFile}
sed -i '1s/^/GeneID\tChr\tStart\tEnd\tStrand\n/' ${outputSAFFile}
# if you're not confident in your sed version (ie using a mac), change the above line to this:
#echo -e "GeneID\tChr\tStart\tEnd\tStrand" | cat - ${outputSAFFile} > temp && mv temp ${outputSAFFile}

featureCounts -T ${threads} -a <input saf file> -F SAF -p --byReadGroup -O -o ${RGID}.tsv ${RGID}.bam


# UNUSED CODE

# Get alignement stats
#samtools flagstat -@ ${threads} ${RGID}.bam

# Split by BC (WILL CRASH YOUR SYSTEM creates too many files)
#samtools split -@ ${threads} --max-split '-1' -u NoBarcode.sam -f "%*_%\!" -d 'CB' -v ${RGID}.bam

#Get list of unique BCs in SAM file
#samtools view temp.sam | cut -f 12- | tr "\t" "\n"  | grep  "^CB:Z:"  | cut -d ':' -f 3 | sort | uniq > uniq_BCs

# TESTING CODE
macs3 callpeak --treatment TESTDONE.bam \
--format BAMPE \
--gsize hs \
--nomodel \
--nolambda \
--keep-dup auto \
--qvalue 0.01 \
--call-summits \
--outdir callpeak_TEST \
--name merged \
--verbose 2

awk 'OFS="\t" {print $1"."$2"."$3"."$4, $1, $2+1, $3, "."}' TEST.bed > TEST.saf
sed -i '1s/^/GeneID\tChr\tStart\tEnd\tStrand\n/' TEST.saf
featureCounts -T 16 -a TEST.saf -F SAF -p --byReadGroup -O -o TESTDONE.tsv TESTDONE.bam
