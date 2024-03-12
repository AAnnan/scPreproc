## scPreproc
### 1) Demultiplexing (01_divmux.codon):
- Have [**codon**](https://github.com/exaloop/codon#install) and [**seq**](https://github.com/exaloop/seq#installation) installed.

- Run script like this:

`codon run -plugin seq -release -D BC_LEN=16 -D BC_START=29 -D MIN_READS_IN_BC=200 01_divmux.codon <Path_I1_gzipped_or_not> <Path_R1_unzipped> <Path_R2_unzipped> <Path_whitelist> <Path_output_stats.tsv> <Path_reads_per_barcode_stats.tsv> <sam_header.txt> <Hamming_distance_allowed>`

e.g.:

`codon run -plugin seq -release -D BC_LEN=16 -D BC_START=29 -D MIN_READS_IN_BC=200 01_divmux.codon ScKDMA_S1_R2_001.fastq.gz ScKDMA_S1_R1_001.fastq ScKDMA_S1_R3_001.fastq whitelist_10xATAC_737K-cratac-v1.txt ScKDMA_S1_BC_stats.tsv ScKDMA_S1_reads_per_barcode.tsv sam_header.txt 1`

Values for BC_LEN, BC_START and MIN_READS_IN_BC need to be passed when invoking codon!

- Output:
`<Path_R1_unzipped>_correct.fastq`
`<Path_R3_unzipped>_correct.fastq`
`<Path_output_stats.tsv>`
`<Path_reads_per_barcode_stats.tsv>`
`<sam_header.txt>`

### 2) Alignment (02_divalign.sh):
- Uses [**bwa**](https://github.com/lh3/bwa), or faster yet, [**bwa-mem2**](https://github.com/bwa-mem2/bwa-mem2) binaries, [**samtools**](https://github.com/samtools/samtools) and [**GATK**](https://github.com/broadinstitute/gatk/releases/latest).

- Run script like this:

`./02_divalign.sh <exp_type> <modality> <sample_name> <Path_R1_correct> <Path_R2_correct> <threads> <path_bwa> <path_bwarefDB> <PathGATK> <RemDups> <PathSamtools> <PathOutputBam> <PathOutputPicardDupStats> <sam_header>  <min_good_reads_in_cells>`

e.g.:

`./02_divalign.sh nanoCNT modA ScKDMA_S1 R1_correct.fq R2_correct.fq 8 /home/bwa-mem2-2.2.1_x64-linux/bwa-mem2 /home/refBWAmem2/hg19.fa /home/gatk-4.5.0.0/gatk true /home/micromamba/envs/ali/bin/samtools /home/testing/TEST.bam /home/testing/TEST_DupMetrics.txt /home/testing/sam_header.txt 200`

- Output:
`<PathOutputBam>`
`<PathOutputBam>.bai`
`<PathOutputPicardDupStats>`
In script directory:
`<sample_name>-<modality>-<exp_type>_ProperPairedMapped_reads_per_barcode.tsv`

