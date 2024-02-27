## scPreproc
### 1) Demultiplexing (01_divmux.codon):
- Have [**codon**](https://github.com/exaloop/codon#install) and [**seq**](https://github.com/exaloop/seq#installation) installed.

- Run script like this:

`codon run -plugin seq -release -D BC_LEN=16 -D BC_START=29 -D MIN_READS_IN_BC=200 01_divmux.codon <Path_R2_gzipped> <Path_R1_unzipped> <Path_R3_unzipped> <Path_whitelist> <Path_output_stats.tsv> <Path_reads_per_barcode_stats.tsv> <Hamming_distance_allowed>`

e.g.:

`codon run -plugin seq -release -D BC_LEN=16 -D BC_START=29 -D MIN_READS_IN_BC=200 01_divmux.codon ScKDMA_S1_R2_001.fastq.gz ScKDMA_S1_R1_001.fastq ScKDMA_S1_R3_001.fastq whitelist_10xATAC_737K-cratac-v1.txt ScKDMA_S1_BC_stats.tsv ScKDMA_S1_reads_per_barcode.tsv 1`

Values for BC_LEN, BC_START and MIN_READS_IN_BC need to be passed when invoking codon!

- Output (in script dir):
`<Path_R1_unzipped>_correct.fastq`
`<Path_R3_unzipped>_correct.fastq`
`<Path_output_stats.tsv>`
`<Path_reads_per_barcode_stats.tsv>`

### 2) Alignment (02_divalign.sh):
- Have [**bwa**](https://github.com/lh3/bwa), or faster yet, [**bwa-mem2**](https://github.com/bwa-mem2/bwa-mem2) binaries and [**samtools**](https://github.com/samtools/samtools) installed.

- Run script like this:

`./02_divalign.sh <exp_type> <modality> <sample_name> <Path_R1_correct> <Path_R3_correct> <threads> <path_bwa> <path_bwarefDB`

e.g.:

`./02_divalign.sh nanoCNT modA ScKDMA_S1 ScKDMA_S1_R1_001_correct.fastq ScKDMA_S1_R3_001_correct.fastq 20 /home/bwa-mem2-2.2.1_x64-linux/bwa-mem2 /home/refBWAmem2/hg19.fa`

- Output (in script dir):
`<sample_name>-<modality>-<exp_type>_MarkedDup.bam`
`<sample_name>-<modality>-<exp_type>_DupMetrics.txt`

