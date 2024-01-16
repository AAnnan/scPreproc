## scPreproc
### 1) Demultiplexing (01_divmux.codon):
- Have [**codon**](https://github.com/exaloop/codon#install) and [**seq**](https://github.com/exaloop/seq#installation) installed.

- Run script like this:

`codon run -plugin seq -release -D BC_LEN=16 -D BC_START=29 divmux.codon R2.fastq.gz R1.fastq R3.fastq whitelist.txt stats.tsv 1`

Values for BC_LEN and BC_START need to be passed when invoking codon!

### 2) Alignment (02_divali.sh):
- Have [**bwa**](https://github.com/lh3/bwa), or faster yet, [**bwa-mem2**](https://github.com/bwa-mem2/bwa-mem2) binaries and [**samtools**](https://github.com/samtools/samtools) installed.

- Run script like this:

`./divali.sh exp_type modality sample_name R1.fastq R3.fastq threads path_bwa path_refDB`
