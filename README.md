## Studying the Epistasis Effect within H3N2 Italy 20 HA Receptor Binding Site using Combinatorial Mutaional Scanning

### Dependencies
* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [R](https://www.r-project.org/) (version 4.1)
* [PEAR](https://github.com/tseemann/PEAR)

### Input files
* [./fasta/Italy20HA_mutlib_ref.fasta](./fasta/Italy20HA_multilib_ref.fasta): Reference amino acid seqeunce of Italy20HA regions of interests (contains L194P backgroud)
* Raw read files in fastq format from NIH SRA database [BioProject PRJNA883249](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA883249)

### Merging the Illumina NGS raw reads using PEAR
1. Merge overlapping paired-end reads using [PEAR](https://github.com/tseemann/PEAR)   
``pear -f [FASTQ FILE FOR FORWARD READ] -r [FASTQ FILE FOR FORWARD READ] -o [OUTPUT FASTQ FILE]``   
    - Output files should be placed in a folder named fastq_merged/

### Counting and calculate enrichment for both replicates
1. Trim off 5' and 3' flanking region, count variants based on nucleotide sequences, then translate the nucleotide seqeunces into amino acid seqeunces, and finally identify mutations
``python3 script/Italy20_HA_fastq2enrich.py``   
    - Input files:
      - Merged read files in fastq_merged/ folder
      - [./fasta/Italy20HA_mutlib_ref.fasta](./fasta/Italy20HA_mutlib_ref.fasta)
    - Output files:
      - [./results/Ita20HA_MultiMutLib.tsv](/results/Ita20HA_MultiMutLib.tsv)
