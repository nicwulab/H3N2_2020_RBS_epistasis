## Mapping mutations that can confer H3N2 A/Italy/11871/2020 (Italy20) with L194P compatibility

### Obejctives
We introduced 10 amino acid substitutions found betweeen H3N2 A/Italy/11871/2020 (Italy20) and H3N2 A/Singapore/INFIMH-16-0019/2016 (Sing16) to generate all possible combinations (1024 variants) and probe for mutations that can confer compatibility with egg-adaptive mutation L194P.

### Dependencies
* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [R](https://www.r-project.org/) (version 4.1)
* [PEAR](https://github.com/tseemann/PEAR) (Zhang et al., PMID: 24142950)

### Input files
* [./fasta/Italy20HA_mutlib_ref.fasta](./fasta/Italy20HA_mutlib_ref.fasta): Reference amino acid seqeunce of Italy20HA regions of interests (contains L194P backgroud)
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

### Filtering the enrichment results
1. Filter out the variants that contains L194P
``python3 script/Italy20_HA_filter.py``   
    - Input files:
      - [./results/Ita20HA_MultiMutLib.tsv](/results/Ita20HA_MultiMutLib.tsv)
    - Output files:
      - [./results/Ita20HA_MultiMutLib_filtered.tsv](/results/Ita20HA_MultiMutLib_filtered.tsv)

2. Find variants that are highly enriched across both replicates
``Rscript script/Italy20_HA_CutOff.R``   
    - Input files:
      - [./results/Ita20HA_MultiMutLib_filtered.tsv](/results/Ita20HA_MultiMutLib_filtered.tsv)
    - Output files:
      - [./results/mutation_count_size4.tsv](/results/mutation_count_size4.tsv)
    - [./results/mutation_count_size43.tsv](/results/mutation_count_size43.tsv)

### Ploting the enrichment data
1. Plot the enrichment data across replicates
``Rscript script/Italy20_plot_compare_rep.R``   
    - Input files:
      - [./results/Ita20HA_MultiMutLib_filtered.tsv](/results/Ita20HA_MultiMutLib_filtered.tsv)
    - Output files:
      - [graph/Italy20_mutlib_rep_compare.png](/graph/Italy20_mutlib_rep_compare.png)
