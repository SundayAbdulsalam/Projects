## Getting raw data from NCBI
The prefetch command was used to download the samples from NCBI while the fastq-dump --split-files was used to separate the samples into the forward and reverse reads.
  ## DMSO Samples
```bash
prefetch SRR37636080
fastq-dump --split-files SRR37636080

prefetch SRR37636083
fastq-dump --split-files SRR37636083

prefetch SRR37636086
fastq-dump --split-files SRR37636086
```
## 50TBT Samples
```bash
prefetch SRR37636110
fastq-dump --split-files SRR37636110

prefetch SRR37636113
fastq-dump --split-files SRR37636113

prefetch SRR37636116
fastq-dump --split-files SRR37636116
```

## Quality control
Quality control was performed on each samples forward and reverse reads, the results from the quality control is then stored in qc_report file
  ## DMSO Samples
```bash  
fastqc SRR37636080_1.fastq SRR37636080_2.fastq -o qc_report
fastqc SRR37636083_1.fastq SRR37636083_2.fastq -o qc_report
fastqc SRR37636086_1.fastq SRR37636086_2.fastq -o qc_report

  ## 50TBT Samples
fastqc SRR37636110_1.fastq SRR37636110_2.fastq -o qc_report
fastqc SRR37636113_1.fastq SRR37636113_2.fastq -o qc_report
fastqc SRR37636116_1.fastq SRR37636116_2.fastq -o qc_report
```
## ALIGNMENT
Mus_musculus.GRCm39.dna.primary_assembly.fa.gz used as the reference genome for the study was collected from ensembl database  (https://ftp.ensembl.org/pub/current/fasta/mus_musculus/dna/), reference index file was built from the reference genome using hisat2, then subsequent alignment was performed on the samples.
  ## Downloading the reference genome
```bash
wget https://ftp.ensembl.org/pub/current/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
unzip Mus_musculus.GRCm39.dna.primary_assembly.zip
```
  ## Building the index files
```bash
hisat2-build Mus_musculus.GRCm39.dna.primary_assembly.fa hisat2_indexFiles
```
  ## Performing alignment
  ## DMSO Samples
```bash
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636080_1.fastq -2 SRR37636080_2.fastq | samtools sort -@ 2 -o female_gwat_dmso_1.sorted.ba  
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636083_1.fastq -2 SRR37636083_2.fastq | samtools sort -@ 2 -o female_gwat_dmso_1.sorted.ba  
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636086_1.fastq -2 SRR37636086_2.fastq | samtools sort -@ 2 -o female_gwat_dmso_1.sorted.ba
```
   ## 50TBT Samples
```bash
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636110_1.fastq -2 SRR37636110_2.fastq | samtools sort -@ 2 -o female_gwat_50tbt_1.sorted.ba
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636113_1.fastq -2 SRR37636113_2.fastq | samtools sort -@ 2 -o female_gwat_50tbt_1.sorted.ba
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636116_1.fastq -2 SRR37636116_2.fastq | samtools sort -@ 2 -o female_gwat_50tbt_1.sorted.ba
```
