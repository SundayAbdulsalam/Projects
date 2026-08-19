## Getting raw data from NCBI
The prefetch command was used to download the samples from NCBI while the fastq-dump --split-files was used to separate the samples into the forward and reverse reads.
  ### DMSO Samples
```bash
prefetch SRR37636080
fastq-dump --split-files SRR37636080

prefetch SRR37636083
fastq-dump --split-files SRR37636083

prefetch SRR37636086
fastq-dump --split-files SRR37636086
```
### 50TBT Samples
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
  ### DMSO Samples
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
  ### Downloading the reference genome
```bash
wget https://ftp.ensembl.org/pub/current/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
unzip Mus_musculus.GRCm39.dna.primary_assembly.zip
```
  ### Building the index files
```bash
hisat2-build Mus_musculus.GRCm39.dna.primary_assembly.fa hisat2_indexFiles
```
  ### Performing alignment
  After HISAT2 aligns the RNA-seq reads to the mouse reference genome, the alignment information can initially be represented in SAM (Sequence Alignment/Map) format. However, SAM files are plain-text files and can become very large for RNA-seq datasets.
For this project, the SAM output from HISAT2 was therefore piped directly into samtools sort, which converts and sorts the alignments into BAM (Binary Alignment/Map) format.
  ### DMSO Samples
```bash
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636080_1.fastq -2 SRR37636080_2.fastq | samtools sort -@ 2 -o female_gwat_dmso_1.sorted.ba  
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636083_1.fastq -2 SRR37636083_2.fastq | samtools sort -@ 2 -o female_gwat_dmso_2.sorted.ba  
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636086_1.fastq -2 SRR37636086_2.fastq | samtools sort -@ 2 -o female_gwat_dmso_3.sorted.ba
```
   ### 50TBT Samples
```bash
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636110_1.fastq -2 SRR37636110_2.fastq | samtools sort -@ 2 -o female_gwat_50tbt_1.sorted.ba
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636113_1.fastq -2 SRR37636113_2.fastq | samtools sort -@ 2 -o female_gwat_50tbt_2.sorted.ba
hisat2 -p 2 -x hisat2_indexfiles/ref_genome_index -1 SRR37636116_1.fastq -2 SRR37636116_2.fastq | samtools sort -@ 2 -o female_gwat_50tbt_3.sorted.ba
```
  ### BAM File Indexing
After generating the coordinate-sorted BAM files, each BAM file was
indexed using `samtools index`.

```bash
for bam in alignment_bamfile/*.sorted.bam
do
    samtools index "$bam"
done
```
### Gene-Level Read Quantification with featureCounts

Following alignment QC, gene-level read counts were generated from the
sorted BAM files using `featureCounts`.

```bash
featureCounts -T 2 \
-p \
--countReadPairs \
-B \
-C \
-t exon \
-g gene_id \
-a annotation/Mus_musculus.GRCm39.116.gtf \
-o gene_counts_paired.txt \
alignment_bamfile/female_gwat_50tbt_1.sorted.bam \
alignment_bamfile/female_gwat_50tbt_2.sorted.bam \
alignment_bamfile/female_gwat_50tbt_3.sorted.bam \
alignment_bamfile/female_gwat_dmso_1.sorted.bam \
alignment_bamfile/female_gwat_dmso_2.sorted.bam \
alignment_bamfile/female_gwat_dmso_3.sorted.bam
```
  ##Generation of DESeq2 Count Matrix

Following gene-level read quantification with featureCounts, the relevant columns containing gene identifiers and read counts for all six RNA-seq samples were extracted to create a clean count matrix for downstream differential expression analysis.
```bash
cut -f1,7-12 gene_counts_paired.txt > counts_matrix.tsv
```
  ##
Creation of Sample Metadata
A sample metadata file was created to define the sample names and experimental conditions for the six RNA-seq samples. This metadata is required by DESeq2 to correctly assign samples to the treatment and control groups during differential expression analysis.
```bash
cat > sample_metadata.csv << 'EOF'
sample,condition
50TBT_1,50TBT
50TBT_2,50TBT
50TBT_3,50TBT
DMSO_1,DMSO
DMSO_2,DMSO
DMSO_3,DMSO
EOF
```
