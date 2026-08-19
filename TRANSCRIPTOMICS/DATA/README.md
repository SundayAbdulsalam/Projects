## Getting raw data from NCBI
The prefetch command was used to download the samples from NCBI while the fastq-dump --split-files was used to separate the samples into the forward and reverse reads.
  ## DMSO Samples
prefetch SRR37636080
fastq-dump --split-files SRR37636080
prefetch SRR37636083
fastq-dump --split-files SRR37636083
prefetch SRR37636086
fastq-dump --split-files SRR37636083
  ## 50TBT Samples
prefetch SRR37636110
fastq-dump --split-files SRR37636110
prefetch SRR37636113
fastq-dump --split-files SRR37636113
prefetch SRR37636116
fastq-dump --split-files SRR37636116

## Quality control
Quality control was performed on each samples forward and reverse reads, the results from the quality control is then stored in qc_report file
  ## DMSO Samples  
fastqc SRR37636080_1.fastq SRR37636080_2.fastq -o qc_report
fastqc SRR37636083_1.fastq SRR37636083_2.fastq -o qc_report
fastqc SRR37636086_1.fastq SRR37636086_2.fastq -o qc_report

  ## 50TBT Samples
fastqc SRR37636110_1.fastq SRR37636110_2.fastq -o qc_report
fastqc SRR37636113_1.fastq SRR37636113_2.fastq -o qc_report
fastqc SRR37636116_1.fastq SRR37636116_2.fastq -o qc_report







