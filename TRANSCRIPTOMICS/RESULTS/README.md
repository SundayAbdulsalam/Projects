# 🧬 RNA-seq Analysis of Preconception Tributyltin Exposure in Mouse Gonadal White Adipose Tissue

## Project Overview
This project investigates the transcriptomic consequences of maternal preconception exposure to metabolism-disrupting environmental factors in first-generation (F1) mouse offspring.
The study uses RNA-seq data to investigate whether exposure before conception is associated with persistent changes in gene expression and biological pathways related to metabolism and disease.

The analysis focuses on metabolically important tissues:
- Gonadal white adipose tissue (GWAT)

Samples were obtained from 11-week-old F1 female and male mice.

---
## Research Question
Does maternal preconception exposure to environmental metabolism disruptors alter gene expression and biological pathways in F1 offspring?

The analysis aims to identify:
- Differentially expressed genes
- Altered biological processes
- Dysregulated molecular pathways
- Tissue-specific transcriptomic responses
- Sex-specific transcriptional responses
- Potential metabolic and disease-related signatures
---

**Female GWAT:**

- 3 DMSO control samples
- 3 TBT-treated samples

The analysis therefore focuses specifically on the transcriptional difference between the TBT and DMSO groups.

---

# Experimental Comparison

| Group | Treatment | Replicates | Tissue | Sex |
|---|---|---:|---|---|
| Control | DMSO | 3 | GWAT | Female |
| Treatment | 50 nM TBT | 3 | GWAT | Female |

### Sample accessions

| Sample | Group |
|---|---|
| SRR37636080 | DMSO |
| SRR37636083 | DMSO |
| SRR37636086 | DMSO |
| SRR37636110 | 50 nM TBT |
| SRR37636113 | 50 nM TBT |
| SRR37636116 | 50 nM TBT |

---

# Analysis Workflow

```text
FASTQ
  ↓
FastQC
  ↓
HISAT2
  ↓
SAM
  ↓
Sorted BAM
  ↓
BAM indexing/QC
  ↓
featureCounts
  ↓
Gene count matrix
  ↓
DESeq2
  ↓
Differential expression
  ↓
PCA / MA plot
  ↓
Volcano plot
  ↓
Heatmap
  ↓
DEG annotation
  ↓
GO / KEGG enrichment
  ↓
Biological interpretation
________________________________________
Tools and Technologies
Operating System
•	Ubuntu Linux 
•	Windows Subsystem for Linux (WSL) 
Command-line tools
•	FastQC 
•	HISAT2 
•	Samtools 
•	featureCounts 
Reference
•	Mus musculus GRCm39 
Statistical analysis
•	R 
•	DESeq2 
•	Bioconductor 
Visualization
•	ggplot2 
•	pheatmap 
Planned functional analysis
•	clusterProfiler 
•	org.Mm.eg.db 
•	enrichplot 
________________________________________
1. Quality Control
FastQC was used to evaluate the quality of the raw paired-end sequencing reads.
Quality metrics examined included:
•	Per-base sequence quality 
•	Per-sequence quality 
•	Per-base sequence content 
•	GC content 
•	Sequence duplication 
•	Adapter contamination 
•	Sequence length distribution 
Example command:
fastqc SRR37636083_1.fastq \
SRR37636083_2.fastq \
-o qc_report \
-t 2
________________________________________
2. Reference Genome
The mouse reference genome used for alignment was:
Mus musculus GRCm39
Reference files:
Mus_musculus.GRCm39.dna.primary_assembly.fa
Mus_musculus.GRCm39.116.gtf
The FASTA file provides the genomic sequence used by HISAT2.
The GTF annotation provides gene and exon coordinates used by featureCounts.
________________________________________
3. HISAT2 Indexing
The reference genome was indexed using HISAT2:
hisat2-build \
Mus_musculus.GRCm39.dna.primary_assembly.fa \
hisat2_indexfiles/ref_genome_index
________________________________________
4. RNA-seq Alignment
Paired-end reads were aligned to the GRCm39 reference genome using HISAT2.
Example:
hisat2 -p 2 \
-x hisat2_indexfiles/ref_genome_index \
-1 SRR37636083_1.fastq \
-2 SRR37636083_2.fastq \
-S female_gwat_dmso_2.sam
Six samples were aligned:
3 DMSO
3 TBT
________________________________________
5. SAM/BAM Processing
SAM files were converted to sorted BAM files using Samtools.
Example:
samtools sort -@ 2 \
-o female_gwat_dmso_2.sorted.bam \
female_gwat_dmso_2.sam
BAM files were indexed:
samtools index female_gwat_dmso_2.sorted.bam
Sorted and indexed BAM files were used for downstream quantification.
________________________________________
6. Gene Quantification
Gene-level read counts were generated using featureCounts.
The analysis used the mouse GRCm39.116 GTF annotation.
Important featureCounts settings included:
-p
--countReadPairs
-B
-C
-t exon
-g gene_id
Command:
featureCounts -T 2 -p --countReadPairs -B -C \
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
________________________________________
7. Differential Expression Analysis
Differential expression was performed using DESeq2.
The experimental design was:
~ condition
with:
DMSO = control
TBT = treatment
The comparison was:
TBT vs DMSO
Genes were considered significantly differentially expressed when:
adjusted p-value < 0.05
________________________________________
Results
Differential Expression
The analysis tested:
12,093 genes
and identified:
24 significant DEGs
Distribution:
Upregulated:     4
Downregulated:  20
This indicates that the selected TBT exposure condition was associated with predominantly decreased expression among the statistically significant genes detected in this small dataset.
________________________________________
Selected Differentially Expressed Genes
Gene	log2FC	Direction
Ccl21a	+2.67	Upregulated
Slc39a8	-1.64	Downregulated
Mt1	-2.11	Downregulated
Prg4	-1.61	Downregulated
Tead1	-1.71	Downregulated
Mt2	-3.07	Downregulated
Pdk4	-2.08	Downregulated
Anpep	+1.01	Upregulated
Cyp2f2	+1.28	Upregulated
Ank3	+2.91	Upregulated
________________________________________
DEG Annotation
The 24 significant DEGs were annotated with:
•	Ensembl gene ID 
•	Gene symbol 
•	Gene biotype 
•	Base mean expression 
•	log2 fold change 
•	p-value 
•	adjusted p-value 
•	Direction of regulation 
All 24 significant genes were successfully annotated.
Example:
ENSMUSG00000094686 → Ccl21a
ENSMUSG00000053897 → Slc39a8
ENSMUSG00000031765 → Mt1
ENSMUSG00000006014 → Prg4
ENSMUSG00000055320 → Tead1
________________________________________
Visualization
The analysis includes:
PCA
Used to evaluate sample-level expression patterns and assess whether biological replicates cluster according to experimental condition.
MA Plot
Used to visualize the relationship between mean expression and log2 fold change.
Volcano Plot
Used to visualize:
•	Magnitude of differential expression 
•	Statistical significance 
•	Upregulated genes 
•	Downregulated genes 
Heatmap
Used to visualize expression patterns of significant DEGs across the six samples.
________________________________________
Output Files
Important generated files include:
gene_counts_paired.txt

DEG_TBT_vs_DMSO.csv

DEG_TBT_vs_DMSO_significant.csv

visualizations/PCA_TBT_vs_DMSO.pdf

visualizations/MAplot_TBT_vs_DMSO.pdf

visualizations/Volcano_TBT_vs_DMSO.pdf

visualizations/DEG_heatmap.pdf

visualizations/Significant_DEGs_final.csv

annotation/DEG_TBT_vs_DMSO_annotated.csv

annotation/DEG_TBT_vs_DMSO_portfolio.csv
________________________________________
Limitations
This project intentionally uses a small subset of the original experimental dataset for computational portfolio purposes.
Important limitations include:
1.	Only six samples were analyzed. 
2.	The analysis contains three biological replicates per condition. 
3.	Only one tissue was investigated. 
4.	Only female GWAT was analyzed. 
5.	The comparison is limited to DMSO versus 50 nM TBT. 
6.	The analysis does not represent the complete experimental design. 
7.	Differential expression findings should therefore be interpreted as exploratory rather than definitive biological conclusions. 
8.	Functional enrichment results may be sensitive to the relatively small number of significant DEGs. 
________________________________________
Current Status
Completed
•	FASTQ acquisition 
•	FastQC 
•	Reference genome preparation 
•	HISAT2 indexing 
•	HISAT2 alignment 
•	SAM/BAM processing 
•	BAM indexing 
•	BAM QC 
•	featureCounts 
•	Gene count matrix 
•	DESeq2 
•	Differential expression 
•	PCA 
•	MA plot 
•	Volcano plot 
•	Heatmap 
•	DEG annotation

Skills Demonstrated
This project demonstrates practical experience with:
Linux / Command Line
•	File management 
•	Bash scripting 
•	Linux pipelines 
•	Process management 
•	Bioinformatics command-line tools 
NGS Data Analysis
•	FASTQ processing 
•	Quality control 
•	RNA-seq alignment 
•	SAM/BAM manipulation 
•	Genome annotation 
•	Read quantification 
Computational Biology
•	Reference genome indexing 
•	Transcriptomic analysis 
•	Differential expression 
•	Functional annotation 
•	Biological interpretation 
R / Bioconductor
•	DESeq2 
•	ggplot2 
•	pheatmap 
•	Bioconductor package management 
•	Statistical analysis 
•	Data visualization 
________________________________________
Conclusion
This project demonstrates an end-to-end bulk RNA-seq workflow from raw sequencing reads to differential gene expression and gene annotation.
The focused analysis of female mouse GWAT identified 24 statistically significant genes following comparison of 50 nM TBT-treated samples with DMSO controls, with 20 genes showing decreased expression and 4 showing increased expression.
The next stage of the project will investigate whether these differentially expressed genes converge on biological processes or pathways relevant to metabolism, adipose tissue biology, environmental exposure responses, and related cellular processes.
________________________________________
Author
Sunday Abdulsalam
Bioinformatics | Computational Biology | Data Analytics
Interested in:
•	Transcriptomics 
•	Cancer genomics 
•	Computational biology 
•	AI for health 
•	Data science 
•	Molecular biology 
________________________________________
Disclaimer
This repository represents a computational portfolio analysis using a selected subset of publicly available RNA-seq data. The analysis is intended to demonstrate bioinformatics methodology and should not be interpreted as a complete reproduction of the original experimental study or as definitive evidence of biological causation.

