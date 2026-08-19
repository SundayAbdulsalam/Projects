# ============================================================
# DESeq2: 50 nM TBT vs DMSO
# Female GWAT RNA-seq
# ============================================================

# Personal R library
.libPaths(c("~/R/library", .libPaths()))

# Load DESeq2
library(DESeq2)

# ------------------------------------------------------------
# 1. Read featureCounts output
# ------------------------------------------------------------

counts_raw <- read.delim(
    "gene_counts_paired.txt",
    header = TRUE,
    comment.char = "#",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Extract gene IDs and sample count columns
# ------------------------------------------------------------

gene_ids <- counts_raw$Geneid

counts <- counts_raw[, 7:12]

# Set clean sample names
colnames(counts) <- c(
    "TBT_1",
    "TBT_2",
    "TBT_3",
    "DMSO_1",
    "DMSO_2",
    "DMSO_3"
)

# Use gene IDs as row names
rownames(counts) <- gene_ids

# Convert counts to numeric matrix
counts <- as.matrix(counts)
mode(counts) <- "numeric"

# ------------------------------------------------------------
# 3. Create sample metadata
# ------------------------------------------------------------

metadata <- data.frame(
    row.names = colnames(counts),
    condition = factor(
        c(
            "TBT",
            "TBT",
            "TBT",
            "DMSO",
            "DMSO",
            "DMSO"
        ),
        levels = c("DMSO", "TBT")
    )
)

# Check sample names
print(metadata)

# ------------------------------------------------------------
# 4. Create DESeq2 object
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ condition
)

# ------------------------------------------------------------
# 5. Filter low-count genes
# ------------------------------------------------------------

keep <- rowSums(counts(dds) >= 10) >= 3

dds <- dds[keep, ]

cat(
    "Genes remaining after filtering:",
    nrow(dds),
    "\n"
)

# ------------------------------------------------------------
# 6. Run DESeq2
# ------------------------------------------------------------

dds <- DESeq(dds)

# ------------------------------------------------------------
# 7. TBT vs DMSO
# ------------------------------------------------------------

res <- results(
    dds,
    contrast = c(
        "condition",
        "TBT",
        "DMSO"
    )
)

# ------------------------------------------------------------
# 8. Convert results to data frame
# ------------------------------------------------------------

res_df <- as.data.frame(res)

# IMPORTANT:
# Add gene IDs from the row names
res_df$gene_id <- rownames(res_df)

# Reorder columns
res_df <- res_df[
    ,
    c(
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj"
    )
]

# Sort by adjusted p-value
res_df <- res_df[
    order(res_df$padj),
    ,
    drop = FALSE
]

# ------------------------------------------------------------
# 9. Save ALL DESeq2 results
# ------------------------------------------------------------

write.csv(
    res_df,
    "DEG_TBT_vs_DMSO.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 10. Identify significant DEGs
# ------------------------------------------------------------

sig <- res_df[
    !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1,
    ,
    drop = FALSE
]

# Save significant genes
write.csv(
    sig,
    "DEG_TBT_vs_DMSO_significant.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 11. Print DEG summary
# ------------------------------------------------------------

cat(
    "\nTotal genes tested:",
    nrow(res_df),
    "\n"
)

cat(
    "Significant DEGs:",
    nrow(sig),
    "\n"
)

cat(
    "Upregulated:",
    sum(sig$log2FoldChange >= 1),
    "\n"
)

cat(
    "Downregulated:",
    sum(sig$log2FoldChange <= -1),
    "\n"
)

# ------------------------------------------------------------
# 12. PCA
# ------------------------------------------------------------

vsd <- vst(
    dds,
    blind = FALSE
)

pdf(
    "PCA_TBT_vs_DMSO.pdf",
    width = 8,
    height = 6
)

print(
    plotPCA(
        vsd,
        intgroup = "condition"
    )
)

dev.off()

# ------------------------------------------------------------
# 13. MA plot
# ------------------------------------------------------------

pdf(
    "MAplot_TBT_vs_DMSO.pdf",
    width = 8,
    height = 6
)

plotMA(
    res,
    ylim = c(-5, 5)
)

dev.off()

# ------------------------------------------------------------
# 14. Save normalized counts
# ------------------------------------------------------------

normalized_counts <- counts(
    dds,
    normalized = TRUE
)

write.csv(
    as.data.frame(normalized_counts),
    "normalized_counts_TBT_vs_DMSO.csv"
)

cat("\nDESeq2 analysis completed successfully.\n")
