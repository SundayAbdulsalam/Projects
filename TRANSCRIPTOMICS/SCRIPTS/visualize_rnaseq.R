# ============================================================
# RNA-seq Visualization
# TBT vs DMSO - Female GWAT
# ============================================================

# -----------------------------
# 1. Load packages
# -----------------------------

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# -----------------------------
# 2. Read DESeq2 results
# -----------------------------

res <- read.csv(
  "DEG_TBT_vs_DMSO.csv",
  stringsAsFactors = FALSE
)

# Check columns
print(colnames(res))

# -----------------------------
# 3. Clean gene IDs
# -----------------------------

# If gene_id is missing, use first column
if (!"gene_id" %in% colnames(res)) {
    colnames(res)[1] <- "gene_id"
}

# Remove Ensembl version numbers if present
res$gene_id <- sub(
    "\\..*$",
    "",
    res$gene_id
)

# -----------------------------
# 4. Define significance
# -----------------------------

res$significance <- "Not significant"

res$significance[
    !is.na(res$padj) &
    res$padj < 0.05 &
    res$log2FoldChange > 1
] <- "Upregulated"

res$significance[
    !is.na(res$padj) &
    res$padj < 0.05 &
    res$log2FoldChange < -1
] <- "Downregulated"

# Count categories
print(table(res$significance))

# -----------------------------
# 5. Create output directory
# -----------------------------

dir.create(
    "visualizations",
    showWarnings = FALSE
)

# ============================================================
# 6. VOLCANO PLOT
# ============================================================

res$minus_log10_padj <- -log10(res$padj)

# Prevent Inf values
res$minus_log10_padj[
    is.infinite(res$minus_log10_padj)
] <- NA

volcano <- ggplot(
    res,
    aes(
        x = log2FoldChange,
        y = minus_log10_padj
    )
) +

geom_point(
    aes(color = significance),
    alpha = 0.7,
    size = 2
) +

geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
) +

geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
) +

labs(
    title = "Volcano Plot: TBT vs DMSO",
    subtitle = "Female Gonadal White Adipose Tissue",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Expression"
) +

theme_classic() +

theme(
    plot.title = element_text(
        face = "bold",
        size = 16
    ),
    plot.subtitle = element_text(
        size = 12
    )
)

ggsave(
    "visualizations/Volcano_TBT_vs_DMSO.png",
    volcano,
    width = 8,
    height = 6,
    dpi = 300
)

ggsave(
    "visualizations/Volcano_TBT_vs_DMSO.pdf",
    volcano,
    width = 8,
    height = 6
)

# ============================================================
# 7. LABEL TOP SIGNIFICANT GENES
# ============================================================

# Select top 10 significant genes
top_genes <- res[
    !is.na(res$padj) &
    res$padj < 0.05,
]

top_genes <- top_genes[
    order(top_genes$padj),
]

top_genes <- head(
    top_genes,
    10
)

print(top_genes[, c(
    "gene_id",
    "log2FoldChange",
    "padj"
)])

# ============================================================
# 8. MA PLOT
# ============================================================

ma_plot <- ggplot(
    res,
    aes(
        x = log10(baseMean + 1),
        y = log2FoldChange
    )
) +

geom_point(
    aes(color = significance),
    alpha = 0.6,
    size = 1.8
) +

geom_hline(
    yintercept = 0,
    linetype = "dashed"
) +

labs(
    title = "MA Plot: TBT vs DMSO",
    x = "log10 Mean Expression",
    y = "log2 Fold Change",
    color = "Expression"
) +

theme_classic()

ggsave(
    "visualizations/MA_TBT_vs_DMSO.png",
    ma_plot,
    width = 8,
    height = 6,
    dpi = 300
)

ggsave(
    "visualizations/MA_TBT_vs_DMSO.pdf",
    ma_plot,
    width = 8,
    height = 6
)

# ============================================================
# 9. HEATMAP
# ============================================================

# Read count matrix
counts <- read.delim(
    "gene_counts_paired.txt",
    comment.char = "#",
    check.names = FALSE
)

# Remove annotation columns
count_matrix <- counts[, 7:ncol(counts)]

# Clean sample names
sample_names <- c(
    "TBT_1",
    "TBT_2",
    "TBT_3",
    "DMSO_1",
    "DMSO_2",
    "DMSO_3"
)

colnames(count_matrix) <- sample_names

# Gene IDs
rownames(count_matrix) <- sub(
    "\\..*$",
    "",
    counts$Geneid
)

# Convert to numeric matrix
count_matrix <- as.matrix(count_matrix)

mode(count_matrix) <- "numeric"

# ------------------------------------------------------------
# Select significant genes
# ------------------------------------------------------------

sig_genes <- res$gene_id[
    !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) > 1
]

sig_genes <- intersect(
    sig_genes,
    rownames(count_matrix)
)

# Limit heatmap to top 24 DEGs
sig_genes <- head(
    sig_genes,
    24
)

heatmap_counts <- count_matrix[
    sig_genes,
    sample_names,
    drop = FALSE
]

# Log transformation
log_counts <- log2(
    heatmap_counts + 1
)

# Z-score each gene
z_counts <- t(
    scale(
        t(log_counts)
    )
)

# Sample annotation
annotation_col <- data.frame(
    Treatment = factor(
        c(
            "TBT",
            "TBT",
            "TBT",
            "DMSO",
            "DMSO",
            "DMSO"
        )
    )
)

rownames(annotation_col) <- sample_names

# ------------------------------------------------------------
# Generate heatmap
# ------------------------------------------------------------

pdf(
    "visualizations/Heatmap_DEGs_TBT_vs_DMSO.pdf",
    width = 9,
    height = 10
)

pheatmap(
    z_counts,
    scale = "none",
    annotation_col = annotation_col,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 7,
    main = "Significant DEGs: TBT vs DMSO"
)

dev.off()

png(
    "visualizations/Heatmap_DEGs_TBT_vs_DMSO.png",
    width = 2400,
    height = 2800,
    res = 300
)

pheatmap(
    z_counts,
    scale = "none",
    annotation_col = annotation_col,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 7,
    main = "Significant DEGs: TBT vs DMSO"
)

dev.off()

# ============================================================
# 10. PCA PLOT
# ============================================================

# Build DESeq2 object again
metadata <- data.frame(
    condition = factor(
        c(
            "TBT",
            "TBT",
            "TBT",
            "DMSO",
            "DMSO",
            "DMSO"
        )
    )
)

rownames(metadata) <- sample_names

dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = metadata,
    design = ~ condition
)

# Filter low-count genes
dds <- dds[
    rowSums(counts(dds) >= 10) >= 3,
]

dds <- DESeq(dds)

# VST transformation
vsd <- vst(
    dds,
    blind = FALSE
)

# PCA data
pca_data <- plotPCA(
    vsd,
    intgroup = "condition",
    returnData = TRUE
)

percent_var <- round(
    100 * attr(
        pca_data,
        "percentVar"
    )
)

# PCA plot
pca_plot <- ggplot(
    pca_data,
    aes(
        x = PC1,
        y = PC2,
        color = condition,
        label = name
    )
) +

geom_point(
    size = 4
) +

geom_text(
    vjust = -1,
    size = 3
) +

xlab(
    paste0(
        "PC1: ",
        percent_var[1],
        "% variance"
    )
) +

ylab(
    paste0(
        "PC2: ",
        percent_var[2],
        "% variance"
    )
) +

ggtitle(
    "PCA: TBT vs DMSO"
) +

theme_classic()

ggsave(
    "visualizations/PCA_TBT_vs_DMSO.png",
    pca_plot,
    width = 8,
    height = 6,
    dpi = 300
)

ggsave(
    "visualizations/PCA_TBT_vs_DMSO.pdf",
    pca_plot,
    width = 8,
    height = 6
)

# ============================================================
# 11. SAMPLE CORRELATION HEATMAP
# ============================================================

cor_matrix <- cor(
    assay(vsd),
    method = "pearson"
)

pdf(
    "visualizations/Sample_Correlation.pdf",
    width = 8,
    height = 7
)

pheatmap(
    cor_matrix,
    annotation_col = annotation_col,
    annotation_row = annotation_col,
    main = "Sample-to-Sample Correlation"
)

dev.off()

png(
    "visualizations/Sample_Correlation.png",
    width = 2400,
    height = 2100,
    res = 300
)

pheatmap(
    cor_matrix,
    annotation_col = annotation_col,
    annotation_row = annotation_col,
    main = "Sample-to-Sample Correlation"
)

dev.off()

# ============================================================
# 12. Save significant DEG table
# ============================================================

sig <- res[
    !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) > 1,
]

sig <- sig[
    order(sig$padj),
]

write.csv(
    sig,
    "visualizations/Significant_DEGs_final.csv",
    row.names = FALSE
)

# ============================================================
# Finished
# ============================================================

cat("\n============================================\n")
cat("RNA-seq visualization completed successfully\n")
cat("============================================\n")

cat(
    "\nSignificant DEGs:",
    nrow(sig),
    "\n"
)

cat(
    "Upregulated:",
    sum(sig$log2FoldChange > 0),
    "\n"
)

cat(
    "Downregulated:",
    sum(sig$log2FoldChange < 0),
    "\n"
)

cat(
    "\nResults saved in: visualizations/\n"
)
