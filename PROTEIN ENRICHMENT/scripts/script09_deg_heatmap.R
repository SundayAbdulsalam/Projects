############################################################
# Script 09 — 09_deg_heatmap.R
# Purpose:
#   Generate a heatmap of normalized expression for the
#   significant DEGs used in the STRING protein analysis.
#
# Output:
#   protein_enrichment/figures/20_DEG_expression_heatmap.pdf
############################################################


############################
# 1. Load required packages
############################

library(pheatmap)
library(RColorBrewer)


############################
# 2. Create output directory
############################

dir.create(
    "protein_enrichment/figures",
    recursive = TRUE,
    showWarnings = FALSE
)


############################
# 3. Read normalized counts
############################

normalized_counts <- read.csv(
    "normalized_counts_TBT_vs_DMSO.csv",
    row.names = 1,
    check.names = FALSE
)

cat("Normalized count matrix dimensions:\n")
print(dim(normalized_counts))


############################
# 4. Read DEG results
############################

deg_results <- read.csv(
    "DEG_TBT_vs_DMSO_significant.csv",
    check.names = FALSE
)

cat("\nNumber of significant DEGs:\n")
print(nrow(deg_results))


############################
# 5. Identify gene column
############################

print(colnames(deg_results))


############################################################
# 6. Determine the gene identifiers used in the DEG results
############################################################

# The DEG table may contain SYMBOL or gene_id.
# Prefer SYMBOL when available.

if ("SYMBOL" %in% colnames(deg_results)) {

    deg_genes <- unique(
        deg_results$SYMBOL[
            !is.na(deg_results$SYMBOL) &
            deg_results$SYMBOL != ""
        ]
    )

} else if ("gene_id" %in% colnames(deg_results)) {

    deg_genes <- unique(
        deg_results$gene_id[
            !is.na(deg_results$gene_id) &
            deg_results$gene_id != ""
        ]
    )

} else {

    stop(
        "Could not find SYMBOL or gene_id in DEG results."
    )
}


cat("\nNumber of DEG identifiers found:\n")
print(length(deg_genes))


############################
# 7. Match DEGs to expression matrix
############################

expression_genes <- rownames(normalized_counts)

matched_genes <- intersect(
    deg_genes,
    expression_genes
)

cat("\nDEGs matched to normalized expression matrix:\n")
print(length(matched_genes))


if (length(matched_genes) == 0) {

    stop(
        "No DEG identifiers matched the normalized count matrix."
    )
}


############################
# 8. Extract DEG expression
############################

deg_expression <- normalized_counts[
    matched_genes,
    ,
    drop = FALSE
]


############################
# 9. Log2 transformation
############################

log_expression <- log2(
    deg_expression + 1
)


############################
# 10. Z-score normalization
############################

# Standardize each gene across samples.

heatmap_matrix <- t(
    scale(
        t(log_expression)
    )
)

heatmap_matrix <- as.matrix(
    heatmap_matrix
)


############################
# 11. Check for NA values
############################

heatmap_matrix[
    is.na(heatmap_matrix)
] <- 0


############################
# 12. Create sample metadata
############################

sample_names <- colnames(heatmap_matrix)

sample_group <- ifelse(
    grepl(
        "TBT|tbt|50TBT",
        sample_names,
        ignore.case = TRUE
    ),
    "50TBT",
    "DMSO"
)

annotation_col <- data.frame(
    Treatment = factor(
        sample_group,
        levels = c(
            "DMSO",
            "50TBT"
        )
    )
)

rownames(annotation_col) <- sample_names


############################
# 13. Order samples
############################

sample_order <- order(
    annotation_col$Treatment
)

heatmap_matrix <- heatmap_matrix[
    ,
    sample_order,
    drop = FALSE
]

annotation_col <- annotation_col[
    sample_order,
    ,
    drop = FALSE
]


############################
# 14. Order genes by variance
############################

gene_variance <- apply(
    heatmap_matrix,
    1,
    var,
    na.rm = TRUE
)

heatmap_matrix <- heatmap_matrix[
    order(
        gene_variance,
        decreasing = TRUE
    ),
    ,
    drop = FALSE
]


############################
# 15. Limit to 20 DEGs
############################

# Keep at most the 20 most variable significant DEGs.

if (nrow(heatmap_matrix) > 20) {

    heatmap_matrix <- heatmap_matrix[
        1:20,
        ,
        drop = FALSE
    ]

}


############################
# 16. Generate PDF
############################

pdf(
    "protein_enrichment/figures/20_DEG_expression_heatmap.pdf",
    width = 10,
    height = 12
)


############################
# 17. Generate heatmap
############################

pheatmap(
    heatmap_matrix,

    scale = "none",

    clustering_distance_rows = "euclidean",

    clustering_distance_cols = "euclidean",

    clustering_method = "complete",

    cluster_rows = TRUE,

    cluster_cols = TRUE,

    show_rownames = TRUE,

    show_colnames = TRUE,

    annotation_col = annotation_col,

    fontsize = 9,

    fontsize_row = 8,

    fontsize_col = 10,

    border_color = NA,

    main = "Expression Heatmap of Significant DEGs",

    angle_col = 45
)


############################
# 18. Close PDF
############################

dev.off()


############################
# 19. Verify output
############################

output_file <-
    "protein_enrichment/figures/20_DEG_expression_heatmap.pdf"

if (file.exists(output_file)) {

    cat(
        "\n============================================\n"
    )

    cat(
        "DEG HEATMAP COMPLETED SUCCESSFULLY\n"
    )

    cat(
        "============================================\n"
    )

    cat(
        "DEGs plotted:",
        nrow(heatmap_matrix),
        "\n"
    )

    cat(
        "Samples:",
        ncol(heatmap_matrix),
        "\n"
    )

    cat(
        "Output:\n",
        output_file,
        "\n"
    )

} else {

    stop(
        "Heatmap PDF was not created."
    )
}
