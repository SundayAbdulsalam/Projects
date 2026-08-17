# ============================================================
# Annotate DEGs using local GRCm39.116 GTF
# TBT vs DMSO - Female GWAT
# ============================================================

# ------------------------------------------------------------
# 1. Read significant DEGs
# ------------------------------------------------------------

deg <- read.csv(
    "DEG_TBT_vs_DMSO_significant.csv",
    stringsAsFactors = FALSE
)

cat("Number of DEGs:", nrow(deg), "\n")

# ------------------------------------------------------------
# 2. Read GTF
# ------------------------------------------------------------

gtf_file <- "annotation/Mus_musculus.GRCm39.116.gtf"

cat("Reading GTF annotation...\n")

gtf <- readLines(gtf_file)

# Keep exon records because they contain gene annotations
gtf <- gtf[grepl("\texon\t", gtf)]

# ------------------------------------------------------------
# 3. Extract gene_id
# ------------------------------------------------------------

gene_id <- sub(
    '.*gene_id "([^"]+)".*',
    "\\1",
    gtf
)

# ------------------------------------------------------------
# 4. Extract gene_name
# ------------------------------------------------------------

gene_name <- sub(
    '.*gene_name "([^"]+)".*',
    "\\1",
    gtf
)

# ------------------------------------------------------------
# 5. Extract gene_biotype
# ------------------------------------------------------------

gene_biotype <- sub(
    '.*gene_biotype "([^"]+)".*',
    "\\1",
    gtf
)

# ------------------------------------------------------------
# 6. Create annotation table
# ------------------------------------------------------------

annotation <- data.frame(
    gene_id = gene_id,
    gene_symbol = gene_name,
    gene_biotype = gene_biotype,
    stringsAsFactors = FALSE
)

# Remove duplicated gene IDs
annotation <- annotation[
    !duplicated(annotation$gene_id),
]

# ------------------------------------------------------------
# 7. Match DEG IDs to GTF annotations
# ------------------------------------------------------------

deg_annotated <- merge(
    deg,
    annotation,
    by = "gene_id",
    all.x = TRUE
)

# Restore DEG order
deg_annotated <- deg_annotated[
    match(deg$gene_id, deg_annotated$gene_id),
]

# ------------------------------------------------------------
# 8. Add regulation status
# ------------------------------------------------------------

deg_annotated$significance <- ifelse(
    deg_annotated$log2FoldChange > 0,
    "Upregulated",
    "Downregulated"
)

# ------------------------------------------------------------
# 9. Create output directory
# ------------------------------------------------------------

dir.create(
    "annotation",
    showWarnings = FALSE
)

# ------------------------------------------------------------
# 10. Save annotated DEG table
# ------------------------------------------------------------

write.csv(
    deg_annotated,
    "annotation/DEG_TBT_vs_DMSO_annotated.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 11. Create simplified portfolio table
# ------------------------------------------------------------

portfolio <- deg_annotated[
    ,
    c(
        "gene_id",
        "gene_symbol",
        "gene_biotype",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj",
        "significance"
    )
]

write.csv(
    portfolio,
    "annotation/DEG_TBT_vs_DMSO_portfolio.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 12. Report results
# ------------------------------------------------------------

cat("\n============================================\n")
cat("DEG ANNOTATION COMPLETED\n")
cat("============================================\n")

cat(
    "\nSuccessfully annotated:",
    sum(!is.na(deg_annotated$gene_symbol) &
        deg_annotated$gene_symbol != ""),
    "/",
    nrow(deg_annotated),
    "genes\n"
)

cat("\nAnnotated genes:\n")

print(
    deg_annotated[
        ,
        c(
            "gene_id",
            "gene_symbol",
            "log2FoldChange",
            "padj",
            "significance"
        )
    ]
)

cat("\nFiles created:\n")
cat("annotation/DEG_TBT_vs_DMSO_annotated.csv\n")
cat("annotation/DEG_TBT_vs_DMSO_portfolio.csv\n")
