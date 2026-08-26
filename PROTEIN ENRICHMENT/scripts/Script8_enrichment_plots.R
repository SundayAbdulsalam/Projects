library(ggplot2)

top_go <- head(
    go_process,
    15
)

p1 <- ggplot(
    top_go,
    aes(
        x = number_of_genes,
        y = reorder(
            description,
            number_of_genes
        ),
        size = number_of_genes,
        color = -log10(fdr)
    )
) +
    geom_point() +
    labs(
        title =
            "STRING GO Biological Process Enrichment",
        x = "Number of Proteins",
        y = "Biological Process",
        color = "-log10(FDR)"
    ) +
    theme_minimal()

ggsave(
    "protein_enrichment/figures/STRING_GO_enrichment_dotplot.pdf",
    p1,
    width = 10,
    height = 7
)

top_reactome <- head(
    reactome,
    15
)

p2 <- ggplot(
    top_reactome,
    aes(
        x = number_of_genes,
        y = reorder(
            description,
            number_of_genes
        ),
        size = number_of_genes,
        color = -log10(fdr)
    )
) +
    geom_point() +
    labs(
        title =
            "STRING Reactome Pathway Enrichment",
        x = "Number of Proteins",
        y = "Reactome Pathway",
        color = "-log10(FDR)"
    ) +
    theme_minimal()

ggsave(
    "protein_enrichment/figures/STRING_Reactome_enrichment_dotplot.pdf",
    p2,
    width = 10,
    height = 7
)
