library(dplyr)
library(ggplot2)

final_results <- final_string_map %>%
    dplyr::select(
        SYMBOL,
        ENTREZID,
        STRING_id,
        baseMean,
        log2FoldChange,
        pvalue,
        padj
    ) %>%
    left_join(
        giant_hubs %>%
            dplyr::select(
                STRING_id,
                degree,
                betweenness,
                closeness,
                hub_score
            ),
        by = "STRING_id"
    )

final_results$Direction <- ifelse(
    final_results$log2FoldChange > 0,
    "Upregulated",
    "Downregulated"
)

final_results <- final_results[
    order(
        final_results$hub_score,
        decreasing = TRUE
    ),
]

candidate_hubs <- final_results[
    final_results$hub_score >=
        quantile(
            final_results$hub_score,
            0.75,
            na.rm = TRUE
        ),
]

write.csv(
    final_results,
    "protein_enrichment/results/FINAL_DEG_STRING_INTEGRATED_RESULTS.csv",
    row.names = FALSE
)

write.csv(
    candidate_hubs,
    "protein_enrichment/results/FINAL_CANDIDATE_HUB_PROTEINS.csv",
    row.names = FALSE
)
