library(dplyr)

network_summary <- data.frame(
    Total_DEGs = length(hits),

    STRING_mapped = length(
        unique(hits)
    ),

    STRING_graph_DEGs = length(
        graph_hits
    ),

    Expanded_network_nodes =
        vcount(clean_graph),

    Expanded_network_edges =
        ecount(clean_graph),

    Giant_component_nodes =
        vcount(giant_graph),

    Giant_component_edges =
        ecount(giant_graph),

    Giant_component_DEGs =
        sum(
            V(giant_graph)$name %in%
                hits
        )
)

write.csv(
    network_summary,
    "protein_enrichment/results/FINAL_NETWORK_SUMMARY.csv",
    row.names = FALSE
)

pathway_summary <- enrichment_sig[
    ,
    c(
        "category",
        "description",
        "number_of_genes",
        "p_value",
        "fdr",
        "preferredNames"
    )
]

write.csv(
    pathway_summary,
    "protein_enrichment/results/FINAL_PATHWAY_SUMMARY.csv",
    row.names = FALSE
)

cat(
    "\n====================================\n"
)

cat(
    "PROTEIN ENRICHMENT ANALYSIS COMPLETE\n"
)

cat(
    "====================================\n\n"
)

cat(
    "Total DEGs: ",
    length(hits),
    "\n"
)

cat(
    "STRING graph DEGs: ",
    length(graph_hits),
    "\n"
)

cat(
    "Expanded network nodes: ",
    vcount(clean_graph),
    "\n"
)

cat(
    "Expanded network edges: ",
    ecount(clean_graph),
    "\n"
)

cat(
    "Giant component nodes: ",
    vcount(giant_graph),
    "\n"
)

cat(
    "Significant enrichment terms: ",
    nrow(enrichment_sig),
    "\n"
)
