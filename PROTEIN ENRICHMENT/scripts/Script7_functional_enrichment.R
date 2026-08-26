enrichment <- string_db$get_enrichment(
    hits
)

enrichment_sig <- enrichment[
    enrichment$fdr < 0.05,
]

enrichment_sig <- enrichment_sig[
    order(
        enrichment_sig$fdr
    ),
]

go_process <- enrichment_sig[
    enrichment_sig$category == "Process",
]

reactome <- enrichment_sig[
    enrichment_sig$category == "RCTM",
]

go_process <- go_process[
    order(go_process$fdr),
]

reactome <- reactome[
    order(reactome$fdr),
]

write.csv(
    enrichment_sig,
    "protein_enrichment/results/STRING_significant_enrichment.csv",
    row.names = FALSE
)

write.csv(
    go_process,
    "protein_enrichment/results/STRING_GO_Biological_Process.csv",
    row.names = FALSE
)

write.csv(
    reactome,
    "protein_enrichment/results/STRING_Reactome_enrichment.csv",
    row.names = FALSE
)
