library(STRINGdb)

ppi_enrichment <- string_db$get_ppi_enrichment(
    hits
)

ppi_result <- data.frame(
    enrichment = ppi_enrichment$enrichment,
    observed_edges = ppi_enrichment$edges,
    lambda = ppi_enrichment$lambda
)

print(ppi_result)

write.csv(
    ppi_result,
    "protein_enrichment/results/STRING_PPI_enrichment.csv",
    row.names = FALSE
)
