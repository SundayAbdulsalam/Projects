library(STRINGdb)
library(dplyr)

string_db <- STRINGdb$new(
    version = "12.0",
    species = 10090,
    score_threshold = 400,
    input_directory = ""
)

test_mapping <- string_db$map(
    data.frame(
        ENTREZID = final_protein_map$ENTREZID
    ),
    "ENTREZID",
    removeUnmappedRows = FALSE
)

string_proteins <- string_db$get_proteins()

protein_annotation <- string_proteins[
    ,
    c(
        "protein_external_id",
        "preferred_name"
    )
]

colnames(protein_annotation) <- c(
    "STRING_id",
    "preferred_name"
)

hits <- unique(
    final_string_map$STRING_id
)

write.csv(
    final_string_map,
    "protein_enrichment/results/STRING_protein_mapping.csv",
    row.names = FALSE
)
