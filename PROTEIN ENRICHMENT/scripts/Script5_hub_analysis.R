library(igraph)
library(dplyr)

components_info <- components(
    clean_graph
)

largest_component <- which.max(
    components_info$csize
)

largest_nodes <- V(clean_graph)[
    components_info$membership ==
        largest_component
]

giant_graph <- induced_subgraph(
    clean_graph,
    vids = largest_nodes
)

giant_degree <- degree(
    giant_graph,
    mode = "all"
)

giant_betweenness <- betweenness(
    giant_graph,
    directed = FALSE,
    normalized = TRUE
)

giant_closeness <- closeness(
    giant_graph,
    normalized = TRUE
)

giant_hubs <- data.frame(
    STRING_id = names(giant_degree),
    degree = as.numeric(giant_degree),
    betweenness = as.numeric(giant_betweenness),
    closeness = as.numeric(giant_closeness)
)

giant_hubs <- giant_hubs %>%
    left_join(
        protein_annotation,
        by = "STRING_id"
    )

giant_hubs$DEG <- ifelse(
    giant_hubs$STRING_id %in% hits,
    "DEG",
    "First_shell"
)

giant_hubs$degree_norm <-
    giant_hubs$degree /
    max(
        giant_hubs$degree,
        na.rm = TRUE
    )

giant_hubs$betweenness_norm <-
    giant_hubs$betweenness /
    max(
        giant_hubs$betweenness,
        na.rm = TRUE
    )

giant_hubs$closeness_norm <-
    giant_hubs$closeness /
    max(
        giant_hubs$closeness,
        na.rm = TRUE
    )

giant_hubs$hub_score <- rowMeans(
    giant_hubs[
        ,
        c(
            "degree_norm",
            "betweenness_norm",
            "closeness_norm"
        )
    ],
    na.rm = TRUE
)

giant_hubs <- giant_hubs[
    order(
        giant_hubs$hub_score,
        decreasing = TRUE
    ),
]

giant_deg_hubs <- giant_hubs[
    giant_hubs$DEG == "DEG",
]

giant_deg_hubs <- giant_deg_hubs[
    order(
        giant_deg_hubs$hub_score,
        decreasing = TRUE
    ),
]

giant_deg_hubs_full <- giant_deg_hubs %>%
    left_join(
        final_string_map[
            ,
            c(
                "STRING_id",
                "SYMBOL",
                "baseMean",
                "log2FoldChange",
                "pvalue",
                "padj"
            )
        ],
        by = "STRING_id"
    )

write.csv(
    giant_hubs,
    "protein_enrichment/results/STRING_giant_component_hub_analysis.csv",
    row.names = FALSE
)

write.csv(
    giant_deg_hubs_full,
    "protein_enrichment/results/STRING_DEG_giant_component_hubs.csv",
    row.names = FALSE
)
