library(STRINGdb)
library(igraph)

string_db <- STRINGdb$new(
    version = "12.0",
    species = 10090,
    score_threshold = 400
)

graph <- string_db$get_graph()

hits <- unique(
    final_string_map$STRING_id
)

graph_hits <- intersect(
    hits,
    V(graph)$name
)

ppi_graph <- induced_subgraph(
    graph,
    vids = graph_hits
)

write.csv(
    as_data_frame(
        ppi_graph,
        what = "edges"
    ),
    "protein_enrichment/results/STRING_DEG_PPI_edges.csv",
    row.names = FALSE
)
