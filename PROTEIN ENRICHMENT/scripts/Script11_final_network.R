library(igraph)

V(clean_graph)$DEG <- ifelse(
    V(clean_graph)$name %in% hits,
    "DEG",
    "First_shell"
)

name_map <- setNames(
    protein_annotation$preferred_name,
    protein_annotation$STRING_id
)

V(clean_graph)$preferred_name <-
    name_map[V(clean_graph)$name]

V(clean_graph)$size <- ifelse(
    V(clean_graph)$DEG == "DEG",
    10,
    6
)

V(clean_graph)$label <- ifelse(
    V(clean_graph)$DEG == "DEG",
    V(clean_graph)$preferred_name,
    NA
)

pdf(
    "protein_enrichment/figures/STRING_DEG_first_shell_network.pdf",
    width = 12,
    height = 10
)

set.seed(123)

plot(
    clean_graph,
    layout = layout_with_fr(
        clean_graph
    ),
    vertex.size =
        V(clean_graph)$size,
    vertex.label =
        V(clean_graph)$label,
    vertex.label.cex = 0.7,
    vertex.label.color = "black",
    main =
        "STRING PPI Network of DEGs and First-Shell Proteins"
)

dev.off()
