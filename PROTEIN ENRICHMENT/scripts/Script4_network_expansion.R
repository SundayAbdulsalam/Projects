expanded_edges_clean <- expanded_edges_annotated[
    !duplicated(
        expanded_edges_annotated[
            ,
            c("from", "to")
        ]
    ),
]

expanded_edges_clean <- expanded_edges_clean[
    expanded_edges_clean$from !=
        expanded_edges_clean$to,
]

clean_graph <- graph_from_data_frame(
    expanded_edges_clean[
        ,
        c("from", "to")
    ],
    directed = FALSE
)

write.csv(
    expanded_edges_clean,
    "protein_enrichment/results/STRING_expanded_PPI_edges_clean.csv",
    row.names = FALSE
)
