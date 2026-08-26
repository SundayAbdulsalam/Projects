# Protein Enrichment and PPI Network Analysis of Differentially Expressed Genes

## Project Overview

This project performs protein-level functional enrichment and protein-protein interaction (PPI) network analysis of differentially expressed genes (DEGs) identified from an RNA-seq experiment investigating the molecular effects of preconception exposure to tributyltin (TBT).

The analysis integrates:

- RNA-seq differential expression
- Gene annotation
- STRING protein mapping
- Protein-protein interaction analysis
- PPI network expansion
- Connected-component analysis
- Network centrality analysis
- Hub protein identification
- STRING PPI enrichment
- GO Biological Process enrichment
- Reactome pathway enrichment
- DEG expression visualization
- Integration of differential expression with network centrality

---

# Biological Question

Which proteins and molecular pathways are associated with the transcriptional response to preconception TBT exposure, and which differentially expressed proteins occupy central positions within the associated protein interaction network?

---

# Experimental Comparison

The underlying RNA-seq experiment compares:

- TBT-exposed samples
- DMSO control samples

The differential-expression analysis identified:

**24 significant DEGs**

after DESeq2 analysis.

Of these, **20 DEGs** were successfully incorporated into the STRING protein-level analysis.

---

# Analysis Workflow

```text
RNA-seq
   |
   v
DESeq2 differential expression
   |
   v
Significant DEGs
   |
   v
Gene annotation
   |
   v
Entrez ID
   |
   v
STRING protein mapping
   |
   v
PPI network
   |
   +-----------------------+
   |                       |
   v                       v
Functional enrichment    Network analysis
   |                       |
   v                       v
GO / Reactome           Network expansion
                           |
                           v
                     Connected components
                           |
                           v
                     Giant component
                           |
                           v
                    Centrality analysis
                           |
                           v
                      Hub proteins
                           |
                           v
              DEG + network integration
```

### Software and Packages
R

### Recommended R version:

R >= 4.3

### Required packages:

``` bash
install.packages(c(
    "dplyr",
    "ggplot2",
    "igraph",
    "pheatmap"
))

Bioconductor:

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(
    c(
        "STRINGdb",
        "ComplexHeatmap"
    )
)
```
Load packages:

library(STRINGdb)
library(dplyr)
library(igraph)
library(ggplot2)
library(ComplexHeatmap)
Directory Setup

Create the project directories:
```
dir.create(
    "protein_enrichment",
    showWarnings = FALSE
)

dir.create(
    "protein_enrichment/input",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "protein_enrichment/results",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "protein_enrichment/figures",
    recursive = TRUE,
    showWarnings = FALSE
)
```
Step 1: Import Differential Expression Results

The analysis begins with the significant DEGs identified from the DESeq2 analysis.

Expected input:

DEG_TBT_vs_DMSO_significant.csv

Example R code:
```
deg_results <- read.csv(
    "protein_enrichment/input/DEG_TBT_vs_DMSO_significant.csv",
    stringsAsFactors = FALSE
)

head(deg_results)

dim(deg_results)
```
The project contains:

24 significant DEGs

Step 2: Prepare STRINGdb

STRINGdb was configured for:

Species: Mus musculus
Taxonomy ID: 10090
STRING version: 12.0

Initialize STRINGdb:

library(STRINGdb)
```
string_db <- STRINGdb$new(
    version = "12.0",
    species = 10090,
    score_threshold = 400,
    input_directory = ""
)
```
Check:
```
string_db$version
string_db$species
```
Expected:

[1] "12.0"

[1] 10090

Step 3: Map Genes to STRING Proteins

The Entrez IDs from the DEG annotation were mapped to STRING protein identifiers.
```
test_mapping <- string_db$map(
    data.frame(
        ENTREZID = final_protein_map$ENTREZID
    ),
    "ENTREZID",
    removeUnmappedRows = FALSE
)
```
Check mapping:
```
table(
    is.na(test_mapping$STRING_id)
)
```
The initial mapping identified:

10 directly mapped identifiers
10 identifiers requiring additional STRING annotation

The final protein mapping contained:
20 proteins

Step 4: Create Final STRING Protein Mapping

The mapped and previously annotated proteins were combined:
```
final_string_map <- rbind(
    mapped_final,
    unmapped_final
)
```
Check:
```
nrow(final_string_map)

sum(
    is.na(final_string_map$STRING_id)
)

length(
    unique(final_string_map$STRING_id)
)
```
Final result:
20 proteins
20 unique STRING IDs
0 missing STRING IDs

Save:
```
write.csv(
    final_string_map,
    "protein_enrichment/results/STRING_protein_mapping.csv",
    row.names = FALSE
)
```

Step 5: Retrieve STRING Protein Information
```
string_proteins <- string_db$get_proteins()

head(string_proteins)
```
Verify that the mapped proteins are present:
```
hits <- unique(
    final_string_map$STRING_id
)

sum(
    hits %in% string_proteins$protein_external_id
)
```
Expected:
20

Step 6: Obtain STRING Interaction Network
```
ppi <- string_db$get_interactions(
    hits
)
```
Because some mapped proteins were not represented as vertices in the STRING interaction graph, the direct interaction table may not contain all proteins.

Therefore the complete STRING graph was examined:
```
library(igraph)

graph <- string_db$get_graph()

vcount(graph)

ecount(graph)
```

The STRING graph contained:

10,640 vertices
34,626 interactions

### Step 7: Identify STRING Proteins Present in the Graph
```
graph_hits <- intersect(
    hits,
    V(graph)$name
)

length(graph_hits)

graph_hits
```
Result:
12 DEG proteins were represented as vertices in the STRING graph.

### Step 8: Construct the DEG PPI Subnetwork
```
ppi_graph <- induced_subgraph(
    graph,
    vids = graph_hits
)

vcount(ppi_graph)
ecount(ppi_graph)
```
The direct DEG-only network contained:
12 DEG proteins
0 direct edges

This indicates that the mapped DEGs did not form a sufficiently connected direct STRING subnetwork at the selected STRING graph threshold.
Therefore, first-shell network expansion was performed.

### Step 9: Expand the Network

The STRING network was expanded around the DEG proteins to identify interacting first-shell proteins.
The expanded network was then annotated with protein names:
```
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
```
Network annotation:
```
expanded_edges_annotated <- expanded_edges %>%
    left_join(
        protein_annotation,
        by = c("from" = "STRING_id")
    ) %>%
    rename(
        from_gene = preferred_name
    ) %>%
    left_join(
        protein_annotation,
        by = c("to" = "STRING_id")
    ) %>%
    rename(
        to_gene = preferred_name
    )
```
### Step 10: Remove Duplicate and Self Interactions
```
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
```
Save:

write.csv(
    expanded_edges_clean,
    "protein_enrichment/results/STRING_expanded_PPI_edges_clean.csv",
    row.names = FALSE
)
Step 11: Construct the Clean Network
clean_graph <- graph_from_data_frame(
    expanded_edges_clean[
        ,
        c("from", "to")
    ],
    directed = FALSE
)

Network size:

vcount(clean_graph)

ecount(clean_graph)

Final expanded network:

47 proteins
128 unique interactions
Step 12: Identify Connected Components
components_info <- components(
    clean_graph
)

table(
    components_info$csize
)

The network contained:

43-node component
2-node component
2-node component

Therefore, the largest connected component contained:

43 proteins
Step 13: Extract the Giant Component
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

Check:

vcount(giant_graph)

ecount(giant_graph)

This produces the principal 43-protein interaction network.

Step 14: Calculate Network Centrality

Degree:

giant_degree <- degree(
    giant_graph,
    mode = "all"
)

Betweenness:

giant_betweenness <- betweenness(
    giant_graph,
    directed = FALSE,
    normalized = TRUE
)

Closeness:

giant_closeness <- closeness(
    giant_graph,
    normalized = TRUE
)

Create hub table:

giant_hubs <- data.frame(
    STRING_id = names(giant_degree),
    degree = as.numeric(giant_degree),
    betweenness = as.numeric(giant_betweenness),
    closeness = as.numeric(giant_closeness)
)
Step 15: Annotate Network Proteins
giant_hubs <- giant_hubs %>%
    left_join(
        protein_annotation,
        by = "STRING_id"
    )

Identify DEGs:

giant_hubs$DEG <- ifelse(
    giant_hubs$STRING_id %in% hits,
    "DEG",
    "First_shell"
)
Step 16: Calculate Hub Score

Centrality measures were normalized:

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

Hub score:

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

Rank:

giant_hubs <- giant_hubs[
    order(
        giant_hubs$hub_score,
        decreasing = TRUE
    ),
]

Save:

write.csv(
    giant_hubs,
    "protein_enrichment/results/STRING_giant_component_hub_analysis.csv",
    row.names = FALSE
)
Step 17: Identify DEG Network Hubs
giant_deg_hubs <- giant_hubs[
    giant_hubs$DEG == "DEG",
]

giant_deg_hubs <- giant_deg_hubs[
    order(
        giant_deg_hubs$hub_score,
        decreasing = TRUE
    ),
]
Step 18: Integrate Network Results with DESeq2
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

Save:

write.csv(
    giant_deg_hubs_full,
    "protein_enrichment/results/STRING_DEG_giant_component_hubs.csv",
    row.names = FALSE
)
Step 19: PPI Enrichment

STRING PPI enrichment was calculated using:

ppi_enrichment <- string_db$get_ppi_enrichment(
    hits
)

The result was:

Observed interactions = 4
Expected interactions ≈ 1
Enrichment statistic = 0.0225
Lambda = 1

Save:

write.csv(
    data.frame(
        enrichment = ppi_enrichment$enrichment,
        edges = ppi_enrichment$edges,
        lambda = ppi_enrichment$lambda
    ),
    "protein_enrichment/results/STRING_PPI_enrichment.csv",
    row.names = FALSE
)
Step 20: Functional Enrichment
enrichment <- string_db$get_enrichment(
    hits
)

Significant terms:

enrichment_sig <- enrichment[
    enrichment$fdr < 0.05,
]

enrichment_sig <- enrichment_sig[
    order(
        enrichment_sig$fdr
    ),
]

Save:

write.csv(
    enrichment_sig,
    "protein_enrichment/results/STRING_significant_enrichment.csv",
    row.names = FALSE
)
Step 21: GO Biological Process
go_process <- enrichment_sig[
    enrichment_sig$category == "Process",
]

go_process <- go_process[
    order(
        go_process$fdr
    ),
]

Save:

write.csv(
    go_process,
    "protein_enrichment/results/STRING_GO_Biological_Process.csv",
    row.names = FALSE
)

Significant GO Biological Processes included:

Cellular zinc ion homeostasis
Cellular response to chemical stimulus
Step 22: Reactome Enrichment

STRING Reactome terms are identified using the RCTM category.

reactome <- enrichment_sig[
    enrichment_sig$category == "RCTM",
]

reactome <- reactome[
    order(
        reactome$fdr
    ),
]

Save:

write.csv(
    reactome,
    "protein_enrichment/results/STRING_Reactome_enrichment.csv",
    row.names = FALSE
)

Important enriched pathways included:

Phase I functionalization of compounds
Cytochrome P450-related pathways
CYP2E1 reactions
Fatty-acid-related processes
Cellular responses to stimuli
Metallothionein-associated processes
YAP1/WWTR1-associated transcriptional regulation
Step 23: GO Enrichment Plot
library(ggplot2)

ggplot(
    go_process,
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
        title = "STRING GO Biological Process Enrichment",
        x = "Number of Proteins",
        y = "Biological Process",
        color = "-log10(FDR)",
        size = "Protein Count"
    ) +
    theme_minimal()

Save:

ggsave(
    "protein_enrichment/figures/STRING_GO_enrichment_dotplot.pdf",
    width = 10,
    height = 6
)
Step 24: Reactome Plot
top_reactome <- head(
    reactome[
        order(reactome$fdr),
    ],
    10
)

ggplot(
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
        title = "STRING Reactome Pathway Enrichment",
        x = "Number of Proteins",
        y = "Reactome Pathway",
        color = "-log10(FDR)",
        size = "Protein Count"
    ) +
    theme_minimal()

Save:

ggsave(
    "protein_enrichment/figures/STRING_Reactome_enrichment_dotplot.pdf",
    width = 10,
    height = 7
)
Step 25: DEG Expression Heatmap

Load normalized counts:

normalized_counts <- read.csv(
    "protein_enrichment/input/normalized_counts_TBT_vs_DMSO.csv",
    row.names = 1,
    check.names = FALSE
)

Create DEG expression matrix:

deg_gene_ids <- final_string_map$SYMBOL

symbol_to_geneid <- setNames(
    final_protein_map$gene_id,
    final_protein_map$SYMBOL
)

deg_ensembl <- symbol_to_geneid[
    deg_gene_ids
]

deg_counts <- normalized_counts[
    rownames(normalized_counts) %in%
        deg_ensembl,
    ,
    drop = FALSE
]

Rename:

rownames(deg_counts) <-
    final_protein_map$SYMBOL[
        match(
            rownames(deg_counts),
            final_protein_map$gene_id
        )
    ]

Log transform:

log_counts <- log2(
    deg_counts + 1
)

Z-score:

heatmap_matrix <- t(
    scale(
        t(log_counts)
    )
)
Step 26: Heatmap
library(ComplexHeatmap)

sample_annotation <- data.frame(
    Treatment = c(
        "50TBT",
        "50TBT",
        "50TBT",
        "DMSO",
        "DMSO",
        "DMSO"
    )
)

rownames(sample_annotation) <-
    colnames(heatmap_matrix)

pdf(
    "protein_enrichment/figures/20_DEG_expression_heatmap.pdf",
    width = 9,
    height = 10
)

Heatmap(
    heatmap_matrix,
    name = "Z-score",
    top_annotation =
        HeatmapAnnotation(
            Treatment =
                sample_annotation$Treatment
        ),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = "TBT vs DMSO",
    row_title = "Significant DEGs"
)

dev.off()
Step 27: Integrate DEG Statistics and Network Centrality
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

Add expression direction:

final_results$Direction <- ifelse(
    final_results$log2FoldChange > 0,
    "Upregulated",
    "Downregulated"
)

Rank:

final_results <- final_results[
    order(
        final_results$hub_score,
        decreasing = TRUE
    ),
]

Save:

write.csv(
    final_results,
    "protein_enrichment/results/FINAL_DEG_STRING_INTEGRATED_RESULTS.csv",
    row.names = FALSE
)
Step 28: Candidate Hub Proteins

Candidate hubs were ranked according to network centrality.

The strongest DEG hubs identified were:

Rank	Protein	log2FC	FDR	Hub Score
1	Cdkn1a	-1.53	0.00475	0.590
2	Anpep	+1.01	0.0170	0.505
3	Mt1	-2.11	0.000238	0.347
4	Mtarc1	-1.88	0.0453	0.347
5	Ncoa1	-1.15	0.0310	0.238
6	Pdk4	-2.08	0.00475	0.218
7	Ank3	+2.91	0.0453	0.214

Save:

candidate_hubs <- final_results[
    final_results$hub_score >=
        quantile(
            final_results$hub_score,
            0.75,
            na.rm = TRUE
        ),
]

write.csv(
    candidate_hubs,
    "protein_enrichment/results/FINAL_CANDIDATE_HUB_PROTEINS.csv",
    row.names = FALSE
)
Biological Interpretation
1. Overall molecular response

The protein-level analysis suggests that preconception TBT exposure is associated with coordinated changes in proteins involved in:

chemical stimulus response
xenobiotic metabolism
metal and zinc homeostasis
transcriptional regulation
fatty-acid metabolism
cellular stress responses

The enrichment of these functional categories suggests that environmental chemical exposure may alter multiple interacting molecular systems rather than producing an isolated transcriptional response.

2. Chemical and xenobiotic response

The most prominent enrichment was:

Cellular response to chemical stimulus

This involved 11 of the mapped proteins.

Several proteins also occurred in:

Phase I - Functionalization of compounds
Cytochrome P450 pathways
CYP2E1 reactions

The presence of:

Cyp2f2
Cyp2d22
Ncoa1
Pdk4

among the enriched proteins suggests alterations in metabolic and xenobiotic-response processes.

These findings are consistent with activation or disruption of molecular systems involved in processing environmental chemical exposures.

However, enrichment alone does not demonstrate increased enzymatic activity.

3. Metal and zinc homeostasis

A significant biological theme was:

Cellular zinc ion homeostasis

The proteins contributing to this pathway included:

Mt1
Mt2
Slc39a8

Both Mt1 and Mt2 are metallothionein-associated proteins, while Slc39a8 is associated with metal ion transport.

This suggests that TBT exposure may influence cellular metal-handling mechanisms.

The strong differential expression of Mt1:

log2FC = -2.106
FDR = 0.000238

makes Mt1 an especially interesting candidate for further investigation.

4. Metallothionein biology

The enrichment analysis identified:

Metallothionein
Metal-thiolate cluster

and related metal-binding processes.

Mt1 and Mt2 were repeatedly represented in these terms.

This provides evidence for disruption of cellular metal-binding/homeostasis pathways.

Because metallothioneins participate in cellular responses to oxidative and chemical stress, this finding may indicate an altered stress-response state.

However, functional validation would be required to determine whether this represents impaired metal buffering, compensatory regulation, or another biological mechanism.

5. Transcriptional regulation

The analysis also identified pathways associated with:

TEA domain
YAP1/WWTR1 transcription
RUNX3/YAP1-mediated transcription

The major contributing proteins were:

Tead1
Tead4

TEAD transcription factors are important regulators of gene expression and cellular responses.

Their presence suggests that transcriptional regulatory networks may be altered following TBT exposure.

6. Network hub proteins
Cdkn1a

Cdkn1a was the strongest network-central DEG:

Degree = 16
Betweenness = 0.170
Closeness = 0.494
Hub score = 0.590
log2FC = -1.53

Its high degree indicates many network connections, while its relatively high betweenness indicates that it occupies an important position for communication between network regions.

Therefore, Cdkn1a represents the strongest candidate network hub in this analysis.

Anpep

Anpep had:

Degree = 7
Betweenness = 0.211
Hub score = 0.505
log2FC = +1.01

Although its degree was lower than Cdkn1a, its high betweenness suggests that Anpep may occupy an important bridging position in the network.

Mt1

Mt1 showed:

log2FC = -2.106
FDR = 0.000238

and was involved in multiple metal/zinc-related enrichment categories.

Its network degree was low, meaning that it should not be interpreted as a classical high-connectivity hub.

Instead, Mt1 is a biologically important DEG candidate because of its strong differential expression and repeated representation in functional enrichment.

Pdk4

Pdk4 was strongly downregulated:

log2FC = -2.078
FDR = 0.00475

and appeared in the chemical/metabolic response landscape.

It may therefore represent a candidate connecting metabolic regulation with the environmental-exposure response.

Important Interpretation Limitation

The current analysis is computational and observational.

Therefore, the results demonstrate:

association

rather than:

causation

The identified hub proteins should be described as:

candidate network hubs

rather than confirmed therapeutic targets or causal regulators.

Experimental validation would be required using approaches such as:

qPCR
Western blotting
targeted proteomics
immunohistochemistry
functional assays
pathway perturbation experiments
Main Conclusions

The integrated RNA-seq and protein-network analysis indicates that preconception TBT exposure is associated with coordinated changes in molecular pathways involving:

Chemical stimulus response
Xenobiotic metabolism
Cytochrome P450 activity
Zinc and metal homeostasis
Metallothionein biology
Transcriptional regulation
Cellular stress responses
Metabolic regulation

The PPI analysis identified a 43-protein giant component within the expanded 47-protein network.

Among the DEGs, Cdkn1a and Anpep showed the highest network centrality, while Mt1, Mtarc1, Pdk4, Ncoa1, and Ank3 also emerged as important candidate network-associated DEGs.

The combination of strong differential expression, network centrality, and functional enrichment identifies Cdkn1a, Anpep, Mt1, Mtarc1, Ncoa1, and Pdk4 as particularly interesting candidates for downstream biological validation.

Reproducibility

All analysis scripts are organized sequentially:
