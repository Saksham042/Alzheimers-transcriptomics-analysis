## ===================================================================
## Alzheimer's Disease GEO Dataset — DEG, Pathway Enrichment
## & PPI Network Analysis
## Dataset: GSE118553 (GEO, Illumina HumanHT-12 V4.0)
##          Entorhinal cortex subset: 37 AD / 24 Control
##          (AsymAD samples and other brain regions excluded)
## ===================================================================
##
## NOTE ON REPRODUCIBILITY:
## The differential expression step (GEO2R group assignment + DEG
## export) was performed manually through NCBI's web-based GEO2R tool,
## selecting only entorhinal cortex AD and Control samples (excluding
## AsymAD and all other brain regions). The steps from DEG list
## onward are exactly what was run in this project.
## ===================================================================


## ---- Step 0: Install required packages (run once) -----------------
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db",
                        "enrichplot", "DOSE", "pathview"))


## ---- Step 1: Load & filter DEG list (from GEO2R) -------------------
## GEO2R export columns for this platform: ID, adj.P.Val, P.Value, t,
## B, logFC, GI, Gene.symbol, Gene.title
## Cutoff applied: adj.P.Val < 0.05 AND |logFC| > 1 -> 163 probes,
## deduplicated to 138 unique genes (kept lowest adj.P.Val probe per
## gene where multiple probes mapped to the same symbol).

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

deg <- read.csv("alzheimer_deg_list_full.csv")
head(deg)
str(deg)
nrow(deg)          # 138
table(deg$Direction)   # 23 Upregulated, 115 Downregulated


## ---- Step 2: Convert gene symbols to Entrez IDs --------------------
gene_list <- deg$Gene.symbol
gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
length(gene_list)   # 138

gene_conversion <- bitr(gene_list,
                         fromType = "SYMBOL",
                         toType   = "ENTREZID",
                         OrgDb    = org.Hs.eg.db)

nrow(gene_conversion)   # 131/138 mapped successfully (95%)
entrez_genes <- gene_conversion$ENTREZID


## ---- Step 3: KEGG pathway enrichment --------------------------------
kegg_result <- enrichKEGG(gene = entrez_genes,
                           organism = "hsa",
                           pvalueCutoff = 0.05)

nrow(as.data.frame(kegg_result))   # 5 significant pathways
head(as.data.frame(kegg_result))

## Check which genes are actually driving these KEGG hits
kegg_genes_entrez <- unique(unlist(strsplit(as.data.frame(kegg_result)$geneID, "/")))
kegg_gene_symbols <- bitr(kegg_genes_entrez,
                           fromType = "ENTREZID",
                           toType = "SYMBOL",
                           OrgDb = org.Hs.eg.db)
kegg_gene_symbols
## -> ITGB4, GJA1, ITGA6, DSG2, ITGB5, SSPN, SNTB1, CD36, CD44, FBLN1
##    (cell adhesion / cytoskeleton-ECM linkage genes; KEGG labels these
##    pathways by cardiac-muscle context, see report for interpretation)


## ---- Step 4: GO Biological Process enrichment -----------------------
go_result <- enrichGO(gene = entrez_genes,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

nrow(as.data.frame(go_result))   # 25 significant terms

go_simplified <- simplify(go_result, cutoff = 0.7,
                           by = "p.adjust", select_fun = min)
nrow(as.data.frame(go_simplified))   # 16 non-redundant terms

png("figures/go_dotplot.png", width = 861, height = 708)
dotplot(go_simplified, showCategory = 16)
dev.off()

png("figures/go_barplot.png", width = 861, height = 708)
barplot(go_simplified, showCategory = 16)
dev.off()


## ---- Step 5: Pathway-level visualization (fold-change overlay) -----
library(pathview)

fc_values <- deg$logFC
names(fc_values) <- gene_conversion$ENTREZID[match(deg$Gene.symbol,
                                                     gene_conversion$SYMBOL)]
fc_values <- fc_values[!is.na(names(fc_values))]

## ECM-receptor interaction - most biologically interpretable KEGG hit
pathview(gene.data = fc_values,
          pathway.id = "hsa04512",
          species = "hsa")
## Produces: hsa04512.pathview.png in the working directory


## ---- Step 6: Export gene list for STRING ---------------------------
## Take deg$Gene.symbol (138 genes) to string-db.org:
##   STRING -> "Multiple proteins" -> paste gene list ->
##   organism: Homo sapiens -> Search
## Export interactions as TSV ("string_interactions.tsv"),
## included in figures/ for reference.


## ---- Step 7: Hub gene identification (degree centrality) ------------
## Network was imported into Cytoscape from string_interactions.tsv,
## analyzed via Tools -> NetworkAnalyzer -> Analyze Network
## (undirected), styled by degree (node size + color, min 25-30/max
## 90-100 given wider degree range than Project 1), with a Prefuse
## Force Directed layout weighted by combined_score.
##
## Equivalent degree calculation done in R for reference:

interactions <- read.delim("figures/string_interactions.tsv",
                            stringsAsFactors = FALSE)

degree_table <- c(interactions$node1, interactions$node2) |>
  table() |>
  sort(decreasing = TRUE)

head(degree_table, 10)
## Top hubs in this project: CD44 (20), ERBB2 (15), CCL2 (14),
## VCAM1 (14), ARMC3 (13, hub of the separate ciliary-gene module)

## ===================================================================
## End of pipeline. See Alzheimers_DEG_Report.pdf for full write-up,
## discussion, and biological interpretation.
## ===================================================================
