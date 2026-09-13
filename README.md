# Alzheimer's Disease Gene Expression Analysis (GSE118553)

Differential gene expression, pathway enrichment, and protein-protein
interaction network analysis of entorhinal cortex tissue from a public
Alzheimer's disease transcriptomics dataset.

## Dataset
- **Source:** NCBI GEO, accession [GSE118553](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118553)
- **Platform:** Illumina HumanHT-12 V4.0 expression beadchip
- **Samples used:** Entorhinal cortex subset only — 37 AD, 24 control (asymptomatic AD samples and other brain regions excluded to keep a clean two-group comparison)

## Workflow
1. **Differential expression** — GEO2R (limma), filtered to adj. p-value < 0.05 and |logFC| > 1 → 163 significant probes → **138 unique DEGs** after deduplication (23 upregulated, 115 downregulated)
2. **Functional enrichment** — `clusterProfiler`: KEGG (5 pathways) and GO Biological Process (25 terms, simplified to 16)
3. **Pathway visualization** — `pathview`, logFC mapped onto the KEGG ECM-receptor interaction pathway (hsa04512)
4. **PPI network** — STRING-db → Cytoscape, node degree analysis to identify hub genes

## Key finding
**CD44** was the dominant network hub (degree = 20), alongside **ERBB2, CCL2, VCAM1**. All five KEGG-hit genes (ITGA6, ITGB4, ITGB5, CD44, CD36) were downregulated and cluster around laminin-binding integrins — consistent with blood-brain barrier / basement membrane integrity loss, a documented feature of AD. CCL2 and VCAM1 point to neuroinflammatory involvement. A structurally separate ciliary-gene module (ARMC3, RSPH4A, CFAP43/53) was also identified. Full discussion, including caveats, in the report.

## Repo contents

| File | Description |
|---|---|
| `Alzheimers_DEG_Report.pdf` | Full write-up: methods, results, figures, discussion, limitations |
| `Alzheimers_DEG_Report.docx` | Same report, editable Word format |
| `analysis_pipeline.R` | Full R pipeline: gene ID conversion → KEGG/GO enrichment → pathview → hub gene calc |
| `alzheimer_deg_list_full.csv` | Full DEG table — 138 genes with logFC, adj.P.Val, direction, etc. |
| `figures/go_dotplot.png` | Top 16 non-redundant enriched GO terms (dot plot) |
| `figures/go_barplot.png` | Same 16 GO terms as a bar plot |
| `figures/pathway_diagram_hsa04512.png` | KEGG ECM-receptor interaction pathway with logFC overlay |
| `figures/ppi_network_cytoscape.png` | STRING/Cytoscape interaction network, styled by hub degree |
| `figures/string_interactions.tsv` | Raw STRING interaction edge list |

## Note on group definition
The full GSE118553 series profiles **4 brain regions** (cerebellum, entorhinal cortex, temporal cortex, frontal cortex) and includes an **asymptomatic AD (AsymAD)** category in addition to AD and control. This analysis uses **only the entorhinal cortex subset**, and **excludes AsymAD samples entirely**, to keep a clean, clinically-defined two-group comparison rather than mixing brain regions or diluting the contrast with an intermediate disease state. This mirrors the approach taken in published re-analyses of this dataset.

## Tools used
R (clusterProfiler, org.Hs.eg.db, enrichplot, pathview) · GEO2R · STRING-db · Cytoscape

---
*Companion project: [asthma-transcriptomics-analysis](https://github.com/Saksham042/asthma-transcriptomics-analysis) — same pipeline applied to an asthma GEO dataset.*
