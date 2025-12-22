# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(stringr)

# ------------------------------
# Load DEG list
# 
deg_genes <- read.table("/Users/tbehr/Desktop/SanRaffaele/Projects/snord116_Simone/Results/GOanalysis_03_11/Top100_foldEnrichment.txt", header = FALSE, stringsAsFactors = FALSE)
gene_symbols <- deg_genes$V1

# ------------------------------
# Convert gene symbols to Entrez IDs
gene_df <- bitr(gene_symbols, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

# Remove duplicates (if any)
gene_entrez <- unique(gene_df$ENTREZID)

# ------------------------------
# GO Enrichment Analysis (Biological Process)
go_bp <- enrichGO(gene = gene_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL",
                  keyType = "ENTREZID",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# ------------------------------
# KEGG Pathway Enrichment
kegg_res <- enrichKEGG(gene = gene_entrez,
                       organism = 'hsa',  # mouse KEGG code
                       keyType = "kegg",
                       pvalueCutoff = 0.05)

# Convert Entrez back to gene symbols for readability
kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# ------------------------------
# Visualizations

# GO Dotplot
dotplot(go_bp, showCategory = 20) + ggtitle(NULL) + 
  theme(text=element_text(size=24), axis.text.x = element_text(size=24), axis.text.y = element_text(size=24), axis.title.x = element_text(size=24)) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))

# KEGG Barplot
barplot(kegg_res, showCategory = 20) + ggtitle("Clusterprofiler: KEGG Pathway Enrichment - Top 100")

# Enrichment Map (optional)
emapplot(pairwise_termsim(go_bp))







