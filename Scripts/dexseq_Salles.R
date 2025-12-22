# (Copied and adapted from dexseq_diffExon (copied on 15/12/2025))

library(DEXSeq)
library(ggplot2)
library(dplyr)


N1_ctrl <- read.table('data/salles_counts_clean/N1_ctrl.dexeq_counts', sep = '\t')
N2_ctrl <- read.table('data/salles_counts_clean/N2_ctrl.dexeq_counts', sep = '\t')
N4_ctrl <- read.table('data/salles_counts_clean/N4_ctrl.dexeq_counts', sep = '\t')
N5_ctrl <- read.table('data/salles_counts_clean/N5_ctrl.dexeq_counts', sep = '\t')
N1_md <- read.table('data/salles_counts_clean/N1_md.dexeq_counts', sep = '\t')
N2_md <- read.table('data/salles_counts_clean/N2_md.dexeq_counts', sep = '\t')
N4_md <- read.table('data/salles_counts_clean/N4_md.dexeq_counts', sep = '\t')
N5_md <- read.table('data/salles_counts_clean/N5_md.dexeq_counts', sep = '\t')


count_files = c('data/salles_counts_clean/N1_ctrl.dexeq_counts',
                'data/salles_counts_clean/N2_ctrl.dexeq_counts',
                'data/salles_counts_clean/N4_ctrl.dexeq_counts',
                'data/salles_counts_clean/N5_ctrl.dexeq_counts',
                'data/salles_counts_clean/N1_md.dexeq_counts',
                'data/salles_counts_clean/N2_md.dexeq_counts',
                'data/salles_counts_clean/N4_md.dexeq_counts',
                'data/salles_counts_clean/N5_md.dexeq_counts'
)


data_list <- list(N1_ctrl, N2_ctrl, N4_ctrl, N5_ctrl, N1_md, N2_md, N4_md, N5_md)
common_features <- Reduce(intersect, lapply(data_list, function(df) df$V1))

count_matrix <- do.call(cbind, lapply(data_list, function(df) df$V2))
rownames(count_matrix) <- common_features
colnames(count_matrix) <- c('N1_ctrl', 'N2_ctrl', 'N4_ctrl', 'N5_ctrl', 'N1_md', 'N2_md', 'N4_md', 'N5_md')

sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = rep(c('CTRL', 'MD'), each = 4)
)


# dxd <- DEXSeqDataSet(
#   countData = count_matrix,
#   sampleData = sample_info,
#   design = ~ sample + exon + condition:exon,
#   featureID = rownames(count_matrix),
#   groupID = gsub(':.*', '', rownames(count_matrix)) # Extract transcript ID
# )

#reate DEXSeq object
#Needed to remove quotation marks in .dexseq_count files (thanks to Arthur (@93355568) on bioconductor forums for solution)
dxd <- DEXSeqDataSetFromHTSeq(
  countfiles = count_files,
  sampleData = sample_info,
  design = ~ sample + exon + condition:exon,
  flattenedfile = 'data/hsap_gencode_v31_dexseq.gff'
)


# Normalization
dxd <- estimateSizeFactors(dxd)

# Dispersion Estimation (of the negative binomial distribution)
dxd <- estimateDispersions(dxd)
#plotDispEsts(dxd)

dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")
results <- DEXSeqResults(dxd)  # have RData object saved at this step (results.RData)

###################
results <- readRDS('salles_dexseq_results_saved.RData')
###################

# count differentially expressed EXONS
table(results$padj < 0.05 & abs(results$log2fold_MD_CTRL) > 1 )
# count differentially affected genes
table(tapply(results$padj < 0.05, results$groupID, any))


# filtered_results <- results[!is.na(results$padj) & !is.infinite(results$log2fold_PW1_EDO), ]

# plot a specific gene
plotDEXSeq(results, "ENSG00000120742.11", legend=TRUE, displayTranscripts = T, cex.axis=1.2, cex=1.3, lwd=2)


results_df <- results %>% 
  as.data.frame %>%
  mutate(logp=-log10(padj)) %>% 
  mutate(highlight=logp > -log10(0.05) & (log2fold_MD_CTRL < -1 | log2fold_MD_CTRL > 1))

# Volcano Plot
ggplot(as.data.frame(results_df), aes(x=log2fold_MD_CTRL, y=logp, color=highlight)) +
  geom_point(size=2) +
  theme_bw() +
  ylab('-log10(adjusted p-value)') + xlab('log2(fold-change)') +
  theme(text=element_text(size=20))
#scale_color_manual(values=c('#7f7f7f','#cd1e05'))

# Volcano Plot with Enhanced Volcano
df_rownames <- rownames(results_df)
ensembl_ids <- sapply(strsplit(df_rownames, "\\."), `[`, 1)

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # change species if needed
mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

mapping_vect <- mapping$external_gene_name
names(mapping_vect) <- mapping$ensembl_gene_id

results_df$gene_name <- mapping_vect[ensembl_ids]

library(dplyr)
library(EnhancedVolcano)
# results_df_grouped <- results_df %>%
#   group_by(groupID) %>%
#   summarise(exonBaseMean = mean(exonBaseMean),
#             min_padj = min(padj, na.rm = T),
#             max_log2fold = max(log2fold_MD_CTRL, na.rm = T))
results_df_grouped <- results_df %>%
  group_by(groupID) %>%
  slice_min(padj, with_ties = FALSE) 

genes_to_label <- c('ANKRD11', 'H3C6','PRUNE2','NDUFS5','CBS','ZNF718','PWRN1','CCZ1','HINT1','TOR1AIP1','CTNNB1','PLD5P1')
EnhancedVolcano(results_df_grouped, x='log2fold_MD_CTRL', y='padj', lab=results_df_grouped$gene_name,
                selectLab = genes_to_label, boxedLabels = T,
                drawConnectors = T, widthConnectors = 1.5, typeConnectors = 'open', arrowheads = F,
                pCutoff = 0.05)


## Compare with OG Results

results_og <- readRDS('dexseq_results_saved.RData')
table(tapply(results_og$padj < 0.05, results_og$groupID, any))
results_df_og <- results_og %>% 
  as.data.frame %>%
  mutate(logp=-log10(padj)) %>% 
  mutate(highlight=logp > -log10(0.05) & (log2fold_PW1_CTRL < -1 | log2fold_PW1_CTRL > 1))

df_rownames_og <- rownames(results_df_og)
ensembl_ids_og <- sapply(strsplit(df_rownames_og, "\\."), `[`, 1)
mapping_og <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensembl_ids_og,
  mart = ensembl
)
mapping_vect_og <- mapping_og$external_gene_name
names(mapping_vect_og) <- mapping_og$ensembl_gene_id
results_df_og$gene_name <- mapping_vect_og[ensembl_ids_og]
results_df_grouped_og <- results_df_og %>%
  group_by(groupID) %>%
  slice_min(padj, with_ties = FALSE) 

sig_genes_og <- results_df_grouped_og$gene_name[results_df_grouped_og$padj < 0.1 & !is.na(results_df_grouped_og$padj)]
sig_genes_salles <- results_df_grouped$gene_name[results_df_grouped$padj < 0.1 & !is.na(results_df_grouped$padj)]

sig_genes_shared <- na.exclude(intersect(sig_genes_og, sig_genes_salles))
sig_genes_shared <- sig_genes_shared[sig_genes_shared!='']


DEXSeqHTML(results, FDR=0.1, color=c("#FF000080", "#0000FF80"))


sig_genes_shared_ensembl <- names(mapping_vect)[mapping_vect %in% sig_genes_shared]
sig_genes_shared_ensembl_v <- results_df_grouped$groupID[results_df_grouped$gene_name %in% sig_genes_shared]


