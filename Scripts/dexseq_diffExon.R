library(DEXSeq)
library(ggplot2)
library(dplyr)


EDO_1 <- read.table('data/DexSeq_counts_clean/EDO_1_clean.dexeq_counts', sep = '\t')
EDO_2 <- read.table('data/DexSeq_counts_clean/EDO_2_clean.dexeq_counts', sep = '\t')
EDO_3 <- read.table('data/DexSeq_counts_clean/EDO_3_clean.dexeq_counts', sep = '\t')
ND1_1 <- read.table('data/DexSeq_counts_clean/ND1_1_clean.dexeq_counts', sep = '\t')
ND1_2 <- read.table('data/DexSeq_counts_clean/ND1_2_clean.dexeq_counts', sep = '\t')
ND1_3 <- read.table('data/DexSeq_counts_clean/ND1_3_clean.dexeq_counts', sep = '\t')
PW1_1 <- read.table('data/DexSeq_counts_clean/PW1_1_clean.dexeq_counts', sep = '\t')
PW1_2 <- read.table('data/DexSeq_counts_clean/PW1_2_clean.dexeq_counts', sep = '\t')
PW1_3 <- read.table('data/DexSeq_counts_clean/PW1_3_clean.dexeq_counts', sep = '\t')


count_files = c('data/DexSeq_counts_clean/EDO_1_clean.dexeq_counts',
          'data/DexSeq_counts_clean/EDO_2_clean.dexeq_counts',
          'data/DexSeq_counts_clean/EDO_3_clean.dexeq_counts',
          'data/DexSeq_counts_clean/ND1_1_clean.dexeq_counts',
          'data/DexSeq_counts_clean/ND1_2_clean.dexeq_counts',
          'data/DexSeq_counts_clean/ND1_3_clean.dexeq_counts',
          'data/DexSeq_counts_clean/PW1_1_clean.dexeq_counts',
          'data/DexSeq_counts_clean/PW1_2_clean.dexeq_counts',
          'data/DexSeq_counts_clean/PW1_3_clean.dexeq_counts'
          )

data_list <- list(EDO_1, EDO_2, EDO_3, ND1_1, ND1_2, ND1_3, PW1_1, PW1_2, PW1_3)
common_features <- Reduce(intersect, lapply(data_list, function(df) df$V1))

count_matrix <- do.call(cbind, lapply(data_list, function(df) df$V2))
rownames(count_matrix) <- common_features
colnames(count_matrix) <- c('EDO_1', 'EDO_2', 'EDO_3', 'ND1_1', 'ND1_2', 'ND1_3', 'PW1_1', 'PW1_2', 'PW1_3')

sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = rep(c('CTRL', 'CTRL', 'PW1'), each = 3)
)

# 
# dxd <- DEXSeqDataSet(
#   countData = count_matrix,
#   sampleData = sample_info,
#   design = ~ sample + exon + condition:exon,
#   featureID = rownames(count_matrix),
#   groupID = gsub(':.*', '', rownames(count_matrix)) # Extract transcript ID
# )

# Create DEXSeq object
# Needed to remove quotation marks in .dexseq_count files (thanks to Arthur (@93355568) on bioconductor forums for solution)
dxd <- DEXSeqDataSetFromHTSeq(
  countfiles = count_files,  
  sampleData = sample_info, 
  design = ~ sample + exon + condition:exon,
  flattenedfile = 'data/gencode.v31.basic.annotation.gff'
)


# Normalization
dxd <- estimateSizeFactors(dxd)

# Dispersion Estimation (of the negative binomial distribution)
dxd <- estimateDispersions(dxd)
#plotDispEsts(dxd)


dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")
results <- DEXSeqResults(dxd)

# count differentially expressed EXONS
table(results$padj < 0.05 )
# count differentially affected genes
table(tapply(results$padj < 0.1, results$groupID, any))


# filtered_results <- results[!is.na(results$padj) & !is.infinite(results$log2fold_PW1_EDO), ]

# plot a specific gene
plotDEXSeq(results, "ENSG00000224078.15", legend=TRUE, displayTranscripts = T, cex.axis=1.2, cex=1.3, lwd=2)


results_df <- results %>% 
  as.data.frame %>%
  mutate(logp=-log10(padj)) %>% 
  mutate(highlight=logp > -log10(0.05) & (log2fold_PW1_CTRL < -1 | log2fold_PW1_CTRL > 1))

# Volcano Plot
ggplot(as.data.frame(results_df), aes(x=log2fold_PW1_CTRL, y=logp, color=highlight)) +
  geom_point(size=2) +
  theme_bw() +
  ylab('-log10(adjusted p-value)') + xlab('log2(fold-change)') +
  theme(text=element_text(size=20))
  #scale_color_manual(values=c('#7f7f7f','#cd1e05'))












