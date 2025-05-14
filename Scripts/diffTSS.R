library(DESeq2)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(readxl)

tss_count_df <- read.csv('data/diffTSS/tss_notDeep_counts_paired_s2.txt', sep='\t', skip = 1)

#colnames(tss_count_df)[7:15] <- substr(colnames(tss_count_df)[7:15], 69, 73)
colnames(tss_count_df)[7:21] <- substr(colnames(tss_count_df)[7:21], 64, 68)


tss_count_matrix <- tss_count_df[,7:21]
rownames(tss_count_matrix) <- tss_count_df$Geneid


# Create metadata (sample info)
col_data <- data.frame(
  row.names = colnames(tss_count_df)[7:21],
  condition = c(rep('CTRL',6), rep('PW',9))
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = tss_count_matrix, colData = col_data, design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# View results
head(res)

res_df <- as.data.frame(res)

gene_transcript_split <- strsplit(rownames(res_df), split='_')

res_df <- data.frame(gene=sapply(gene_transcript_split, function(x) x[1]), transcript=sapply(gene_transcript_split, function(x) x[2]), baseMean=res_df$baseMean, log2FoldChange=res_df$log2FoldChange,
                     lfcSE=res_df$lfcSE, stat=res_df$stat, pvalue=res_df$pvalue, padj=res_df$padj)

res_df$gene <- res_df$gene %>% strsplit('\\.') %>% sapply('[', 1)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= res_df$gene, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
res_df$geneID <- geneIDs$SYMBOL[match(res_df$gene, geneIDs$GENEID)]

res_significant <- res_df[res_df$padj<0.05 & !is.na(res_df$padj),]


# Save as column to excel file
ase1_excel <- read_excel('data/ase1_to_peaks_complementarity_ALL.xlsx')

ase1_excel$deseq_notDeep_padj <- res_df$padj[match(ase1_excel$geneSymbol, res_df$geneID)]

write.csv(ase1_excel, '/Users/tbehr/Desktop/ase1_excel_updated.csv')
