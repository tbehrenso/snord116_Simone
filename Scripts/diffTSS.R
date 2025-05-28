library(DESeq2)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(readxl)
library(biomaRt)

tss_count_df <- read.csv('data/diffTSS/tss_counts_firstexon_extended_s2.txt', sep='\t', skip = 1)

colnames(tss_count_df)[7:15] <- substr(colnames(tss_count_df)[7:15], 69, 73)
#colnames(tss_count_df)[7:21] <- substr(colnames(tss_count_df)[7:21], 64, 68)


tss_count_matrix <- tss_count_df[,7:15]
rownames(tss_count_matrix) <- tss_count_df$Geneid


# Create metadata (sample info)
col_data <- data.frame(
  row.names = colnames(tss_count_df)[7:15],
  condition = c(rep('CTRL',3), rep('PW',6))
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
res_df$geneID <- geneIDs$SYMBOL[match(res_df$gene, geneIDs$GENEID, nomatch=NA)]

res_significant <- res_df[res_df$padj<0.05 & !is.na(res_df$padj),]


# Save as column to excel file
ase1_excel <- read_excel('data/ase1_to_peaks_complementarity_ALL.xlsx')

ase1_excel$deseq_padj_new <- res_df$padj[match(ase1_excel$geneSymbol, res_df$geneID)]

ase1_excel$deseq_padj_new <- unlist(lapply(ase1_excel$geneSymbol, function(x) min(res_df$padj[which(res_df$geneID %in% x)], na.rm=T)), use.names = F)

write.csv(ase1_excel, '/Users/tbehr/Desktop/ase1_excel_updated.csv')


# ------------------------------------
#   Differential TSS per Gene
# ------------------------------------

gene_transcript_split <- strsplit(rownames(tss_count_matrix), split='_')
matrix_genes <- sapply(gene_transcript_split, function(x) x[1])
matrix_transcripts <- sapply(gene_transcript_split, function(x) x[2])

matrix_genes <- matrix_genes %>% strsplit('\\.') %>% sapply('[', 1)

# Convert ensembl to gene symbol
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = matrix_genes,
  mart = ensembl
)

matrix_geneSymbols <- gene_info$hgnc_symbol[match(matrix_genes, gene_info$ensembl_gene_id)]

gene_pval_df <- data.frame(gene = ase1_excel$geneSymbol, min_pval = numeric(length=length(ase1_excel$geneSymbol)))

# for(gene in ase1_excel$geneSymbol){
#   count_matrix_subset <- tss_count_matrix[matrix_geneSymbols == gene,]
#   
#   if(dim(count_matrix_subset)[1] > 1 & 
#      sum(count_matrix_subset) != 0 & 
#      sum(!apply(count_matrix_subset, MARGIN = 1, function(x) any(x==0)))>2){
#     
#     #count_matrix_subset <- count_matrix_subset[rowSums(count_matrix_subset)!=0,]
#     # Create DESeq2 dataset
#     gene_dds <- DESeqDataSetFromMatrix(countData = count_matrix_subset, colData = col_data, design = ~ condition)
#     
#     # Run DESeq2
#     gene_dds <- DESeq(gene_dds, fitType = 'mean')   # added fitType='mean' here to avoid an error. Maybe better to try removing all lines with all zeros?
#     gene_res <- results(gene_dds)
#     
#     lowest_pval <- min(gene_res$padj, na.rm=T)
#     
#     gene_pval_df$min_pval[gene_pval_df$gene==gene] <- lowest_pval
#   } else {
#     gene_pval_df$min_pval[gene_pval_df$gene==gene] <- NA
#   }
# }

# two-way ANOVA method
for(gene in ase1_excel$geneSymbol){
  
  count_matrix_subset <- tss_count_matrix[matrix_geneSymbols == gene,]
  count_matrix_subset <- count_matrix_subset[rowSums(count_matrix_subset)!=0,]
  
  if(dim(count_matrix_subset)[1] > 1){
    count_matrix_subset$geneID <- rownames(count_matrix_subset)
    
    count_matrix_long <- as.data.frame(pivot_longer(count_matrix_subset, colnames(count_matrix_subset)[-10]))
    
    count_matrix_long$geneID <- factor(count_matrix_long$geneID)
    count_matrix_long$name <- factor(count_matrix_long$name)
    count_matrix_long$condition <- ifelse(count_matrix_long$name %in% c('PW1_1','PW1_2','PW1_3'), 'PW', 'CTRL')
    count_matrix_long$condition <- factor(count_matrix_long$condition)
    
    model <- lm(value ~ condition * geneID, data = count_matrix_long)
    
    count_anova <- anova(model)
    
    gene_pval_df$min_pval[gene_pval_df$gene==gene] <- count_anova$`Pr(>F)`[3]
    
  } else {
    gene_pval_df$min_pval[gene_pval_df$gene==gene] <- NA
  }
}

# Save as column to excel file
ase1_excel <- read_excel('data/ase1_to_peaks_complementarity_ALL.xlsx')

ase1_excel$tss_anova_pval <- gene_pval_df$min_pval

write.csv(ase1_excel, '/Users/tbehr/Desktop/ase1_excel_updated.csv')








