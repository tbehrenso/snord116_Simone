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
  condition = rep(c('CTRL', 'CTRL', 'PW1'), each = 3),
  patient = rep(c('EDO', 'ND1', 'PW1'), each = 3)
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
  design = ~ sample + patient + exon + condition:exon,
  flattenedfile = 'data/gencode.v31.basic.annotation.gff'
)


# Normalization
dxd <- estimateSizeFactors(dxd)

# Dispersion Estimation (of the negative binomial distribution)
dxd <- estimateDispersions(dxd)
#plotDispEsts(dxd)


dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")
results <- DEXSeqResults(dxd)  # have RData object saved at this step (results.RData)

# count differentially expressed EXONS
table(results$padj < 0.05 & abs(results$log2fold_PW1_CTRL) > 1 )
# count differentially affected genes
table(tapply(results$padj < 0.05, results$groupID, any))


# filtered_results <- results[!is.na(results$padj) & !is.infinite(results$log2fold_PW1_EDO), ]

# plot a specific gene
plotDEXSeq(results, 'ENSG00000149970.16', legend=TRUE, displayTranscripts = T, cex.axis=1.2, cex=1.3, lwd=2)


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
#             max_log2fold = max(log2fold_PW1_CTRL, na.rm = T))
results_df_grouped <- results_df %>%
  group_by(groupID) %>%
  slice_min(padj, with_ties = FALSE) 

genes_to_label <- c('ANKRD11', 'H3C6','PRUNE2','NDUFS5','CBS','ZNF718','PWRN1','CCZ1','HINT1','TOR1AIP1','CTNNB1','PLD5P1')
EnhancedVolcano(results_df_grouped, x='log2fold_PW1_CTRL', y='padj', lab=results_df_grouped$gene_name,
                selectLab = genes_to_label, boxedLabels = T,
                drawConnectors = T, widthConnectors = 1.5, typeConnectors = 'open', arrowheads = F,
                pCutoff = 0.05)


# ------------------------------------------
#         5' / 3' Exon Usage Analysis
# ------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(cowplot) 

#dex <- read.delim("path_to_dexseq_results.tsv", stringsAsFactors=FALSE)
dex <- as.data.frame(results)
# required columns: gene_id, bin_id, log2FoldChange, padj, chr, start, end

dex$groupID <- vapply(strsplit(dex$groupID, "\\."), `[`, character(1), 1)
dex <- dplyr::rename(dex,
              gene_id = 'groupID',
              bin_id = 'featureID',
              log2FoldChange = 'log2fold_PW1_CTRL',
              chr = 'genomicData.seqnames',
              start = 'genomicData.start',
              end = 'genomicData.end')

# Compute bin genomic order and relative position within gene:
# We'll need gene-level coordinates from the GTF/TxDb
txdb <- makeTxDbFromGFF('/Users/tbehr/Desktop/Homo_sapiens.GRCh38.81.gtf', format="gtf") 
exons_by_gene <- exonsBy(txdb, by="gene")

# Make a helper df of gene start/end
gene_coords <- lapply(names(exons_by_gene), function(g){
  gr <- reduce(exons_by_gene[[g]])
  data.frame(gene_id=g,
             chr=as.character(seqnames(gr)[1]),
             gene_start=min(start(gr)),
             gene_end=max(end(gr)),
             strand=as.character(strand(gr)[1]),
             stringsAsFactors=FALSE)
}) %>% bind_rows()

# Join coordinates to dex
dex2 <- dex %>%
  left_join(gene_coords, by = c("gene_id"))

# If dex has start/end, compute bin_mid and order by position
dex2 <- dex2 %>%
  mutate(bin_mid = (start + end)/2) %>%
  group_by(gene_id) %>%
  arrange(if_else(strand=="+", bin_mid, -bin_mid)) %>% # order respecting strand
  mutate(bin_index = row_number(),
         n_bins = n(),
         rel_pos = (bin_index - 1) / (n_bins - 1) ) %>% # 0 = first bin, 1 = last bin
  ungroup()

####    Alternative to the block of code above (explicitly reorder, but seems unnecessary)
# dex2 <- dex2 %>%
#   group_by(gene_id) %>%
#   arrange(case_when(
#     strand == "+" ~ start,
#     strand == "-" ~ desc(start)
#   ), .by_group = TRUE) %>%
#   mutate(
#     bin_index = row_number(),
#     rel_pos = (bin_index - 1) / (n() - 1)
#   ) %>%
#   ungroup()



plot_gene_bins <- function(gene, df = dex2, sig_cut = 0.05, min_bins = 3){
  d <- df %>% filter(gene_id == gene)
  if(nrow(d) < min_bins) {
    message("Too few bins for gene: ", gene); return(NULL)
  }
  p <- ggplot(d, aes(x = rel_pos, y = log2FoldChange)) +
    geom_line(aes(group=gene_id), alpha=0.6) +
    geom_point(aes(color = padj < sig_cut), size=2) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    labs(title = paste0(gene, " (n_bins=", unique(d$n_bins), ")"),
         x = "Relative position (0 = 5' / 1 = 3')",
         y = "log2 fold-change (DEXSeq)") +
    theme_minimal()
  return(p)
}


# Create bins across relative position, e.g., 20 bins
dex_summary <- dex2 %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(pos_bin = cut(rel_pos, breaks = seq(0,1,length.out=21), include.lowest = TRUE, labels = FALSE)) %>%
  group_by(pos_bin) %>%
  summarise(mean_log2FC = mean(log2FoldChange, na.rm=TRUE),
            median_log2FC = median(log2FoldChange, na.rm=TRUE),
            n = n())

ggplot(dex_summary, aes(x = pos_bin, y = mean_log2FC)) +
  geom_point() +
  geom_smooth(method="loess", se=TRUE) +
  labs(x="Position bin across gene (0=5',20=3')", y="Mean log2FoldChange",
       title="Aggregate exon-bin log2FC across gene body") +
  theme_minimal()


# Aggregate only significant bins
# dex_summary_sig <- dex2 %>%
#   filter(!is.na(log2FoldChange), padj < 0.05) %>%
#   mutate(pos_bin = cut(rel_pos, breaks = seq(0,1,length.out=21),
#                        include.lowest = TRUE, labels = FALSE)) %>%
#   group_by(pos_bin) %>%
#   summarise(
#     mean_log2FC   = mean(log2FoldChange, na.rm=TRUE),
#     median_log2FC = median(log2FoldChange, na.rm=TRUE),
#     n = n()
#   )
# 
# ggplot(dex_summary_sig, aes(x = pos_bin, y = mean_log2FC)) +
#   geom_point(aes(size = n)) +
#   geom_smooth(method = "loess", se = TRUE) +
#   labs(
#     x = "Position bin across gene (0 = 5′, 20 = 3′)",
#     y = "Mean log2FoldChange (significant exons only)",
#     title = "Aggregate exon-bin log2FC (padj < 0.05)"
#   ) +
#   theme_minimal()

# Create with 95% CI
dex_summary_sig <- dex2 %>%
  filter(padj < 0.05, !is.na(log2FoldChange), log2FoldChange > log2(1.5)  | log2FoldChange < -log2(1.5)) %>%
  #filter(padj > 0.05, !is.na(log2FoldChange), log2FoldChange < log2(1.5)  & log2FoldChange > -log2(1.5)) %>%
  mutate(pos_bin = cut(rel_pos, breaks = seq(0,1,length.out=11),
                       include.lowest = TRUE, labels = FALSE)) %>%
  group_by(pos_bin) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    sd_log2FC   = sd(log2FoldChange, na.rm = TRUE),
    n           = n(),
    se_log2FC   = sd_log2FC / sqrt(n)
  )

dex_summary_sig <- dex_summary_sig %>%
  mutate(
    ci_low  = mean_log2FC - 1.96 * se_log2FC,
    ci_high = mean_log2FC + 1.96 * se_log2FC
  )

ggplot(dex_summary_sig,
       aes(x = pos_bin, y = mean_log2FC)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.2, fill = "skyblue") +
  geom_line(color = "steelblue") +
  geom_point(size = 2) +
  labs(x = "Position bin (0 = 5′, 20 = 3′)",
       y = "Mean log2FC ± 95 % CI",
       title = "Aggregate exon-bin log2FC with 95 % CI") +
  theme_minimal()

# Create 95% CI with bootstrapping
library(boot)
# Bootstrap mean for each position bin
boot_mean <- function(data, indices) mean(data[indices], na.rm=TRUE)
ci_boot <- function(x) {
  b <- boot(x, boot_mean, R=1000)
  quantile(b$t, c(0.025, 0.975), na.rm=TRUE)
}

dex_summary_boot <- dex2 %>%
  filter(padj < 0.05, !is.na(log2FoldChange), log2FoldChange > log2(1.5)  | log2FoldChange < -log2(1.5),
         abs(log2FoldChange) < 50) %>%
  #filter(padj > 0.05, !is.na(log2FoldChange), log2FoldChange < log2(1.5)  & log2FoldChange > -log2(1.5)) %>%
  mutate(pos_bin = cut(rel_pos, breaks = seq(0,1,length.out=11),
                       include.lowest = TRUE, labels = FALSE, right = T)) %>%
  group_by(pos_bin) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange, na.rm=TRUE),
    ci_low  = ci_boot(log2FoldChange)[1],
    ci_high = ci_boot(log2FoldChange)[2],
    n = n(),
    exonBaseMean = mean(exonBaseMean),
    median_log2FC = median(log2FoldChange, na.rm=T),
    padj = mean(padj),
    standDev = sd(log2FoldChange)
  )

ggplot(dex_summary_boot, aes(x = pos_bin, y = mean_log2FC)) +
  # geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
  #             alpha = 0.2, fill = "skyblue") +
  geom_ribbon(aes(ymin = mean_log2FC - standDev, ymax = mean_log2FC + standDev),
              alpha = 0.2, fill = "skyblue") +
  geom_line(color = "steelblue") +
  geom_point(aes(size = n)) +
  labs(x = "Position bin (0 = 5′, 20 = 3′)",
       y = "Mean log2FC",
       title = "Aggregate exon-bin log2FC with bootstrapped CI") +
  theme_minimal()




#      Statistical Test

# Classify bins as first / internal / last based on relative position
dex2 <- dex2 %>%
  mutate(
    bin_region = case_when(
      rel_pos <= 0.1 ~ "first",
      rel_pos >= 0.9 ~ "last",
      TRUE ~ "internal"
    ),
    sign = case_when(
      log2FoldChange > 0 ~ "up",
      log2FoldChange < 0 ~ "down",
      TRUE ~ "zero"
    )
  )

# Make contingency table for first vs last bins (ignore internal and zero)
tab <- dex2 %>%
  filter(sign != "zero", bin_region %in% c("first", "last")) %>%
  count(bin_region, sign) %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0)
# Convert cleanly to a matrix
tab_mat <- as.matrix(tab[, c("up", "down")])
rownames(tab_mat) <- tab$bin_region
# Print and test
print(tab_mat)
fisher.test(tab_mat)

# Spearman Correlation Test for General Downward Trend
cor_test <- cor.test(dex2$rel_pos, dex2$log2FoldChange,
                     method = "spearman", use = "complete.obs")
cor_test

# Linear Model
model <- lm(log2FoldChange ~ rel_pos, data = dex2)
summary(model)

ggplot(dex2, aes(rel_pos, log2FoldChange)) +
  geom_hex(bins = 50) +
  geom_smooth(method = "lm", color = "red") +
  labs(x = "Relative exon position (0=5′,1=3′)",
       y = "log₂ fold-change (DEXSeq)",
       title = "Global 5′→3′ trend in exon usage") +
  theme_minimal()

# Mixed Linear Model
library(lme4)
m <- lmer(log2FoldChange ~ rel_pos + (rel_pos | gene_id), data = dex2_filtered)
summary(m)
confint(m, parm = "rel_pos")




# ------------
# Filter dex2 to exclude genes without a significantly differential used exon

sig_bins <- dex2 %>%
  filter(!is.na(padj) & padj < 0.05, !is.na(log2FoldChange), log2FoldChange > log2(1.5)  | log2FoldChange < -log2(1.5))
genes_with_sig_bins <- unique(sig_bins$gene_id)
dex2_filtered <- dex2 %>%
  filter(gene_id %in% genes_with_sig_bins)



# ------------

library(dplyr)

dex2 <- dex2 %>%
  group_by(gene_id) %>%
  arrange(
    case_when(
      genomicData.strand == "+" ~ start,
      genomicData.strand == "-" ~ desc(start)
    ),
    .by_group = TRUE
  ) %>%
  mutate(
    bin_index = row_number(),
    rel_pos = (bin_index - 1) / (n() - 1)
  ) %>%
  ungroup()



# ----------------------
#     Plot by exon coordinates instead of relative exon position
n_bins <- 10

dex2 <- dex %>%
  left_join(gene_coords, by = c("gene_id"))

dex2_binned <- dex2 %>%
  group_by(gene_id) %>%
  mutate(
    gene_start = min(start),
    gene_end   = max(end),
    gene_len   = gene_end - gene_start,
    exon_mid   = (start + end) / 2,
    rel_coord  = ifelse(strand=='+', (exon_mid - gene_start) / gene_len, (gene_end - exon_mid) / gene_len),  # 0 = 5', 1 = 3'
    coord_bin  = floor(rel_coord * n_bins) + 1,
    n_exons = length(start),
    n = n()
  ) %>%
  ungroup()

dex2_binned <- dex2_binned %>%
  mutate(
    ci_low  = mean_log2FC - 1.96 * se_log2FC,
    ci_high = mean_log2FC + 1.96 * se_log2FC
  )

#     Combination of the relative exon position and the relative coordinates methods (to compare)
# dex2_binned <- dex2 %>%
#   mutate(bin_mid = (start + end)/2) %>%
#   group_by(gene_id) %>%
#   arrange(if_else(strand=="+", bin_mid, -bin_mid)) %>% # order respecting strand
#   mutate(
#     gene_start = min(start),
#     gene_end   = max(end),
#     gene_len   = gene_end - gene_start,
#     exon_mid   = (start + end) / 2,
#     rel_coord  = (exon_mid - gene_start) / gene_len,  # 0 = 5', 1 = 3'
#     coord_bin  = floor(rel_coord * n_bins) + 1,
#     n_exons = length(start),
#     bin_index = row_number(),
#     n_bins = n(),
#     rel_pos = (bin_index - 1) / (n_bins - 1)
#   ) %>%
#   ungroup()


dex2_binned <- dex2_binned %>%
  filter(padj < 0.05, !is.na(log2FoldChange), log2FoldChange > log2(1.5)  | log2FoldChange < -log2(1.5)) %>%
  group_by(coord_bin) %>%
  summarise(mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
            se_log2FC = sd(log2FoldChange, na.rm = TRUE) / sqrt(n()))



ggplot(dex2_summary, aes(x = coord_bin, y = mean_log2FC)) +
  geom_point() +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),alpha = 0.2, fill = "skyblue") +
  labs(x="Position bin across gene (0=5',20=3')", y="Mean log2FoldChange",title="Aggregate exon-bin log2FC across gene body") +
  theme_minimal()

  
# ------------------------------------
#     Separate long genes (gene locus > 100kb)
# ------------------------------------

# If dex has start/end, compute bin_mid and order by position
dex2 <- dex %>%
  left_join(gene_coords, by = c("gene_id")) %>%
  mutate(bin_mid = (start + end)/2) %>%
  group_by(gene_id) %>%
  arrange(if_else(strand == "+", bin_mid, -bin_mid)) %>%
  mutate(
    bin_index = row_number(),
    n_bins = n(),
    rel_pos = (bin_index - 1) / (n_bins - 1)
  ) %>%
  # compute total gene length per gene
  mutate(
    gene_start = min(start, na.rm = TRUE),
    gene_end   = max(end, na.rm = TRUE),
    gene_length = gene_end - gene_start + 1
  ) %>%
  ungroup()


library(boot)
# Bootstrap mean for each position bin
boot_mean <- function(data, indices) mean(data[indices], na.rm=TRUE)
ci_boot <- function(x) {
  b <- boot(x, boot_mean, R=1000)
  quantile(b$t, c(0.025, 0.975), na.rm=TRUE)
}

dex_summary_boot <- dex2 %>%
  filter(padj < 0.05, !is.na(log2FoldChange), 
         log2FoldChange > log2(1.5)  | log2FoldChange < -log2(1.5),
         gene_length > 100000) %>%
  #filter(padj > 0.05, !is.na(log2FoldChange), log2FoldChange < log2(1.5)  & log2FoldChange > -log2(1.5)) %>%
  mutate(pos_bin = cut(rel_pos, breaks = seq(0,1,length.out=11),
                       include.lowest = TRUE, labels = FALSE, right = T)) %>%
  group_by(pos_bin) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange, na.rm=TRUE),
    ci_low  = ci_boot(log2FoldChange)[1],
    ci_high = ci_boot(log2FoldChange)[2],
    n = n(),
    exonBaseMean = mean(exonBaseMean),
    median_log2FC = median(log2FoldChange, na.rm=T),
    padj = mean(padj)
  )

ggplot(dex_summary_boot, aes(x = pos_bin, y = mean_log2FC)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.2, fill = "skyblue") +
  geom_line(color = "steelblue") +
  geom_point(aes(size = n)) +
  labs(x = "Position bin (0 = 5′, 20 = 3′)",
       y = "Mean log2FC",
       title = "Aggregate exon-bin log2FC with bootstrapped CI") +
  theme_minimal()


# Mixed Linear Model
sig_bins <- dex2 %>%
  filter(!is.na(padj) & padj < 0.05, !is.na(log2FoldChange),
         log2FoldChange > log2(1.5)  | log2FoldChange < -log2(1.5),
         gene_length > 100000)
genes_with_sig_bins <- unique(sig_bins$gene_id)
dex2_filtered <- dex2 %>%
  filter(gene_id %in% genes_with_sig_bins)

library(lme4)
m <- lmer(log2FoldChange ~ rel_pos + (rel_pos | gene_id), data = dex2_filtered)
summary(m)
confint(m, parm = "rel_pos")


# ------------------------------------
#     Plot by raw gene location
# ------------------------------------

dex2 <- dex2 %>%
  mutate(exon_number = as.numeric(sub('.', '', dex2$bin_id)))

dex_summary_boot <- dex2 %>%
  filter(padj < 0.05, !is.na(log2FoldChange),
         log2FoldChange > log2(1.5)  | log2FoldChange < -log2(1.5)) %>%
  group_by(exon_number) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange, na.rm=TRUE),
    ci_low  = ci_boot(log2FoldChange)[1],
    ci_high = ci_boot(log2FoldChange)[2],
    n = n(),
    exonBaseMean = mean(exonBaseMean),
    median_log2FC = median(log2FoldChange, na.rm=T),
    padj = mean(padj)
  )

ggplot(dex_summary_boot, aes(x = exon_number, y = mean_log2FC)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.2, fill = "skyblue") +
  geom_line(color = "steelblue") +
  geom_point(aes(size = n)) +
  labs(x = "Position bin (0 = 5′, 20 = 3′)",
       y = "Mean log2FC",
       title = "Aggregate exon-bin log2FC with bootstrapped CI") +
  theme_minimal()


















