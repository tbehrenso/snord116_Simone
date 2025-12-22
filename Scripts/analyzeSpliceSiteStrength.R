library(DEXSeq)
library(ggplot2)
library(dplyr)

SE_MATS_JCEC <- read.csv('/Users/tbehr/Desktop/SE.MATS.JCEC.txt', sep = '\t')


se_mats_bed <- data.frame(chromosome = SE_MATS_JCEC$chr,
                          start = SE_MATS_JCEC$exonStart_0base,
                          end = SE_MATS_JCEC$exonEnd,
                          strand = SE_MATS_JCEC$strand)

write_tsv(se_mats_bed, '/Users/tbehr/Desktop/se_mats_jcec.bed')



# extract BED data from DEXseq results
results <- readRDS('salles_dexseq_results_saved.RData')

sig_results <- results[results$padj<0.05 & !is.na(results$padj) & 
                         abs(results$log2fold_MD_CTRL)>=1 & 
                         !is.na(results$log2fold_MD_CTRL),] 
nonsig_results <- results[!is.na(results$padj) & !is.na(results$log2fold_MD_CTRL) &
                         !(abs(results$log2fold_MD_CTRL)>=1 & results$padj<0.05)
                         ,]
extra_nonsig_results <- results[results$padj>0.05 & !is.na(results$padj) & 
                                  abs(results$log2fold_MD_CTRL)<1 & 
                                  !is.na(results$log2fold_MD_CTRL),] 

nonsig_results_126 <- extra_nonsig_results[sample(nrow(extra_nonsig_results), 1000), ]

if(F){
  sig_results <- nonsig_results_126
}


granges_data <- sig_results$genomicData
granges_data_frame <- as.data.frame(granges_data)
granges_data_frame <- granges_data_frame %>%
  mutate(gene_name = sapply(strsplit(rownames(granges_data_frame), split='\\.'), `[`, 1),
         exon_number = sapply(strsplit(rownames(granges_data_frame), split='\\:'), `[`, 2)
           )

ensembl_ids <- sapply(strsplit(sig_results$groupID, split = '\\.'), `[`, 1)
name_combo <- paste0(ensembl_ids, '|', sig_results$featureID, '|', 'donor')

# donor 9-mer
granges_data_frame_donor <- granges_data_frame %>%
  mutate(ninemer_start = ifelse(strand=='+', end - 3, start - 7),
         ninemer_end = ifelse(strand=='+', end + 6, start + 2))

dexseq_sig_bed_donor <- data.frame(chromosome = granges_data_frame_donor$seqnames,
                          start = granges_data_frame_donor$ninemer_start,
                          end = granges_data_frame_donor$ninemer_end,
                          name = name_combo,
                          score = 0,
                          strand = granges_data_frame_donor$strand)

# acceptor 23-mer
granges_data_frame_acceptor <- granges_data_frame %>%
  mutate(twentythreemer_start = ifelse(strand=='+', start - 21, end - 3),
         twentythreemer_end = ifelse(strand=='+', start + 2, end + 20))

dexseq_sig_bed_acceptor <- data.frame(chromosome = granges_data_frame_acceptor$seqnames,
                               start = granges_data_frame_acceptor$twentythreemer_start,
                               end = granges_data_frame_acceptor$twentythreemer_end,
                               name = name_combo,
                               score = 0,
                               strand = granges_data_frame_acceptor$strand)


# write bed file
write.table(dexseq_sig_bed_donor, '/Users/tbehr/Desktop/nonsig_1000_donor.bed', row.names = F, col.names = F, sep = '\t', quote = F)
write.table(dexseq_sig_bed_acceptor, '/Users/tbehr/Desktop/nonsig_1000_acceptor.bed', row.names = F, col.names = F, sep = '\t', quote=F)







