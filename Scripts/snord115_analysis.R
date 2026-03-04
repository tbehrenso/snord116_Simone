library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)


diff_peaks_broad <- read.csv('/Users/tbehr/Desktop/snord115_IGV/snord115 macs2 peaks/macs2_PE_narrow_q0.1/differential_peaks.narrowPeak', sep = '\t', header = F)

colnames(diff_peaks_broad) <- c('chr','start', 'end', 'peakID', 'score', 'strand', 'signalValue', 'pval', 'qval')

# convert to granges
diff_peaks_broad_granges <- makeGRangesFromDataFrame(diff_peaks_broad, keep.extra.columns = TRUE)


# ------------------
# Annotate with Genes
library(rtracklayer)

h38_gtf <- import('/Users/tbehr/Desktop/SanRaffaele/Projects/REFERENCE/Homo_sapiens.GRCh38.115.gtf')
h38_gtf <- h38_gtf[h38_gtf$type == 'gene']

if(!startsWith(as.character(h38_gtf@seqnames@values[1]),'chr')){
  seqlevels(h38_gtf) <- paste0("chr", seqlevels(h38_gtf))
}


annotation_matches <- findOverlaps(diff_peaks_broad_granges, h38_gtf, maxgap = 1000)

# diff_peaks_broad_granges$gene_id <- NA_character_
# diff_peaks_broad_granges$gene_id[queryHits(annotation_matches)] <- h38_gtf$gene_name[subjectHits(annotation_matches)]
# 
# diff_peaks_broad_granges$gene_type <- NA_character_
# diff_peaks_broad_granges$gene_type[queryHits(annotation_matches)] <- h38_gtf$gene_biotype[subjectHits(annotation_matches)]

h38_gene_id <- mcols(h38_gtf)$gene_name
h38_gene_type <- mcols(h38_gtf)$gene_biotype

diff_peaks_broad_granges$gene_id <- NA_character_
diff_peaks_broad_granges$gene_type <- NA_character_

# collapse per query
gene_id_collapsed <- tapply(
  h38_gene_id[subjectHits(annotation_matches)],
  queryHits(annotation_matches),
  function(x) paste(unique(x), collapse = ";")
)

gene_type_collapsed <- tapply(
  h38_gene_type[subjectHits(annotation_matches)],
  queryHits(annotation_matches),
  function(x) paste(unique(x), collapse = ";")
)

# assign using names (query indices)
diff_peaks_broad_granges$gene_id[as.integer(names(gene_id_collapsed))] <-
  gene_id_collapsed
diff_peaks_broad_granges$gene_type[as.integer(names(gene_type_collapsed))] <-
  gene_type_collapsed


diff_peaks_broad_granges$gene_id[grep('snoRNA', diff_peaks_broad_granges$gene_type)]
diff_peaks_broad_granges$gene_id[grep('snRNA', diff_peaks_broad_granges$gene_type)]


diff_peaks_df <- as.data.frame(diff_peaks_broad_granges)
write.csv(diff_peaks_df, '/Users/tbehr/Desktop/diff_peaks_broad.csv', quote = F)








