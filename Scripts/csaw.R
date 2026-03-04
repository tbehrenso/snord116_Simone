library(csaw)


# tf.data <- NFYAData()
# tf.data <- head(tf.data, -1) # skip the input.
# bam.files <- tf.data$Path

design_df <- data.frame(Sample = c('Ctrl','V5'),
                     Type = as.factor(c('Ctrl','V5')))

design <- model.matrix(~design_df$Type)
rownames(design) <- design_df$Sample
colnames(design) <- c('Intercept','Group')

# Load BAMs
bam.files <- c(
  '/Users/tbehr/Desktop/snord115_IGV/Ctrl_PE.dedup.bam',
  '/Users/tbehr/Desktop/snord115_IGV/V5_PE.dedup.bam'
)
param <- readParam(minq=20)
data <- windowCounts(bam.files, ext=53, width=15, param=param)


# Filter out uninteresting regions
binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
keep <- filterWindowsGlobal(data, binned)$filter > log2(4)
data <- data[keep,]

# Calc normalization factors
data <- normFactors(binned, se.out=data)

# Identify differentially bound windows
library(edgeR)
y <- asDGEList(data)
y <- estimateDisp(y, design)

fit.norep <- glmFit(y, design, dispersion=0.05)
results.norep <- glmLRT(fit.norep, contrast=c(0, 1))

# Correct for multiple testing
merged <- mergeResults(data, results.norep$table, tol=1000L)

# ------------------------------------------------------------------------------------------------------

merged <- readRDS('/Users/tbehr/Desktop/merged.RData')

results_best <- merged$best

results_best_sig <- results_best[results_best$FDR < 0.1,]
results_best_sig$grange <- merged$regions[results_best$FDR < 0.1,]

# -----------------------------------------------
# Annotate with Genes
library(rtracklayer)

h38_gtf <- import('/Users/tbehr/Desktop/SanRaffaele/Projects/REFERENCE/Homo_sapiens.GRCh38.115.gtf')
h38_gtf <- h38_gtf[h38_gtf$type == 'gene']

if(!startsWith(as.character(h38_gtf@seqnames@values[1]),'chr')){
  seqlevels(h38_gtf) <- paste0("chr", seqlevels(h38_gtf))
}

annotation_matches <- findOverlaps(results_best_sig$grange, h38_gtf, maxgap = 1000)


h38_gene_id <- mcols(h38_gtf)$gene_name
h38_gene_type <- mcols(h38_gtf)$gene_biotype

results_best_sig$gene_id <- NA_character_
results_best_sig$gene_type <- NA_character_

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
results_best_sig$gene_id[as.integer(names(gene_id_collapsed))] <-
  gene_id_collapsed
results_best_sig$gene_type[as.integer(names(gene_type_collapsed))] <-
  gene_type_collapsed


all_types <- character()
for(i in 1:nrow(results_best_sig)){
  id_vector <- strsplit(results_best_sig$gene_id[i], split=';')[[1]]
  type_vector <- strsplit(results_best_sig$gene_type[i], split=';')[[1]]
  
  all_types <- unique(c(all_types, type_vector))
  
}

# create columns for PRIMARY gene and genetype, for prioritization
nearest_hits <- distanceToNearest(results_best_sig$grange, h38_gtf)
nearest_gene_names <- mcols(h38_gtf)$gene_name[subjectHits(nearest_hits)]
nearest_gene_types <- mcols(h38_gtf)$gene_biotype[subjectHits(nearest_hits)]

results_best_sig$gene_primary <- NA_character_
results_best_sig$biotype_primary <- NA_character_
results_best_sig$gene_primary[queryHits(nearest_hits)] <- nearest_gene_names
results_best_sig$biotype_primary[queryHits(nearest_hits)] <- nearest_gene_types


results_best_sig$gene_id[grep('snoRNA', results_best_sig$gene_type)]

table(results_best_sig$gene_type)


# -----------------------------------------------
# Overlap with macs2 peaks

macs2_csaw_overlap <- findOverlaps(diff_peaks_broad_granges, results_best_sig$grange, maxgap = 50)





