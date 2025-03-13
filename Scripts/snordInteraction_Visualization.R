library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

##########################################################
#     Plot snoGloBe interaction frequency across snord116
##########################################################

snord116_track <- GenomeAxisTrack(
  range = GRanges(seqnames = 'chr15',IRanges(start=25051476, end=25051572)), 
  genome = 'hg38',
  chromosome = 'chr15',
  name = 'SNORD116'
)

ase_ranges <- GRanges(seqnames = 'chr15', 
                                   ranges = IRanges(start = c(25051502, 25051546), end = c(25051518, 25051562)))
ase_track <- AnnotationTrack(
  range = ase_ranges,
  genome = 'hg38',
  chromosome = 'chr15',
  name = 'ASEs'
)

# get track for sequence of snord116
strack <- SequenceTrack(Hsapiens, chromosome = 'chr15')

# for plotting the full interaction windows  (WAY slower than it needs to be, should find a better way)
snoglobe_windows <- snoglobe_output[c('sno_window_start','sno_window_end')]
interaction_count <- numeric(length=nchar(SNORD116_SEQUENCE))
for(i in seq(1,97)){
  for(jrow in seq(1,nrow(snoglobe_windows))){
    if(i >= snoglobe_windows[jrow,1] & i < snoglobe_windows[jrow,2]){
      interaction_count[i] <- interaction_count[i] + 1
    }
  }
}

# for plotting the center of the interaction windows
snoglobe_means <- rowMeans(snoglobe_output[c('sno_window_start','sno_window_end')])-0.5
snoglobe_hist_values <- c(tabulate(snoglobe_means), rep(0,7))

# Plot line to show snoGloBe interaction frequency
snoglobe_interaction <- DataTrack(
  range = GRanges(seqnames = 'chr15',IRanges(start=25051476:25051572, end=25051476:25051572)),
  genome = 'hg38',
  chromosome = 'chr15',
  name = 'snoGloBe Interaction Frequency',
  data = interaction_count,
  type='l'
)

# annotate important regions of the snord116
highlight_regions <- GRanges(
  seqnames = 'chr15',
  ranges = IRanges(
    start = c(25051484, 25051502, 25051546),  # Start positions of highlights
    end = c(25051489, 25051518, 25051562)  # End positions of highlights
  ),
  id = c('C Box', 'ASE1', 'ASE2')  # Labels for the boxes
)

snord116_annotation <- AnnotationTrack(
  range = highlight_regions,
  genome = 'hg38',
  chromosome = 'chr15',
  name = 'Regions',
  col = 'red',
  featureAnnotation = 'id',
  fontcolor.feature = 'black',
  group = rep('snord116', 3)
)

ase_highlight <- HighlightTrack(
  trackList = list(snoglobe_interaction),
  range = ase_ranges
)

plotTracks(list(snord116_track, snord116_annotation, strack, snoglobe_interaction))



##########################################################
#     Plot gene exons + peak + snoGloBe interaction
##########################################################

mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

gene_symbol <- 'TP53'  # Replace with your gene of interest
genome <- 'hg38'       # Change if using another genome build

# Connect to Ensembl (for hg38, use 'GRCh38'; for hg19, use 'GRCh37')
mart <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', GRCh = 38)

# Retrieve gene information (exons, introns, TSS)
gene_info <- getBM(
  attributes = c('chromosome_name', 'exon_chrom_start', 'exon_chrom_end', 
                 'transcript_start', 'transcript_end', 'strand', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = gene_symbol,
  mart = mart
)

# Check if gene is found
if (nrow(gene_info) == 0) {
  stop('Gene not found in Ensembl. Check the gene symbol and genome assembly.')
}

#  Convert to GRanges for Gviz plotting
gene_ranges <- GRanges(
  seqnames = paste0('chr', gene_info$chromosome_name),
  ranges = IRanges(start = gene_info$exon_chrom_start, end = gene_info$exon_chrom_end),
  strand = ifelse(gene_info$strand == 1, '+', '-')
)

# Create GeneRegionTrack for plotting exons/introns
geneTrack <- GeneRegionTrack(
  gene_ranges,
  genome = genome,
  chromosome = unique(paste0('chr', gene_info$chromosome_name)),
  name = gene_symbol,
  showId = TRUE,       # Display transcript IDs
  transcriptAnnotation = 'symbol',  # Show gene symbol
  col = 'black',       # Border color
  fill = 'skyblue',    # Exon color
  stacking = 'full'
)

#  Get Transcription Start Site (TSS) and add a marker
tss_start <- min(gene_info$transcript_start)  # Take the first transcript TSS
tssTrack <- AnnotationTrack(
  start = tss_start, end = tss_start + 1,  # Small mark at TSS
  chromosome = unique(paste0('chr', gene_info$chromosome_name)),
  genome = genome,
  name = 'TSS',
  col = 'red', fill = 'red'
)

#  Plot the gene structure
plotTracks(
  list(geneTrack, tssTrack),
  from = min(gene_info$transcript_start) - 2000,  # Show upstream region
  to = max(gene_info$transcript_end) + 2000       # Show downstream region
)



















































