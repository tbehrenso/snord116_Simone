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

# Connect to Ensembl (for hg38, use 'GRCh38'; for hg19, use 'GRCh37')
mart <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl')

gene_symbol <- 'NCL'  
genome <- 'hg38'     


# get attributes in two separate calls (Biomart does not allow calling attributes from different "pages" in the same call)
# (thanks James W. MacDonald)
gene_info_a <- getBM(
  attributes = c('chromosome_name', 'ensembl_transcript_id',
                 'strand', 'hgnc_symbol'),
  filters = 'hgnc_symbol',
  values = gene_symbol,
  mart = mart)

gene_info_b <- getBM(
  attributes = c('ensembl_transcript_id', '3_utr_start','3_utr_end','transcript_start','transcript_end', 'ensembl_exon_id', 'exon_chrom_start', 'exon_chrom_end'),
  filters = 'ensembl_transcript_id',
  values = as.character(gene_info_a[,2]),
  mart = mart)

gene_info <- merge(gene_info_a, gene_info_b)


#  Convert to GRanges for Gviz plotting
gene_ranges <- GRanges(
  seqnames = paste0('chr', gene_info$chromosome_name),
  ranges = IRanges(start = gene_info$exon_chrom_start, end = gene_info$exon_chrom_end),
  strand = ifelse(gene_info$strand == 1, '+', '-'),
  transcript = gene_info$ensembl_transcript_id,
  exonID = gene_info$ensembl_exon_id
)

gene_chromosome <- unique(paste0('chr', gene_info$chromosome_name))

# Create GeneRegionTrack for plotting exons/introns
geneTrack <- GeneRegionTrack(
  gene_ranges,
  genome = genome,
  chromosome = gene_chromosome,
  name = gene_symbol,
  #showId = TRUE,
  #transcriptAnnotation = 'symbol',
  col = 'black',       # Border color
  fill = 'skyblue',    # Exon color
  stacking = 'full',
  transcript = gene_info$ensembl_transcript_id,
  exon = gene_info$ensembl_exon_id
)


tss_start <- min(gene_info$transcript_start)  # Take the first transcript TSS
tssTrack <- AnnotationTrack(
  start = unique(gene_info$transcript_start), end = unique(gene_info$transcript_start)+1,  # Small mark at TSS
  chromosome = unique(paste0('chr', gene_info$chromosome_name)),
  genome = genome,
  name = 'TSS',
  col = 'red', fill = 'red'
)

ideogram_track <- IdeogramTrack(genome = genome, chromosome = gene_chromosome)

# Annotation Track for Peak Location
SYS5_differential <- read_excel('peaksAnno.xlsx', sheet='SYS5_differential')
HEK_differential <- read_excel('peaksAnno.xlsx', sheet='differential_peaks')

if(gene_symbol %in% SYS5_differential$SYMBOL){
  
  gene_SYS5_peaks <- SYS5_differential[SYS5_differential$SYMBOL == gene_symbol & !is.na(SYS5_differential$SYMBOL),]
  
  SYS5_peaks <- GRanges(
    seqnames = gene_chromosome,
    ranges = IRanges(
      start = gene_SYS5_peaks$start,  
      end = gene_SYS5_peaks$end  
    )
  )
  SYS5_annotation <- AnnotationTrack(
    range = SYS5_peaks,
    genome = genome,
    chromosome = gene_chromosome,
    name = 'SYS5 Peaks',
    col = 'red',
    fill = 'orange'
  )
} else {
  gene_SYS5_peaks <- NA
  SYS5_annotation <- NA
  print(paste0(gene_symbol, ' not in SYS5 peaks'))
}


if(gene_symbol %in% HEK_differential$SYMBOL){
  
  gene_HEK_peaks <- HEK_differential[HEK_differential$SYMBOL == gene_symbol & !is.na(HEK_differential$SYMBOL),]
  
  HEK_peaks <- GRanges(
    seqnames = gene_chromosome,
    ranges = IRanges(
      start = gene_HEK_peaks$start,  
      end = gene_HEK_peaks$end  
    )
  )
  HEK_annotation <- AnnotationTrack(
    range = HEK_peaks,
    genome = genome,
    chromosome = gene_chromosome,
    name = 'HEK Peaks',
    col = 'red',
    fill = '#EB5017'
  )
} else {
  gene_HEK_peaks <- NA
  HEK_annotation <- NA
  print(paste0(gene_symbol, ' not in HEK peaks'))
}


# Annotation Track for snoGloBe interaction
all_starts <- c(gene_info$transcript_start, if(any(!is.na(gene_SYS5_peaks))){gene_SYS5_peaks$start}, if(any(!is.na(gene_HEK_peaks))){gene_HEK_peaks$start})
all_ends <- c(gene_info$transcript_end, if(any(!is.na(gene_SYS5_peaks))){gene_SYS5_peaks$end}, if(any(!is.na(gene_HEK_peaks))){gene_HEK_peaks$end})

plot_range_start <- min(all_starts) - 2000
plot_range_end <- max(all_ends) + 2000

snoglobe_output <- read.table('data/snoGloBe/snoglobe_116.tsv', sep='\t', header=T)

snoglobe_ranges <- as_granges(snoglobe_output[c('target_chromo','target_window_start','target_window_end','score')],
                              seqnames = target_chromo, start = target_window_start, end = target_window_end)

snoglobe_gene_peaks <- join_overlap_intersect(GRanges(seqnames = gene_chromosome, ranges = IRanges(start = plot_range_start, end = plot_range_end)), snoglobe_ranges) 

snoglobe_annotation <- AnnotationTrack(
  range = snoglobe_gene_peaks,
  genome = genome,
  chromosome = gene_chromosome,
  name = 'snoGloBe Peaks',
  col = '#524ACC',
  fill = '#A19BDE'
)

# finalize track_list (ie. remove NAs)
#  Plot the gene structure
track_list <- list(geneTrack, SYS5_annotation, HEK_annotation, snoglobe_annotation)
track_list <- track_list[!is.na(track_list)]


# Highlight track to emphasize three consecutive high-scoring windows
snoglobe_output_consecutive <- snoglobe_output %>%
  group_by(target_id) %>%
  mutate(
    high_score = score > 0.98,
    next1_seq = lead(target_window_start, 1) == target_window_start + 1,
    next2_seq = lead(target_window_start, 2) == target_window_start + 2,
    next1_high = lead(high_score, 1, default = FALSE),
    next2_high = lead(high_score, 2, default = FALSE)
  ) %>%
  filter(high_score & next1_high & next2_high & next1_seq & next2_seq)


snoglobe_consecutive_windows_ranges <- as_granges(snoglobe_output_consecutive[c('target_chromo','target_window_start','target_window_end','score')],
                                                  seqnames = target_chromo, start = target_window_start, end = target_window_end)

snoglobe_consecutive_peaks <- join_overlap_intersect(GRanges(seqnames = gene_chromosome, ranges = IRanges(start = plot_range_start, end = plot_range_end)), 
                                                     snoglobe_consecutive_windows_ranges)
# make wider so its visible
snoglobe_consecutive_peaks <- flank(snoglobe_consecutive_peaks, 2, both=T)

consecutive_window_highlight <- HighlightTrack(trackList = track_list, 
                                               range = snoglobe_consecutive_peaks,
                                               chromosome = gene_chromosome,
                                               col = 'red'
                                               )


plotTracks(
  append(ideogram_track, consecutive_window_highlight),
  from = plot_range_start, 
  to = plot_range_end,
  transcriptAnnotation = "transcript"
  #exonAnnotation = 'exon'
)







#### (not used)
# but can auto-plot BioMart tracks?
afrom <- 2960000
ato <- 3160000

bmt <- BiomartGeneRegionTrack(genome = "hg38", chromosome = "chr12",
                              start = afrom, end = ato,
                              filter = list(with_ox_refseq_mrna = TRUE),
                              stacking = "dense")
plotTracks(bmt, from = afrom, to = ato, 
           chromosome = "chr12")



























