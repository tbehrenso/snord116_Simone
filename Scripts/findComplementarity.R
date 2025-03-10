library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
library(tidyverse)
library(biomaRt)
library(rentrez)
library(plyranges)
source('Scripts/complimentarityFunctions.R')

SNORD116_SEQUENCE <- 'TGGATCGATGATGAGTCCCCTATAAAAACATTCCTTGGAAAAGCTGAACAAAATGAGTGAGAACTCATAACGTCATTCTCATCGGAACTGAGGTCCA'
ASE1_SEQUENCE <- 'AACATTCCTTGGAAAAG'
ASE2_SEQUENCE <- 'CGTCATTCTCATCGGAA'
HG38 <- BSgenome.Hsapiens.UCSC.hg38

# get reverse complement of relevant sequences
snord116 <- DNAString(SNORD116_SEQUENCE)
snord116_revcomp <- reverseComplement(snord116)
ase1 <- DNAString(ASE1_SEQUENCE)
ase1_revcomp <- reverseComplement(ase1)
ase2 <- DNAString(ASE2_SEQUENCE)
ase2_revcomp <- reverseComplement(ase2)

### SELECT which sequences to use
query_sequence <- ase2

SYS5_differential <- read_excel('data/peaksAnnoFULL.xlsx', sheet='SYS5_differential')

##########################################################
#     Identify the strand of each gene symbol 
##########################################################
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)

mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

gene_info <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'strand'),
                   filters = 'hgnc_symbol',
                   values = SYS5_differential$SYMBOL,
                   mart = mart,
                   uniqueRows = T)

# Load genome annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes <- genes(txdb)

# Loop to order strand info, and fill in missing info
gene_strands <- rep(NA, length(SYS5_differential$SYMBOL))
for(i in seq(1, length(SYS5_differential$SYMBOL))){
  if (SYS5_differential$SYMBOL[i] %in% gene_info$hgnc_symbol){
    gene_strands[i] <- gene_info$strand[gene_info$hgnc_symbol == SYS5_differential$SYMBOL[i]]
  } else {
    # Convert to GRanges
    query_region <- GRanges(seqnames = SYS5_differential$seqnames[i], 
                            ranges = IRanges(start = SYS5_differential$start[i], end = SYS5_differential$end[i]))
    
    # Find overlaps with known genes
    overlaps <- findOverlaps(query_region, all_genes)
    
    if (length(overlaps) > 0) {
      matched_gene <- all_genes[subjectHits(overlaps)[1]]
      strand_info <- as.character(strand(matched_gene))
      gene_strands[i] <- case_when(strand_info=='1' ~ '+', strand_info=='-1' ~ '-') # convert for consistency (doesn't work)
    } else {
      print(paste0(i, ': Gene in neither ENSEMBL nor UCSC TxDb annotations'))
      gene_strands[i] <- NA
    }
  }
}

gene_strands_corrected <- sapply(gene_strands,
                                 function(x) case_when(x=='1'~'+', x=='-1'~'-', is.na(x)~'+'), USE.NAMES = F)

##########################################################
#     Extract sequences for where peaks are located
##########################################################

SYS5_differential <- read_excel('peaksAnno.xlsx', sheet='SYS5_differential')

peak_regions <- as.data.frame(SYS5_differential[c('seqnames','start','end')])

peak_regions_granges <- GRanges(seqnames=peak_regions$seqnames,
                                ranges=IRanges(start = peak_regions$start, end = peak_regions$end),
                                strand = gene_strands_corrected)

peak_sequences <- reverseComplement(getSeq(HG38, peak_regions_granges))

##########################################################
#     Find matches between given sequence and the peaks' sequences
##########################################################

#matches <- vmatchPattern(snord116_revcomp, peak_sequences, with.indels=TRUE)
alignment_list <- vector('list', length = dim(peak_regions)[1])
score_vector <- rep(NA, dim(peak_regions)[1])

for (i in seq_along(peak_sequences)) {
  if (countPattern('N', peak_sequences[[i]])){
    print(paste0('Warning: DNA string #', i, ' (', peak_regions[i,1], ', ',peak_regions[i,2], '-',peak_regions[i,3],
                 ') contains N values. Remove these values (regardless of position, so careful!).'))
    peak_sequences[[i]] <- gsub('N', '', peak_sequences[[i]])
  }
    
  alignment <- pwalign::pairwiseAlignment(query_sequence, peak_sequences[[i]], 
                                 type = 'global-local', 
                                 substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
                                   match = 2, mismatch = -2, baseOnly = TRUE),
                                 gapOpening = -5, 
                                 gapExtension = -3)

  alignment_list[[i]] <- alignment
  score_vector[i] <- score(alignment)
  if (i%%500==0){print(i)}
}

to_export <- data.frame(
  geneSymbol = SYS5_differential$SYMBOL,
  geneStrand = gene_strands,
  score = score_vector,
  substringLength = unlist(lapply(alignment_list, nchar), use.names = F),
  numberMatches = unlist(lapply(alignment_list, nmatch), use.names = F),
  annotation = SYS5_differential$annotation,
  distanceToTSS = SYS5_differential$distanceToTSS,
  foldEnrichment = SYS5_differential$fold_enrichment
) %>% 
  mutate(percentMatch = (numberMatches / substringLength) * 100, .after = numberMatches)


##########################################################
#     Cross-reference with other datasets
##########################################################

### For exons, identify if peak overlaps with diff. exon expression (taking +-200bp from the edges of the exon)
dex_seq <- read_excel('data/DexSeq_with_coordinates.xlsx')

dexseq_ranges <- as_granges(dex_seq[c('genomicData.seqnames','genomicData.start','genomicData.end','log2fold_PW1_CTRL')],
                            seqnames = genomicData.seqnames, start = genomicData.start-200, end = genomicData.end+200)
peaks_ranges <- as_granges(SYS5_differential[c('seqnames','start','end')],
                           seqnames = seqnames, start = start, end = end)

# get overlap, and filter duplicate rows (keep larger absolute value)
dex_peaks_overlap <- join_overlap_left(peaks_ranges, dexseq_ranges) %>% 
  as.data.frame %>% 
  group_by(seqnames, start, end) %>% 
  slice_max(order_by = log2fold_PW1_CTRL, n=1) %>% 
  ungroup()

to_export$log2fold <- dex_peaks_overlap$log2fold_PW1_CTRL

### Check if peaks overlap with a methylation site
rm_cytosolic <- read_excel('data/RMBase/RMBase_cytosolic.xlsx', skip=1)
rm_other <- read_excel('data/RMBase/RMBase_Other.xlsx', skip=1)
rm_hek <- read_excel('data/RMBase/RMBase_HEK.xlsx', skip=1)
rm_hela <- read_excel('data/RMBase/RMBase_HeLa.xlsx', skip=1)
rm_pa1 <- read_excel('data/RMBase/RMBase_PA1.xlsx', skip=1)

rmbase_all <- bind_rows(lst(rm_cytosolic, rm_other, rm_hek, rm_hela, rm_pa1), .id='rmBase')

rmbase_ranges <- as_granges(rmbase_all[c('ModChr','ModStart','ModEnd', 'Strand', 'rmBase')], 
                            seqnames = ModChr, start = ModStart, end = ModEnd, strand = Strand)

rmbase_peaks_overlap <- join_overlap_left(peaks_ranges, rmbase_ranges) %>% 
  as.data.frame() %>% 
  group_by(seqnames, start, end) %>% 
  summarize(rmBase = paste(rmBase, collapse = ', '), .groups = 'drop')

to_export$rmbaseSite <- rmbase_peaks_overlap$rmBase

### For 3' UTR, determine if present in RBA abundance (do I have this file?)

### For promoter....???????


##########################################################
#     Get a weighted scoring based on position (based on Chen 2007)
##########################################################

weighted_score_vector <- rep(NA, length(alignment_list))

for(i in seq_along(alignment_list)){
  query_aligned_sequence <- toString(alignment_list[[i]])           # Need to use toString (or as.character) directly on alignment_list[[i]] to get the full sequence. However, its WAY slower than pattern()
  subject_aligned_sequence <- as.character(subject(alignment_list[[i]]))
  
  weighted_score_vector[i] <- get_weighted_score(query_aligned_sequence, subject_aligned_sequence)
  
  if (i%%1000==0){print(paste0('Weighted scoring: alignment #', i))}
}

to_export$weightedScore <- weighted_score_vector

filename = 'data/ase2_to_peaks_complementarity_ALL.csv'
if(file.exists(filename)){
  print('Warning: File already exists! Will NOT overwrite')
} else {
  write.csv(to_export, filename, row.names=F)
}

##########################################################
#     Assess complementarity of the ASE sequences with the identified regions with high complementarity with snord116
##########################################################
# The idea here is to see if the identified regions have complementarity with one or both ASE sequences

if(F){
  alignments_highscore <- alignment_list[score_vector>100]
  gene_info_highscore <- to_export[score_vector>100,]
  gene_strands_highscore <- gene_strands[score_vector>100]
  score_vector_highscore <- score_vector[score_vector>100]
  
  alignment_subjects <- sapply(alignments_highscore, function(x) str_replace_all(as.character(subject(x)), '-',''))
  
  alignment_list_ase1 <- vector('list', length = length(alignment_subjects))
  alignment_list_ase2 <- vector('list', length = length(alignment_subjects))
  score_vector_ase1 <- rep(NA, length(alignment_subjects))
  score_vector_ase2 <- rep(NA, length(alignment_subjects))
  distance_between_ases <- numeric(length=length(alignment_subjects))
  
  for(i in seq(1:length(alignment_subjects))){
    
    query_aligned_sequence <- as.character(pattern(alignments_highscore[[i]]))
    subject_aligned_sequence <- as.character(subject(alignments_highscore[[i]]))
    ase1_start_index <- find_index_after_n_chars(query_aligned_sequence, 27)
    ase1_end_index <- find_index_after_n_chars(query_aligned_sequence, 43)
    ase2_start_index <- find_index_after_n_chars(query_aligned_sequence, 71)
    ase2_end_index <- find_index_after_n_chars(query_aligned_sequence, 87)
    distance_between_ases[i] <- ase2_start_index - ase1_end_index
    
    # find complementarity of the ASE sequences within the snord116 alignment
    subject_aligned_with_ase1 <- substring(subject_aligned_sequence, ase1_start_index, ase1_end_index)
    subject_aligned_with_ase2 <- substring(subject_aligned_sequence, ase2_start_index, ase2_end_index)
    
    ase1_subalignment <- pwalign::pairwiseAlignment(ase1, DNAString(str_replace_all(subject_aligned_with_ase1,'-','')), 
                                            type = 'global', 
                                            substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
                                              match = 2, mismatch = -1, baseOnly = TRUE),
                                            gapOpening = 0, 
                                            gapExtension = -1)
    ase2_subalignment <- pwalign::pairwiseAlignment(ase2, DNAString(str_replace_all(subject_aligned_with_ase2,'-','')), 
                                                    type = 'global', 
                                                    substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
                                                      match = 2, mismatch = -1, baseOnly = TRUE),
                                                    gapOpening = 0, 
                                                    gapExtension = -1)
    alignment_list_ase1[[i]] <- ase1_subalignment
    alignment_list_ase2[[i]] <- ase2_subalignment
    score_vector_ase1[i] <- score(ase1_subalignment)
    score_vector_ase2[i] <- score(ase2_subalignment)
    
  }
  
  to_export_fuller <- gene_info_highscore %>% 
    mutate(ase1_score = score_vector_ase1) %>% 
    mutate(ase1_percentMatch = unlist(lapply(alignment_list_ase1, pwalign::pid), use.names = F)) %>%
    mutate(ase1_indels = unlist(lapply(alignment_list_ase1, function(x) pwalign::nindel(pattern(x))[2]), use.names = F) +
             unlist(lapply(alignment_list_ase1, function(x) pwalign::nindel(subject(x))[2]), use.names = F)) %>% 
    mutate(ase2_score = score_vector_ase2) %>% 
    mutate(ase2_percentMatch = unlist(lapply(alignment_list_ase2, pwalign::pid), use.names = F)) %>%
    mutate(ase2_indels = unlist(lapply(alignment_list_ase2, function(x) pwalign::nindel(pattern(x))[2]), use.names = F) +
             unlist(lapply(alignment_list_ase2, function(x) pwalign::nindel(subject(x))[2]), use.names = F)) %>%
    mutate(distanceBetweenAses = distance_between_ases)
    
  
  
  filename = 'data/snord116_and_ases_to_peaks_HIGHSCORE.csv'
  if(file.exists(filename)){
    print('Warning: File already exists! Will NOT overwrite')
  } else {
    write.csv(to_export_fuller, filename, row.names=F)
  }
}





