library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
library(tidyverse)
library(biomaRt)
library(rentrez)

ASE1_SEQUENCE <- 'AACATTCCTTGGAAAAG'
ASE2_SEQUENCE <- 'CGTCATTCTCATCGGAA'
ase1 <- DNAString(ASE1_SEQUENCE)
ase2 <- DNAString(ASE2_SEQUENCE)
HG38 <- BSgenome.Hsapiens.UCSC.hg38

query_sequence <- ase2

peaks_methylDatabase_crossed <- read_excel('data/rmbase_methylBase_crossed.xlsx')
names(peaks_methylDatabase_crossed) <- c('Chromosome', 'Start', 'End', 'snoRNA', 'Gene')

##########################################################
#     Identify the strand of each gene symbol (copied and edited from findComplementarity.R)
##########################################################
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = 'https://useast.ensembl.org')

gene_info <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'strand'),
                   filters = 'hgnc_symbol',
                   values = peaks_methylDatabase_crossed$Gene,
                   mart = mart,
                   uniqueRows = T)

# Load genome annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes <- genes(txdb)

# Loop to order strand info, and fill in missing info
gene_strands <- rep(NA, length(peaks_methylDatabase_crossed$Gene))
for(i in seq(1, length(peaks_methylDatabase_crossed$Gene))){
  if (peaks_methylDatabase_crossed$Gene[i] %in% gene_info$hgnc_symbol){
    gene_strands[i] <- gene_info$strand[gene_info$hgnc_symbol == peaks_methylDatabase_crossed$Gene[i]]
  } else {
    # Convert to GRanges
    query_region <- GRanges(seqnames = peaks_methylDatabase_crossed$Chromosome[i], 
                            ranges = IRanges(start = peaks_methylDatabase_crossed$Start[i], 
                                             end = peaks_methylDatabase_crossed$End[i]))
    
    # Find overlaps with known genes
    overlaps <- findOverlaps(query_region, all_genes)
    
    if (length(overlaps) > 0) {
      matched_gene <- all_genes[subjectHits(overlaps)[1]]
      strand_info <- as.character(strand(matched_gene))
      gene_strands[i] <- case_when(strand_info=='+' ~ '1', strand_info=='-' ~ '-1') # convert to +1 and -1 for consistency
    } else {
      print(paste0(i, ': Gene in neither ENSEMBL nor UCSC TxDb annotations'))
      gene_strands[i] <- NA
    }
  }
}

##########################################################
#     Get sequences 
##########################################################

methylated_regions <- as.data.frame(peaks_methylDatabase_crossed[c('Chromosome','Start','End')]) %>% 
  mutate(strand = gene_strands) %>% 
  mutate(trueStart = ifelse(gene_strands=='-1' | is.na(gene_strands), Start - 12, Start - 4)) %>% 
  mutate(trueEnd = ifelse(gene_strands=='-1' | is.na(gene_strands), Start + 4, Start + 12)) 

methylated_regions_granges <- GRanges(seqnames=methylated_regions$Chromosome, 
                                ranges=IRanges(start = methylated_regions$trueStart, 
                                               end = methylated_regions$trueEnd))

methylated_sequences <- getSeq(HG38, methylated_regions_granges)

# switch to other strand if that gene is on reverse strand
for(i in seq(1:length(gene_strands))){
  if(gene_strands[i]=='1' | is.na(gene_strands[i])){
    methylated_sequences[i] <- reverseComplement(methylated_sequences[i])
  }
}

methylated_sequences[which(peaks_methylDatabase_crossed$Gene=='RNU5A-1')[1]]


##########################################################
#     Search for complementarity
##########################################################

alignment_list <- vector('list', length = dim(peaks_methylDatabase_crossed)[1])
score_vector <- rep(NA, dim(peaks_methylDatabase_crossed)[1])

for (i in seq_along(methylated_sequences)) {
  if (countPattern('N', methylated_sequences[[i]])){
    print(paste0('Warning: DNA string #', i, ' (', peaks_methylDatabase_crossed[i,1], ', ',
                 peaks_methylDatabase_crossed[i,2], '-',peaks_methylDatabase_crossed[i,3],
                 ') contains N values. Remove these values (regardless of position, so careful!).'))
    methylated_sequences[[i]] <- gsub('N', '', methylated_sequences[[i]])
  }
  
  alignment <- pwalign::pairwiseAlignment(query_sequence, methylated_sequences[[i]], 
                                          type = 'global', 
                                          substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
                                            match = 2, mismatch = -2, baseOnly = TRUE),
                                          gapOpening = -5, 
                                          gapExtension = -3)

  alignment_list[[i]] <- alignment
  score_vector[i] <- score(alignment)
  if (i%%100==0){print(i)}
  #print(alignment)
}

to_export <- data.frame(
  geneSymbol = peaks_methylDatabase_crossed$Gene,
  geneStrand = gene_strands,
  methylSite = peaks_methylDatabase_crossed$Start,
  score = score_vector,
  substringLength = unlist(lapply(alignment_list, nchar), use.names = F),
  numberMatches = unlist(lapply(alignment_list, nmatch), use.names = F)
) %>% 
  mutate(percentMatch = (numberMatches / substringLength) * 100, .after = numberMatches)

write.csv(to_export, 'peaks_methylDatabase_crossed_ASE1_complementarity_global_noGap20.csv', row.names=F)


