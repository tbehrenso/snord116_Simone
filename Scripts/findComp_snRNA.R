library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
library(tidyverse)
library(biomaRt)
library(rentrez)

ASE1_SEQUENCE <- 'AACATTCCTTGGAAAAG'
ASE2_SEQUENCE <- 'CGTCATTCTCATCGGAA'
SNORD116_SEQUENCE <- 'TGGATCGATGATGAGTCCCCTATAAAAACATTCCTTGGAAAAGCTGAACAAAATGAGTGAGAACTCATAACGTCATTCTCATCGGAACTGAGGTCCA'
ase1 <- DNAString(ASE1_SEQUENCE)
ase2 <- DNAString(ASE2_SEQUENCE)
snord116 <- DNAString(SNORD116_SEQUENCE)
HG38 <- BSgenome.Hsapiens.UCSC.hg38

query_sequence <- ase2


snRNA_list <- c('RNU5A-1', 'RNVU1-14', 'RNVU1-19', 'RNU5B-1', 'RNU4-2', 'RNVU1-17', 'RNU105C', 'RNU5E-1', 'RNA5S11', 'RNU2-1',
                'RNU1-4', 'RNU5D-1', 'RNA5S7', 'RNA5S15', 'RNA5S17', 'RNA5S4', 'RNA5S5', 'RNA5S6', 'RNA5S16', 'RNA5S9')
rRNA_list <- c('RNA5S1', 'RNA5S2', 'RNA5S3', 'RNA5S4', 'RNA5S5', 'RNA5S6', 'RNA5S7', 'RNA5S8', 'RNA5S9', 'RNA5S10', 'RNA5S11', 
               'RNA5S12','RNA5S13', 'RNA5S14', 'RNA5S15', 'RNA5S16', 'RNA5S17')

GENE_SYMBOLS_LIST <- rRNA_list


##########################################################
#     Identify the strand of each gene symbol (copied and edited from findComplementarity.R)
##########################################################
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = 'https://useast.ensembl.org')

gene_info <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'transcript_start', 'transcript_end', 'strand'),
                   filters = 'hgnc_symbol',
                   values = GENE_SYMBOLS_LIST,
                   mart = mart,
                   uniqueRows = T) %>% 
  filter(chromosome_name %in% c(as.character(1:23)))

# Load genome annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes <- genes(txdb)


##########################################################
#     Get sequences 
##########################################################

snrna_regions_granges <- GRanges(seqnames=sapply(gene_info$chromosome_name, function(x) paste0('chr',x)),
                                      ranges=IRanges(start = gene_info$transcript_start,
                                                     end = gene_info$transcript_end),
                                 strand=gene_info$strand)

snrna_sequences <- reverseComplement(getSeq(HG38, snrna_regions_granges))

##########################################################
#     Search for complementarity
##########################################################

alignment_list <- vector('list', length = dim(gene_info)[1])
score_vector <- rep(NA, dim(gene_info)[1])

for (i in seq_along(snrna_sequences)) {
  if (countPattern('N', snrna_sequences[[i]])){
    print(paste0('Warning: DNA string #', i, ' (', gene_info$hgnc_symbol[i], ', ',
                 gene_info$transcript_start[i], '-',gene_info$transcript_end[i],
                 ') contains N values. Remove these values (regardless of position, so careful!).'))
    snrna_sequences[[i]] <- gsub('N', '', snrna_sequences[[i]])
  }
  
  alignment <- pwalign::pairwiseAlignment(query_sequence, snrna_sequences[[i]], 
                                          type = 'global-local', 
                                          substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
                                            match = 2, mismatch = -1, baseOnly = TRUE),
                                          gapOpening = 0, 
                                          gapExtension = -1)
  
  alignment_list[[i]] <- alignment
  score_vector[i] <- score(alignment)
  if (i%%100==0){print(i)}
  #print(alignment)
}

to_export <- data.frame(
  geneSymbol = gene_info$hgnc_symbol,
  geneStrand = gene_info$strand,
  score = score_vector,
  substringLength = unlist(lapply(alignment_list, nchar), use.names = F),
  numberMatches = unlist(lapply(alignment_list, nmatch), use.names = F)
) %>% 
  mutate(percentMatch = (numberMatches / substringLength) * 100, .after = numberMatches)

write.csv(to_export, 'peaks_methylDatabase_crossed_ASE2_complementarity_local_noGap20.csv', row.names=F)




