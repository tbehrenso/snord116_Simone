library(IsoformSwitchAnalyzeR)
library(dplyr)
library(readxl)
library(sva)

# Import Salmon data
# salmon_quant <- importIsoformExpression(parentDir = 'data/salmon_quants')

# Import Kallisto Data
kallisto_quant <- importIsoformExpression(parentDir = 'data/kallisto_quants_notDeep')

# Import subset of Kallisto Data
kallisto_quant <- importIsoformExpression(parentDir = 'data/kallisto_quants_Gilmore_H9/')

colnames(kallisto_quant[[1]])

# SVA Surrogate Variable Identification
if(F){
  sampleMetadata <- design_matrix[2:3]
  rownames(sampleMetadata) <- design_matrix$sampleID
  
  mod <- model.matrix(~ condition, data=sampleMetadata)
  mod0 <- model.matrix(~ 1, data=sampleMetadata)
  
  countMatrix <- as.matrix(kallisto_quant$counts[-1])
  rownames(countMatrix) <- kallisto_quant$counts$isoform_id
  countMatrix <- countMatrix[rowSums(countMatrix)>0,]
  
  svobj <- sva(countMatrix, mod, mod0)
  
  sv1 <- svobj$sv[,1]
  sampleMetadata$SV1 <- sv1
  
  boxplot(SV1 ~ condition, data = sampleMetadata,
          ylab = "SV1", main = "SV1 across Batches")
}

# Design matrix - Simone's Data
design_matrix <- data.frame(
  sampleID = colnames(kallisto_quant$abundance)[-1],
  condition = c(rep('CTRL',6), rep('PW',3))          # add column for batch effects??
)
# Design Matrix - Gilmore Data
design_matrix <- data.frame(
  sampleID = colnames(kallisto_quant$abundance)[-1],
  condition = c('smDEL','smDEL','smDEL','smDEL','smDEL','CTRL','CTRL','CTRL','CTRL','CTRL','CTRL')
  #condition = c('smDEL','smDEL','smDEL','smDEL','CTRL','smDEL','smDEL','CTRL','CTRL','CTRL','CTRL','CTRL')          # add column for batch effects??
  #condition = c('smDEL','smDEL','smDEL','smDEL','smDEL','CTRL','CTRL','CTRL','CTRL','CTRL','CTRL',
  #              'smDEL','smDEL','smDEL','smDEL','CTRL','smDEL','smDEL','CTRL','CTRL','CTRL','CTRL','CTRL'),
  #cellLine = c(rep('H9',11), rep('CT2',12)),
  #sv1 <- svobj$sv[,1]
  )

# Create switchAnalyzeRlist
aSwitchList <- importRdata(
  isoformCountMatrix   = kallisto_quant$counts,
  isoformRepExpression = kallisto_quant$abundance,
  designMatrix         = design_matrix,
  isoformExonAnnoation = '/Users/tbehr/Desktop/gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz',
  isoformNtFasta       = '/Users/tbehr/Desktop/gencode.v48.transcripts.fa.gz',
  showProgress = FALSE,
  ignoreAfterPeriod = TRUE,
  removeNonConvensionalChr = TRUE
)

# Filter
SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

# test for Isoform switches with DexSeq         
SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE,
  alpha = 0.5
)

# Extract Nucleotide and Amino Acid Sequences (write to file to analyze externally)         # files called "isoformSwitchAnalyzeR_isoform_AA.fasta" and
SwitchListAnalyzed <- extractSequence(                                                     #   "isoformSwitchAnalyzeR_isoform_nt.fasta
  SwitchListAnalyzed, 
  pathToOutput = 'data/',
  writeToFile=TRUE                                                           # RData object is from right AFTER this step (SwitchListAnalyzed.RData)
)


########## External Analyses Here ##############
# Then need to import and add back into the IsoformSwitchAnalyzeR object



## Analyze Alternative Splicing 
SwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchListAnalyzed
)
# Table of Intron Retention (IR)
table(SwitchListAnalyzed$AlternativeSplicingAnalysis$IR)

## Predict Switch Consequences
consequencesOfInterest <- c('tss','tts','exon_number')
SwitchListAnalyzed <- analyzeSwitchConsequences(
  SwitchListAnalyzed,
  consequencesToAnalyze = consequencesOfInterest, 
  dIFcutoff = 0.1,
  showProgress=TRUE
)

extractSwitchSummary(SwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = FALSE)

## Analysis of Individual Isoform Switching

SwitchListAnalyzedSubset <- subsetSwitchAnalyzeRlist(
  SwitchListAnalyzed, 
  SwitchListAnalyzed$isoformFeatures$condition_1 == 'CTRL'
)
SwitchListAnalyzedSubset

top_switches_extracted <- extractTopSwitches(
  SwitchListAnalyzedSubset, 
  filterForConsequences = FALSE, 
  n = NA, 
  sortByQvals = TRUE,
  alpha = 0.5
)

subset(top_switches_extracted, gene_name=='KLF13')

switchPlot(SwitchListAnalyzedSubset, gene = 'HCN1', additionalArguments = list(cex=1.5))

# View raw (?) transcript isoform counts
SwitchListAnalyzed$isoformCountMatrix[SwitchListAnalyzed$isoformCountMatrix$isoform_id=='ENST00000307145',]

gene_symbol_vector <- c("HEATR4", "PABIR2", "GRAMD1A", "PSMB10", "LIN54", "RPLP1", "RBM20", "MKNK2", "PRELID1", 
                        "NLE1", "SLC11A1", "PPA1", "SETDB2", "CCDC86", "EIF3C", "ZNF439", "RPSA", 
                        "GAA", "OAF", "UNC5B", "RPS4X", "HSPA8", "RPL17-C18orf32", "IFITM2", "RNASEH1-DT", "ZNF564",'GIT1','SRSF2','CTNNB1','GABRA2','ACTG1','PLEKHM2','NCL')



# Export as new column on the Excel
ase1_excel <- read_excel('data/ase1_to_peaks_complementarity_ALL.xlsx')

ase1_excel$ISA_GilmoreH9_qval <- top_switches_extracted$gene_switch_q_value[match(ase1_excel$geneSymbol, top_switches_extracted$gene_name)]

write.csv(ase1_excel, '/Users/tbehr/Desktop/ase1_excel_updated.csv')



switchPlotTopSwitches(
  switchAnalyzeRlist = SwitchListAnalyzedSubset, 
  n = 965,                                             # Set to Inf for all
  filterForConsequences = FALSE,
  fileType = "png",                                   # alternative is "png"
  pathToOutput = "/Users/tbehr/Desktop/"
)





















