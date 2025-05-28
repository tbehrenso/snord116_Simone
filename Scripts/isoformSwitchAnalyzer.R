library(IsoformSwitchAnalyzeR)
library(dplyr)


# Import Salmon data
salmon_quant <- importIsoformExpression(parentDir = 'data/salmon_quants')

# Design matrix
design_matrix <- data.frame(
  sampleID = colnames(salmon_quant$abundance)[-1],
  condition = c(rep('CTRL',6), rep('PW',3))          # add column for batch effects??
)

# Create switchAnalyzeRlist
aSwitchList <- importRdata(
  isoformCountMatrix   = salmon_quant$counts,
  isoformRepExpression = salmon_quant$abundance,
  designMatrix         = design_matrix,
  'data/salmon_ensembl.gtf.gz',  # isoformExonAnnoation = '/Users/tbehr/Desktop/gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz',  # 
  #isoformNtFasta       = '/Users/tbehr/Desktop/gencode.v48.transcripts.fa.gz',  # 'data/salmon_ntfasta.fa.gz',  # 
  showProgress = FALSE,
  ignoreAfterPeriod = TRUE,
  removeNonConvensionalChr = TRUE
)




