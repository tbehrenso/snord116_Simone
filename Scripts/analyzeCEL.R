library(affy)
library(limma)
library(org.Mm.eg.db)
library(mogene10sttranscriptcluster.db)


setwd("C:/Users/tbehr/Desktop/SMN_CoIP")

Data <- ReadAffy()

eset <- rma(Data)

pData(eset)

strain <- c("FLSMN","FLSMN","FLSMN","CTRL","CTRL")

design <- model.matrix(~0 + factor(strain))

colnames(design) <- c("CTRL","FLSMN")

fit <- lmFit(eset, design)

fit <- eBayes(fit)



##
contrast.matrix <- makeContrasts(FLSMNvsCTRL = FLSMN - CTRL, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, number=Inf, adjust.method="BH")
##



options(digits=2)

res <- topTable(fit, number=Inf, adjust.method="BH", coef=1)# save your results

#write.table(res,"dif_exp.txt",sep="\t")

gene_symbols <- mapIds(mogene10sttranscriptcluster.db,
                       keys = rownames(res),
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

res$geneName <- gene_symbols

res[na.omit(gene_symbols=='Snord116l1'),]

gene_symbols[grep('Snord', gene_symbols)]
res[grep('Snord', gene_symbols),]










