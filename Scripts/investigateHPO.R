library(DESeq2)
library(AnnotationDbi)
library(AnnotationHub)
library(DBI)


# Load dexseq results and extract gene list

res <- readRDS('/Users/tbehr/Desktop/SanRaffaele/Projects/snord116_Simone/OldMisc/dexseq_results_saved.RData')
res <- res[!is.na(res$log2fold_PW1_CTRL) & !is.na(res$padj),]
#res_sig <- res[res$padj < 0.05 & abs(res$log2fold_PW1_CTRL) > 1,]
res_sig <- res[res$padj < 0.05,]

ensembl_unique <- unique(res_sig$groupID)

ensembl_unique_nodot <- sapply(strsplit(ensembl_unique, "\\."), `[`, 1)

# Convert ensembl to NCBI IDs

library(org.Hs.eg.db)

genelist_ncbi <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_unique_nodot,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)


# Prep HPO database

ah <- AnnotationHub()
query(ah, "human phenotype ontology")
hpo <- ah[["AH117056"]]

ah[['AH117056']] # this part doesn't work right, so need to manually load the database file

library(RSQLite)
filename <- "C:/Users/tbehr/AppData/Local/R/cache/R/AnnotationHub/1b44c446407_123802"
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,
                dbname = filename)
dbListTables(db)
mytable <- dbReadTable(db,"hpoid_gene")
head(dbReadTable(db,"do_term"))

dbGetQuery(db, 'SELECT * FROM hpoid_gene WHERE "gene" = 10')


hpo2gene <- dbReadTable(db,"hpoid_gene")
hpo2name <- dbReadTable(db,"do_term")


library(clusterProfiler)
library(ggplot2)
library(stringr)

hpo_enrichment <- enricher(genelist_ncbi,
         TERM2GENE = hpo2gene,
         TERM2NAME = hpo2name,
         pvalueCutoff = 0.05)
hpo_res <- hpo_enrichment@result


barplot(hpo_enrichment)
hpo_dotplot <- dotplot(hpo_enrichment, showCategory = 30)

# automatically select terms from file (from HPO database) associated with Prader-Willi
phenotypes_for_pw <- read.csv('data/phenotypes_for_OMIM_176270', sep='\t')
highlight_terms <- hpo_res$Description[hpo_res$p.adjust<0.05][hpo_res$Description[hpo_res$p.adjust<0.05] %in% trimws(phenotypes_for_pw$name)]

highlight_terms <- c('Delayed gross motor development', 'Aggressive behaviour', 'High forehead', 'Intellectual disability, severe',
                     'Short attention span', 'Thin vermilion border', 'Long philtrum', 'Attention deficit hyperactivity disorder',
                     'Hypoplasia of the corpus callosum', ' Abnormal nasal tip morphology', 'Downslanted palpebral fissures', 'Thin upper lip vermilion',
                     'Aggressive behavior', 'Narrow forehead')


# Modify label colors
hpo_dotplot + theme(
  axis.text.y = element_text(
    colour = ifelse(levels(hpo_dotplot$data$Description) %in% highlight_terms, "red", "black")
  )
) + scale_y_discrete(labels=function(x) str_wrap(x, width=40))


# Extract genes related to the important HPO terms

intDisSev_genes <- strsplit(hpo_res$geneID[hpo_res$ID == 'HP:0010864'], split = '/')[[1]]
delGroMotDev_genes <- strsplit(hpo_res$geneID[hpo_res$ID == 'HP:0002194'], split = '/')[[1]]
ShoAttSpa_genes <- strsplit(hpo_res$geneID[hpo_res$ID == 'HP:0000736'], split = '/')[[1]]

three_term_genes <- intersect(intersect(intDisSev_genes, delGroMotDev_genes), ShoAttSpa_genes)
two_term_genes_a <- intersect(intDisSev_genes, delGroMotDev_genes)
two_term_genes_a <- setdiff(two_term_genes_a, three_term_genes)
two_term_genes_b <- intersect(intDisSev_genes, ShoAttSpa_genes)
two_term_genes_b <- setdiff(two_term_genes_b, c(three_term_genes,two_term_genes_a))
two_term_genes_c <- intersect(delGroMotDev_genes, ShoAttSpa_genes)
two_term_genes_c <- setdiff(two_term_genes_c, c(three_term_genes,two_term_genes_a,two_term_genes_b))

three_term_genes_ensembl <- names(genelist_ncbi[genelist_ncbi %in% three_term_genes])
two_term_genes_a_ensembl <- names(genelist_ncbi[genelist_ncbi %in% two_term_genes_a])
two_term_genes_b_ensembl <- names(genelist_ncbi[genelist_ncbi %in% two_term_genes_b])
two_term_genes_c_ensembl <- names(genelist_ncbi[genelist_ncbi %in% two_term_genes_c])
two_term_genes_ensembl <- unique(c(two_term_genes_a_ensembl,two_term_genes_b_ensembl,two_term_genes_c_ensembl))

c('Intellectual disability, severe', 'Delayed gross motor development', 'Short attention span')
c('HP:0010864','HP:0002194','HP:0000736')



