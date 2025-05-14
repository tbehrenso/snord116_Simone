library('TSRexploreR')
library('ggplot2')
library(GenomicFeatures)

#--------------------------------------------
#        Prep files and TSR object
#--------------------------------------------

# load annotation (gtf) and assembly (fasta)
annotation <- read.table('/beegfs/scratch/ric.broccoli/behrens.thomas/reference/gencode.v47.basic.annotation.nochr.gtf', sep='\t')
assembly <- '/beegfs/scratch/ric.broccoli/behrens.thomas/reference/GRCh38.primary_assembly.genome.nochr.fa'

# load samples
samples <- data.frame(sample_name = c('EDO_1','EDO_2','EDO_3','ND1', 'ND2', 'ND3','PW1_1','PW1_2','PW1_3'),
                      file_1 = c('/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/Bam_file/EDO_1.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/Bam_file/EDO_2.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/Bam_file/EDO_3.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/Bam_file/ND1_1.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/Bam_file/ND1_2.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/Bam_file/ND1_3.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/Bam_file/PW1_1.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/Bam_file/PW1_2.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/Bam_file/PW1_3.bam'
                      ),
                      file_2 = NA,
                      conditions = c(rep('CTRL',6), rep('PW1', 3))
)


exp_base <- tsr_explorer(sample_sheet = samples, genome_assembly = assembly, genome_annotation = annotation) %>% 
  format_counts(data_type='tss')


#--------------------------------------------
#        Process BAMs
#--------------------------------------------

exp <- import_bams(exp_base, paired=T, proper_pair=T, remove_duplicate=T, remove_secondary=T)

# plots for soft-clipped base assessment
#gg_softclip_histogram <- softclip_histogram(exp)+
#  theme_bw() +
#  scale_fill_viridis_d()

#ggsave('softclip_histogram.png', gg_softclip_histogram, device="png", path='Plots', width=12, height=8)


#gg_softclip_composition <- softclip_composition(exp) +
#  theme_bw() +
#  scale_fill_viridis_d()

#ggsave('softclip_composition.png', gg_softclip_composition, device="png", path='Plots', width=12, height=8)

#print('First Plots Plotted')

# G correction (may be unnecessary, check plots above after first run)
# exp <- G_correction(exp)


#--------------------------------------------
#        TSS Shifting
#--------------------------------------------

# Identify/aggregate TSSs
exp <- tss_aggregate(exp)

saveRDS(exp, 'exp_tssAggregated.RData')



# cluster TSSs and merge replicates
exp <- tss_clustering(exp, threshold=3, n_samples=1) %>% 
  merge_samples(data_type='tss', merge_group='conditions') %>% 
  merge_samples(data_type='tsr', merge_group='conditions')

# calculate shifting
exp <- tss_shift(
  exp,
  sample_1=c(TSS='CTRL', TSR='CTRL'),
  sample_2=c(TSS='PW1', TSR='PW1'),
  comparison_name = 'CTRL_vs_PW1',
  max_distance = 100, min_threshold = 10, n_resamples = 1000L
)

saveRDS(exp, 'exp_tssShifted.RData')

# plot shift scores ranked
#gg_shift_rank <- plot_shift_rank(exp) +
#  scale_fill_viridis_c() +
#  theme_bw() +
#  geom_hline(yintercept = 0)

#ggsave('shift_rank.png', gg_shift_rank, device="png", path='Plots', width=12, height=8)

# plot upstream and downstream shifts
#gg_shift_count <- plot_shift_count(exp) +
#  scale_fill_viridis_d() +
#  theme_bw()

#ggsave('shift_count.png', gg_shift_count, device="png", path='Plots', width=12, height=8)


# Annotate Shifts
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
library('ChIPseeker')

exp <- annotate_features(exp, data_type='shift', feature_type='transcript')

gg_genomic_distribution <- plot_genomic_distribution(exp, data_type='shift') +
  scale_fill_viridis_d(direction=-1)

ggsave('genomic_distribution.png', gg_genomic_distribution, device="png", path='Plots', width=12, height=8)



#--------------------------------------------
#        Differential TSSs
#--------------------------------------------

exp2 <- readRDS('exp_tssAggregated.RData')

exp2 <- exp2 %>% 
  normalize_counts(data_type='tss', method='DESeq2') %>% 
  tss_clustering(threshold=3, max_distance=25)

saveRDS(exp, 'exp2_tssClustered.RData')

print('exp2 Generated and Clustered')

gg_corr_plot <- plot_correlation(
  exp2, data_type='tsr', font_size=12,
  use_normalized=TRUE,
  heatmap_colors=viridis::viridis(100)
)

ggsave('corr_plot.png', gg_corr_plot, device='png', path='Plots', width=12, height=8)


gg_PCA <- plot_reduction(exp2, data_type='tsr', remove_var=0.25, colby='condition')

ggsave('pca_plot.png', gg_PCA, device='png', path='Plots', width=12, height=8)

print('Plots Plotted 2')

exp2 <- fit_de_model(exp2, data_type='tsr', formula=~condition, method='DESeq2')

exp2 <- differential_expression(
  exp2, data_type = 'tsr',
  comparison_name = 'Ctrl vs PW1',
  comparison_type = 'contrast',
  comparison = c('conditions', 'CTRL', 'PW1')
)

exp2 <- annotate_features(exp2, data_type='tsr_diff', feature_type='transcript', upstream=250, downstream=250)

gg_ma_plot <- plot_ma(exp2, data_type='tsr') +
  scale_color_viridis_d() +
  theme_bw()

ggsave('ma_plot.png', gg_ma_plot, device='png', path='Plots', width=12, height=8)

gg_volcano_plot <- plot_volcano(exp2, data_type='tsr') +
  scale_color_viridis_d() +
  theme_bw()

ggsave('volcano_plot.png', gg_volcano_plot, device='png', path='Plots', width=12, height=8)

de_features_plot <- plot_num_de(exp2, data_type='tsr', de_comparisons='Ctrl vs PW1') +
  scale_fill_viridis_d() +
  theme_bw()

ggsave('feature_plot.png', de_features_plot, device='png', path='Plots', width=12, height=8)

print('DE Plots Saved')

# Gene Ontology

enrichment_data <- export_for_enrichment(exp, data_type='tsr', anno_categories='Promoter', keep_unchanged=TRUE)

BiocManager::install("clusterProfiler")
BiocManager::install("org.Sc.sgd.db")
library('clusterProfiler')
library('org.Sc.sgd.db')

go_enrichment <- compareCluster(
  geneId ~ sample + de_status,
  data=enrichment_data,
  fun='enrichGO',
  OrgDb='org.Sc.sgd.db',
  pAdjustMethod='fdr',
  ont='BP',
  keyType='ENSEMBL'
)

print('GOs Enriched')

gg_go_enrichment <- dotplot(go_enrichment) +
  scale_color_viridis_c() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave('go_enrichment.png', gg_go_enrichment, device='png', path='Plots', width=12, height=8)


















