library('TSRexploreR')
library('ggplot2')


#--------------------------------------------
#        Prep files and TSR object
#--------------------------------------------

# load annotation (gtf) and assembly (fasta)
annotation <- '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/gencode.v31.basic.annotation.gff'
assembly <- '/beegfs/scratch/ric.broccoli/behrens.thomas/reference/GRCh38.primary_assembly.genome.fa'

# load samples
samples <- data.frame(sample_name = c('EDO_1','EDO_2','EDO_3','PW1_1','PW1_2','PW1_3'),
                      file_1 = c('/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/EDO_1.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/EDO_2.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/EDO_3.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/PW1_1.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/PW1_2.bam',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/Bam_file/PW1_3.bam'
                      ),
                      file_2 = NA,
                      condition = c(rep('CTRL',3), rep('PW1', 3))
)

exp_base <- tsr_explorer(sample_sheet = samples, genome_assembly = assembly, genome_annotation = annotation)


#--------------------------------------------
#        Process BAMs
#--------------------------------------------

exp <- import_bams(exp_base, paired=T, proper_pair=T, remove_duplicate=T, remove_secondary=T)

print('Bams imported')

# # plots for soft-clipped base assessment
# gg_softclip_histogram <- softclip_histogram(exp)+
#  theme_bw() +
#  scale_fill_viridis_d()
# 
# ggsave('softclip_histogram.png', gg_softclip_histogram, device="png", path='Plots', width=12, height=8)


# gg_softclip_composition <- softclip_composition(exp) +
#   theme_bw() +
#   scale_fill_viridis_d()
# 
# ggsave('softclip_composition.png', gg_softclip_composition, device="png", path='Plots', width=12, height=8)

# print('First Plots Plotted')

# G correction (may be unnecessary, check plots above after first run)
# exp <- G_correction(exp)

# Identify/aggregate TSSs
exp <- tss_aggregate(exp)

saveRDS(exp, 'exp_tssAggregated_return.RData')

print('TSSs aggregated')


#--------------------------------------------
#        Differential Features
#--------------------------------------------

exp <- exp %>%
  format_counts(data_type = 'tss') %>% 
  normalize_counts(data_type = 'tss', method='DESeq2') %>% 
  tss_clustering(threshold=3, max_distance=25) %>% 


saveRDS(exp, 'exp_tssClustered_return.RData')

print('Counts Formatted, Normalized, Clustered')

# exp <- exp %>% 
#   merge_samples(data_type='tss', merge_group='condition') %>% 
#   merge_samples(data_type='tsr', merge_group='condition')
# 
# saveRDS(exp, 'exp_tssMerged_return.RData')
# 
# print('Samples Merged')

# calculate shifting
exp <- tss_shift(
  exp,
  sample_1=c(TSS='PW1', TSR='PW1'),
  sample_2=c(TSS='CTRL', TSR='CTRL'),
  comparison_name = 'CTRL_vs_PW1',
  max_distance = 100, min_threshold = 10, n_resamples = 1000L
)


saveRDS(exp, 'exp_tssShifted_return.RData')

print('TSS Shifts Calculated')


#--------------------------------------------
#        Analysis of TSR Similarity
#--------------------------------------------

# Correlation Matrix
gg_corrplot <- plot_correction(exp, data_type='tsr', font_size=12,
                use_normalized=TRUE,
                heatmap_colors=viridis::viridis(100)
                )

ggsave('correlation_heatmap.png', gg_corrplot, device="png", path='Plots', width=12, height=8)

# PCA
gg_pca <- plot_reduction(exp, data_type='tsr', remove_var=0.25, colby='condition')

ggsave('diffFeat_PCA.png', gg_pca, device="png", path='Plots', width=12, height=8)


#--------------------------------------------
#        Differential TSR Analysis
#--------------------------------------------


exp <- fit_de_model(exp, data_type='tsr', formula=~condition, method='DESeq2')


exp <- differential_expression(exp, 
                               data_type='tsr',
                               comparison_name='EDO_vs_PW1',
                               comparison_type='contrast',
                               comparison=c('condition', 'CTRL', 'PW1'))

saveRDS(exp, 'exp_diffExp_return.RData')

print('Model fitted')

# Annotate Differential Features
library('ChIPseeker')

exp <- annotate_features(exp, data_type='tsr_diff', feature_type='transcript', upstream=250, downstream=250)

# MA Plots
gg_maplot <- plot_ma(exp, data_type='tsr') +
  scale_color_viridis_d() +
  theme_bw()

ggsave('MA_Plot.png', gg_maplot, device="png", path='Plots', width=12, height=8)

# Volcano Plot
gg_volcano <- plot_volcano(exp, data_type='tsr') +
  scale_color_viridis_d() +
  theme_bw()

ggsave('diffFeature_volcano.png', gg_volcano, device="png", path='Plots', width=12, height=8)


