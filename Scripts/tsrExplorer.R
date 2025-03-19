library('TSRexploreR')
library('ggplot2')


# load annotation (gtf) and assembly (fasta)
annotation <- read.table('/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/gencode.v31.basic.annotation.gtf', sep='\t')
assembly <- '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/genomes/hg38/fa/ebi/GRCh38.primary_assembly.genome.fa.gz'



# load samples
samples <- data.frame(sample_name = c('EDO_1','EDO_2','EDO_3','PW1_1','PW1_2','PW1_3'),
                      file_1 = c('/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/DexSeq_counts/EDO_1.dexeq.counts',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/DexSeq_counts/EDO_2.dexeq.counts',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/DexSeq_counts/EDO_3.dexeq.counts'
                                 ),
                      file_2 = c('/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/DexSeq_counts/PW1_1.dexeq.counts',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/DexSeq_counts/PW1_2.dexeq.counts',
                                 '/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/DexSeq_counts/PW1_3.dexeq.counts'
                                 ),
                      conditions = c(rep('EDO',3), rep('PW1', 3))
                      )

exp <- theDATA %>% 
  tsr_explorer(
    genome_annotation = annotation,
    genome_assembly = assembly,
    sample_sheet = samples
  ) %>% 
  format_counts(data_type = 'tss')

# cluster TSSs and merge replicates
exp <- tss_clustering(exp, threshold=3, n_samples=1) %>% 
  merge_samples(data_type='tss', merge_groups='condition') %>% 
  merge_samples(data_type='tsr', merge_groups='condition')

# calculate shifting
exp <- tss_shift(
  exp,
  sample_1=c(TSS='EDO', TSR='EDO'),
  sample_2=c(TSS='PW1', TSR='PW1'),
  comparison_name = 'EDO_vs_PW1',
  max_distance = 100, min_threshold = 10, n_resamples = 1000L
)

# plot shift scores ranked
plot_shift_rank(exp) +
  scale_fill_viridis_c() +
  theme_bw() +
  geom_hline(yintercept = 0)

# plot upstream and downstream shifts
plot_shift_count(exp) +
  scale_fill_viridis_d() +
  theme_bw()





