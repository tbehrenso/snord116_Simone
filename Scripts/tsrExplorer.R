library('TSRexploreR')
library('ggplot2')


# load annotation (gtf) and assembly (fasta)
annotation <- 11
assembly <- 22
  

# load samples

samples <- data.frame(sample_name = c('EDO_1','EDO_2','EDO_3','ND1_1','ND1_2','ND1_3','PW1_1','PW1_2','PW1_3'),
                      file_1 = c('FILES PATHS HERE'),
                      file_2 = c('FILES PATHS ERE MATE'),
                      conditions = c(rep('EDO',3), rep('ND1',3), rep('PW1', 3))
                      )

exp <- theDATA %>% 
  tsr_explorer(
    genome_annotation = annotation,
    genome_assembly = assembly,
    sample_sheet = samples
  ) %>% 
  format_counts(data_type = 'tss')










