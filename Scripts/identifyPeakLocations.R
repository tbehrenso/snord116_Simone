library(readxl)
library(tidyverse)
library(ggrepel)

SYS5_differential <- read_excel('peaksAnno.xlsx', sheet='SYS5_differential')

# SYS5_diff_annotation <- SYS5_differential$annotation
# 
# is_exon <- startsWith(SYS5_diff_annotation, c('Exon'))
# is_intron <- startsWith(SYS5_diff_annotation, c('Intron'))
# 
# SYS5_diff_annotation[is_exon] <- 'Exon'
# SYS5_diff_annotation[is_intron] <- 'Intron'
# 
# unique(SYS5_diff_annotation)

# change annotations to group all that represent promoters, exons, and introns
SYS5_differential <- SYS5_differential %>% 
  mutate(annotation_category = case_when(
    startsWith(annotation, 'Exon') & distanceToTSS >= 0 ~ 'Exon/Intron',
    startsWith(annotation, 'Intron') & distanceToTSS >= 0 ~ 'Exon/Intron',
    startsWith(annotation, 'Promoter') ~ 'Promoter',
    startsWith(annotation, "5'") ~ "5' UTR",
    startsWith(annotation, "3'") ~ "3' UTR",
    TRUE ~ 'Intergenic'
  ))

SYS5_diff_table <- as.data.frame(table(SYS5_differential$annotation_category))
# calculate and create labels for pie chart
SYS5_diff_table <- SYS5_diff_table %>% 
  mutate(percentage = Freq / sum(SYS5_diff_table$Freq)) %>% 
  mutate(percent_label = paste(Var1,'\n',scales::percent(percentage))) %>% 
  mutate(absolute_label = paste(Var1,'\n', Freq))

ggplot(SYS5_diff_table, aes(x="", y=Freq, fill=Var1,label = absolute_label)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  geom_label_repel(aes(), 
             size = 10,
             position = position_stack(vjust = .3),
             show.legend = FALSE) +
  scale_fill_discrete(name = 'Annotation') +
  ggtitle('Frequency of Annotation Types for SHSY Peaks') +
  scale_fill_brewer(palette="Set2")

