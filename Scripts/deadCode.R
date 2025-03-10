
ncbi_result <- entrez_search(db = "gene", term = "LOC729737[Gene Name] AND human[Organism]")

# If a result is found, fetch details
if (length(ncbi_result$ids) > 0) {
  gene_info <- entrez_summary(db = "gene", id = ncbi_result$ids[1])
  print(gene_info$strand)  # Strand information
} else {
  print("Gene not found in NCBI database")
}


####
#   Simulating random sequences to test complementary scores
####

query <- ase1
NUCLEOTIDES <- c('A','G','C','T')
REPS <- 20000

test_scores <- numeric(length=REPS)
test_alignment_list <- vector('list', length = REPS)

for(i in seq(1,REPS)){
  DNA_sim <- DNAString(paste(sample(NUCLEOTIDES, 300, replace=T), collapse=''))
  test_alignment <- pwalign::pairwiseAlignment(query, DNA_sim, 
                             type = 'global-local', 
                             substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
                               match = 2, mismatch = -2, baseOnly = TRUE),
                             gapOpening = -5, 
                             gapExtension = -3)
  test_scores[i] <- score(test_alignment)
  test_alignment_list[[i]] <- test_alignment
}

####
# Find overlap between databases
####

peakfile <- read_excel('data/peaksAnnoFULL.xlsx', sheet='SYS5_differential')

rmbase2 <- read_excel('data/RMBase_cytosolic.xlsx', skip=1)
rmbase3 <- read_excel('data/RMBase_Other.xlsx', skip=1)
rmbase4 <- read_excel('data/RMBase_HEK.xlsx', skip=1)
rmbase5 <- read_excel('data/RMBase_HeLa.xlsx', skip=1)
rmbase6 <- read_excel('data/RMBase_PA1.xlsx', skip=1)

rmbase_list <- list(rmbase2, rmbase3, rmbase4, rmbase5, rmbase6)

relevant_rows <- c()
relevant_binary <- c()

for(i in seq(1,dim(peakfile)[1])){
  test1 <- any(rmbase2$ModStart[rmbase2$ModChr==peakfile$seqnames[i]] > peakfile$start[i] & rmbase2$ModStart[rmbase2$ModChr==peakfile$seqnames[i]] < peakfile$end[i])
  test2 <- any(rmbase3$ModStart[rmbase3$ModChr==peakfile$seqnames[i]] > peakfile$start[i] & rmbase3$ModStart[rmbase3$ModChr==peakfile$seqnames[i]] < peakfile$end[i])
  test3 <- any(rmbase4$ModStart[rmbase4$ModChr==peakfile$seqnames[i]] > peakfile$start[i] & rmbase4$ModStart[rmbase4$ModChr==peakfile$seqnames[i]] < peakfile$end[i])
  test4 <- any(rmbase5$ModStart[rmbase5$ModChr==peakfile$seqnames[i]] > peakfile$start[i] & rmbase5$ModStart[rmbase5$ModChr==peakfile$seqnames[i]] < peakfile$end[i])
  test5 <- any(rmbase6$ModStart[rmbase6$ModChr==peakfile$seqnames[i]] > peakfile$start[i] & rmbase6$ModStart[rmbase6$ModChr==peakfile$seqnames[i]] < peakfile$end[i])
  if(any(test1,test2,test3,test4,test5)){
    relevant_rows <- append(relevant_rows, i)
    relevant_binary <- append(relevant_binary, 1)
  } else {
    relevant_binary <- append(relevant_binary, 0)
  }
}

# loop over the 5 methyl databases and extract rows in which methylation site falls within a peak

rmbase_crossed_compiled <- data.frame(chr = character(), start = numeric(), end = numeric(), snoRNA = character(), gene = character())

for(k in seq(1, length(rmbase_list))){
  relevant_rows <- c()
  for(i in seq(1, dim(rmbase_list[[k]])[1])){
    #rmbase_list[[k]][i,]
    if(any(peakfile$start[peakfile$seqnames==rmbase_list[[k]]$ModChr[i]] < rmbase_list[[k]]$ModStart[i] & 
       peakfile$end[peakfile$seqnames==rmbase_list[[k]]$ModChr[i]] > rmbase_list[[k]]$ModStart[i])){
      rmbase_crossed_compiled[nrow(rmbase_crossed_compiled)+1,] <- c(rmbase_list[[k]]$ModChr[i], rmbase_list[[k]]$ModStart[i], rmbase_list[[k]]$ModEnd[i], 
                                                                  rmbase_list[[k]]$`snoRNA List`[i], peakfile$SYMBOL[peakfile$seqnames==rmbase_list[[k]]$ModChr[i]][which(
                                                                    peakfile$start[peakfile$seqnames==rmbase_list[[k]]$ModChr[i]] < rmbase_list[[k]]$ModStart[i] & 
                                                                      peakfile$end[peakfile$seqnames==rmbase_list[[k]]$ModChr[i]] > rmbase_list[[k]]$ModStart[i])])
    }
  }
}

write.csv(rmbase_crossed_compiled, 'data/rmbase_methylBase_crossed.csv', row.names=F)

## temp code just to manually loop 
j=1
if(T){
  print(j)
  print(to_export[score_vector>100,][j,])
  print(alignment_list[score_vector>100][j])
  
  j <- j + 1
}











