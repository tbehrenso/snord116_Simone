library(dplyr)


# Load the data
hg19_chimeras <- read.csv("/Users/tbehr/Desktop/NOP56_IP.snoRNA.hg19.chimeras.csv")


# ---------------------------------------------------
# Fisher's Exact Test --> doesn't do an apt comparison

# Find rows containing SNORD116
snord116_df <- hg19_chimeras %>% filter(grepl("SNORD116", reference_snoRNA))

# From those, find those that also mention SNORD3
snord116_snord3_df <- snord116_df %>% filter(grepl("SNORD3", reference_snoRNA))

# Get counts
n_total_snord116 <- nrow(snord116_df)
n_snord116_snord3 <- nrow(snord116_snord3_df)

# Print
cat("Total SNORD116 chimeras:", n_total_snord116, "\n")
cat("Chimeras with SNORD116 and SNORD3:", n_snord116_snord3, "\n")

# Total SNORD3 chimeras (regardless of SNORD116)
snord3_df <- hg19_chimeras %>% filter(grepl("SNORD3", reference_snoRNA))
n_total_snord3 <- nrow(snord3_df)

# Define cell values
a <- n_snord116_snord3                           # SNORD116 and SNORD3
b <- n_total_snord116 - a                        # SNORD116 only
c <- n_total_snord3 - a                          # SNORD3 only
d <- nrow(hg19_chimeras) - (a + b + c)                      # Neither

# Build contingency matrix
contingency_matrix <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
colnames(contingency_matrix) <- c("SNORD3", "Not_SNORD3")
rownames(contingency_matrix) <- c("SNORD116", "Not_SNORD116")

print(contingency_matrix)

# Fisher's exact test
fisher_result <- fisher.test(contingency_matrix)

# Show results
cat("Fisher's exact test p-value:", fisher_result$p.value, "\n")
cat("Odds ratio:", fisher_result$estimate, "\n")

# -------------------------------------------------------------------------
# Alternative

snord3_df <- hg19_chimeras %>% filter(grepl("SNORD3[^\\d]", reference_snoRNA, perl = T))

snord3_interactions <- c()

for(chim in snord3_df$reference_snoRNA){
  chim_split <- unlist(strsplit(chim, split='\\|'))
  
  is_snord_ABCD <- grepl("SNORD3[A-Z]", chim_split)
  is_snord_at <- grepl("SNORD3at", chim_split)
  is_u3 <- grepl("^U3", chim_split)
  
  is_interesting_snord <- !(is_snord_ABCD | is_snord_at | is_u3)
  
  chim_interesting <- chim_split[is_interesting_snord]
  
  if(length(chim_interesting) > 0){                                             # NOTE: only one chimera had two different snords (SNODB2088 + SNORD91A)
    print(chim_interesting)                                                     #       The one SNORD91A is not excluded here
    snord3_interactions <- c(snord3_interactions, chim_interesting[1])
  }
}

table(snord3_interactions)






