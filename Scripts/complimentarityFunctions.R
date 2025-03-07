


find_index_after_n_chars <- function(seq, n) {
  count <- 0
  
  # Convert string into individual characters
  chars <- unlist(strsplit(seq, ''))
  
  # Iterate over the characters
  for (i in seq_along(chars)) {
    if (chars[i] != '-') {
      count <- count + 1
    }
    if (count == n) {
      return(i)
    }
  }
  # Return NA if there are fewer than n non-hyphen characters
  return(NA)
}



get_weighted_score <- function(seq1, seq2, match=2, mismatch=-2, gapOpening=-5, gapExtension=-3){
  score <- 0
  for(i in seq(1:nchar(seq1))){
    seq1_char <- substr(seq1, i, i)
    seq2_char <- substr(seq2, i, i)
    if(seq1_char == seq2_char){
      score <- score + match
    } else if(seq1_char == '-' | seq2_char == '-'){
      score <- score - gapOpening
    }
    ####### FIGURE OUT HOW EXACTLY GAP OPENING + EXTENSION WORKS
  }
  
  
  return(score)
}


get_weighted_score('CGTCATTCTCATCGGAA', 'CCTGCCTCCCAGCAGCC')