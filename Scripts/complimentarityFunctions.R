


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



get_weighted_score <- function(seq1, seq2, match=2, mismatch=-2, gapOpening=-5, gapExtension=-3, weighting=T){
  score <- 0
  firstGap1 <- T
  firstGap2 <- T
  seq1_position <- 1 # keep track of seq1 position (assumes seq1 is the ASE, so the weighting can work properly)
  position_weighting <- c(0.3, 0.9, 0.95, 0.99, 1, 0.9, 0.95, 0.95, 0.99, 0.95, 0.9, 0.8, 0.7, 0.45, 0.3, 0.3, 0.3)
  if(!weighting){position_weighting <- (rep(1, 17))}
  
  # correct for unequal sequence lengths
  while(nchar(seq1) > nchar(seq2)){
    seq2 <- paste0(seq2, '-')
  }
  
  for(i in seq(1:nchar(seq1))){
    seq1_char <- substr(seq1, i, i)
    seq2_char <- substr(seq2, i, i)
    # if not a gap, reset the firstGapX counter
    if(seq1_char != '-'){
      firstGap1 <- T
    }
    if(seq2_char != '-'){
      firstGap2 <- T
    }
    # loop over possible values to assign relevant scores
    if(seq1_char == seq2_char){
      score <- score + match*position_weighting[seq1_position]
      seq1_position <- seq1_position + 1
    } else if(seq1_char == '-'){
      if(firstGap1){                                       # gaps make weighting difficult. Uses weight of previous position
        score <- score + (gapOpening + gapExtension)*position_weighting[seq1_position]
        firstGap1 <- F
      } else {
        score <- score + gapExtension*position_weighting[seq1_position]
      }
    } else if(seq2_char == '-'){
      if(firstGap2){
        score <- score + (gapOpening + gapExtension)*position_weighting[seq1_position]
        firstGap2 <- F
      } else {
        score <- score + gapExtension*position_weighting[seq1_position]
      }
      seq1_position <- seq1_position + 1
    } else if(seq1_char != seq2_char){
      score <- score + mismatch*(position_weighting[seq1_position])
      seq1_position <- seq1_position + 1
    }
    #print(score)
  }
  
  return(score)
}

get_weighted_score('CGTCATTCTCATCGGAA', 'CCTGCCTCCCAGCAGCC')

get_weighted_score(seq1,seq2, weighting=F)

