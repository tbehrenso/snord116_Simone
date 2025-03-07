


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

