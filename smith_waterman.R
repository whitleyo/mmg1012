#### SMITH WATERMAN FUNCTION #####

smith_wat <- function(seq1, #character vector of length 1
                      seq2, # character vector of length 2
                      BLOSUM,
                      gap_init_pen = 1.0,
                      gap_extend_pen = 0.5) {
  if (!is.character(seq1) || (length('seq1') != 1)){
    stop()
  } else if (!is.character(seq2) || (length('seq2') != 1)) {
    stop()
  } else if(!is.matrix(BLOSUM)) {
    stop()
  }
  
  AA_set <- rownames(BLOSUM) #amino acids in BLOSUMmatrix
  seq1 <- strsplit(seq1, split = '') # make seq1 indexable character vector
  seq2 <- strsplit(seq2, split = '') # do same for seq2
  valid_aas1 <- all(seq1 %in% AA_set)
  valid_aas2 <- all(seq2 %in% AA_set)
  if (!(valid_aas1 && valid_aas2)) {
    print('valid aas1')
    print(seq1)
    print(valid_aas1)
    print('valid aas2')
    print(seq2)
    print(valid_aas2)
    stop('not all characters in seq 1 or seq 2 valid')
  }
  
  ## Initialize Matrix H
  n <- length(seq1) + 1
  m <- length(seq2) + 1
  H <- matrix(nrow = n, ncol = m)
  H[1:n,] <- 0 # set first column to 0
  H[,1:m] <- 0 # set first row to 0
  rownames(H) <- c('', seq1)
  colnames(H) <- c('', seq2)
  for (I in 1:(n-1)) {
    i <- I + 1 #set index i
    for (J in 1:(m-1)) {
      j <- J + 1 # set index j
      aa1 <- seq1[I]
      aa2 <- seq2[J]
      aa1_ind <- which(AA_set == aa1)
      aa2_ind <- which(AA_set == aa2)
      # Find H[i-1,j-1] + s(ai,bj)
      option1 <- H[i-1,j-1] + BLOSUM[aa1_ind, aa2_ind]
      opt2_vect <- numeric()
      # go along row from i-1 to H[1,j] (zero row) to add gapped scores
      for (k in seq(i-1, 1, -1)) {
        row_gap_pen <- gap_init_pen + gap_extend_pen*(i - k)
        score_k <- H[k,j] - row_gap_pen
        opt2_vect <- c(opt2_vect, score_k)
      }
      option2 <- max(opt2_vect) # take max row gapped score
      # go along col from j-1 to H[i,1] (zero column) to add gapped scores
      opt3_vect <- numeric()
      # go along row from i-1 to H[1,j] (zero row) to add gapped scores
      for (l in seq(j-1, 1, -1)) {
        col_gap_penalty <- gap_init_pen + gap_extend_pen*(j - l)
        score_k <- H[i,l] - col_gap_penalty
        opt3_vect <- c(opt2_vect, score_k)
      }
      option3 <- max(opt3_vect) # take max col gapped score
      
      #assign to H[i,j] the max value of options 1-3
      H[i,j] <- max(c(option1,option2,option3))
    }
  }
  # end for loop, H has been constructed
  
  ## TRACEBACK
  # Locate the max value in matrix
  end_pts <- which(H == max(H), arr.ind = TRUE)
  # in case there are multiple solutions
  solutions <- list()
  
  for (p in nrow(end_pts)) {
    i <- end_pts[p,1]
    j <- end_pts[p,2]
    out_seq1 <- character()
    out_seq2 <- character()
    cont <- TRUE
    while(cont) {
      # end condition is if 'max next' is 0
      opt1 <- H[i - 1, j - 1] # if max, then at i and j, perfect alignment
      opt2 <- H[i - 1, j] # if max, then gap in seq2 at this position
      opt3 <- H[i, j - 1] # if max, then gap in seq1 at this position
      opts <- c(opt1, opt2, opt3)
      choice <- which(opts %in% max(opts))
      if (choice == 1) {
        seq1_append <- rownames(H)[i]
        seq2_append <- colnames(H)[j]
        i <- i - 1
        j <- j - 1
      } else if (choice == 2) {
        seq1_append <- rownames(H)[i]
        seq2_append < '-'
        i <- i - 1
      } else if (choice == 3) {
        seq1_append <- '-'
        seq2_append <- colnames(H)[j]
        j <- j -1
      } else {
        stop('for some reason choice not made')
      }
      out_seq1 <- c(out_seq1, seq1_append)
      out_seq2 <- c(out_seq2, seq2_append)
      # Terminate if the maximum of options is less than or equal to 0, or if we've hi i = 0 or j = 0
      if (max(opts) <= 0) {
        cont <- FALSE
      } else if ((i == 0) || (j == 0)) {
        cont <- FALSE
      }
    } # end while loop
    out_seq1 <- paste(rev(out_seq1), collapse = '')
    out_seq2 <- paste(rev(out_seq2), collapse = '')
    solutions[[p]] <- list(out_seq1, out_seq2)
  }
}