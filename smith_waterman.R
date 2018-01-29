#### SMITH WATERMAN FUNCTION #####

smith_wat <- function(seq1, #character vector of length 1
                      seq2, # character vector of length 2
                      BLOSUM,
                      gap_penalty = 1.0,
                      rev_seq2 = TRUE) {
  if (!is.character(seq1) || (length('seq1') != 1)){
    stop()
  } else if (!is.character(seq2) || (length('seq2') != 1)) {
    stop()
  } else if(!is.matrix(BLOSUM)) {
    stop()
  }
  
  base_set <- rownames(BLOSUM) # nucleotide basese in BLOSUMmatrix
  seq1 <- strsplit(seq1, split = '')[[1]] # make seq1 indexable character vector
  seq2 <- strsplit(seq2, split = '')[[1]] # do same for seq2
  if( rev_seq2) {
    seq2 <- rev(seq2)
  }
  valid_bases1 <- all(seq1 %in% base_set)
  valid_bases2 <- all(seq2 %in% base_set)
  if (!(valid_bases1 && valid_bases2)) {
    print('base set for BLOSUM')
    print(base_set)
    print('valid bases1')
    print(unique(seq1))
    print(valid_bases1)
    print('valid bases2')
    print(unique(seq2))
    print(valid_bases2)
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
      base1 <- seq1[I]
      base2 <- seq2[J]
      base1_ind <- which(base_set == base1)
      base2_ind <- which(base_set == base2)
      # Find H[i-1,j-1] + s(ai,bj)
      option1 <- H[i-1,j-1] + BLOSUM[base1_ind, base2_ind]
      # Find H[i - 1, j] - gap_penalty
      option2 <- H[i - 1,j] - gap_penalty
      # Find H[i, j - 1] - gap_penalty
      option3 <- H[i,j - 1] - gap_penalty
      
      #assign to H[i,j] the max value of options 1-3
      H[i,j] <- max(c(option1,option2,option3, 0))
    }
  }
  # end for loop, H has been constructed
  
  ## TRACEBACK
  # Locate the max value in matrix
  end_pts <- which(H == max(H), arr.ind = TRUE)
  # in case there are multiple solutions
  solutions <- list()
  
  gap_count <- 0
  for (p in 1:nrow(end_pts)) {
    i <- end_pts[p,1]
    j <- end_pts[p,2]
    end_s1 <- i-1
    end_s2 <- j-1
    out_seq1 <- character()
    out_seq2 <- character()
    cont <- TRUE
    while(cont) {
      # end condition is if 'max next' is 0
      base1 <- rownames(H)[i]
      base2 <- colnames(H)[j]
      base1_ind <- which(base_set == base1)
      base2_ind <- which(base_set == base2)
      # Find H[i-1,j-1] + s(ai,bj)
      opt1 <- H[i-1,j-1] + BLOSUM[base1_ind, base2_ind] # if max, then at i and j, perfect alignment
      opt2 <- H[i - 1, j] - gap_penalty # if max, then gap in seq2 at this position
      opt3 <- H[i, j - 1] - gap_penalty # if max, then gap in seq1 at this position
      opts <- c(opt1, opt2, opt3, 0)
      choice <- which(opts %in% max(opts))
      if (length(choice) > 1) {
        choice <- choice[1] # if more than one choice satisfies max score, take first option (i.e. direct alignment)
      }
      
      if (choice == 1) { #if we have a 'match' for given pair
        if(gap_count > 0) {
          seq1_append <- c(seq1_append, rownames(H)[i])
          seq2_append <- c(seq2_append, colnames(H)[j])
        } else if(gap_count == 0) {
          seq1_append <- rownames(H)[i]
          seq2_append <- colnames(H)[j]
        } else {
          stop('this should not be possible')
        }
        gap_count <- 0
        i <- i - 1 #back increment i and j
        j <- j - 1
      } else if (choice == 2) { # else if we have a gap in seq2 (columns)
        if(gap_count > 0) {
          seq1_append <- c(seq1_append, rownames(H)[i])
          seq2_append <- c(seq2_append, '-')
        } else if(gap_count == 0) {
          seq1_append <- rownames(H)[i]
          seq2_append <- '-'
        } else {
          stop('this should not be possible')
        }
        i <- i - 1
        gap_count <- gap_count + 1
      } else if (choice == 3) { # we have a gap in sequence 1
        if(gap_count > 0) {
          seq1_append <- c(seq1_append, '-')
          seq2_append <- c(seq2_append, colnames(H)[j])
        } else if(gap_count == 0) {
          seq1_append <- '-'
          seq2_append <- colnames(H)[j]
        } else {
          stop('this should not be possible')
        }
        j <- j -1
        gap_count <- gap_count + 1
      } else {
        stop('for some reason choice not made')
      }
      # Terminate if the maximum of options is less than or equal to 0, or if we've hi i = 0 or j = 0
      if (max(opts) <= 0) {
        cont <- FALSE
      } else if ((i == 0) || (j == 0)) {
        cont <- FALSE
      } else if (H[i,j] == 0) {
        cont <- FALSE
      }
      if (gap_count == 0){
        out_seq1 <- c(out_seq1, seq1_append)
        out_seq2 <- c(out_seq2, seq2_append)
        }
      

    } # end while loop
    if ((length(out_seq1) > 0) && (length(out_seq2) > 0)) {
      ali_seq <- rep('|', length(out_seq1))
      ali_seq[out_seq1 == '-'] <- ' ' 
      ali_seq[out_seq2 == '-'] <- ' '
      ali_seq <- paste(rev(ali_seq), collapse = '')
      ali_seq <- paste('  ', ali_seq, '  ')
      start_s1 <- i-1
      out_seq1 <- paste(rev(out_seq1), collapse = '')
      out_seq2 <- paste(rev(out_seq2), collapse = '')
      out_seq1 <- paste('5\'', out_seq1, '3\'')
      if(rev_seq2) {
        start_s2 <- length(seq2) - end_s2 + 1
        end_s2 <- length(seq2) - (j- 1) + 1 
        out_seq2 <- paste('3\'', out_seq2, '5\'')
      } else {
        out_seq2 <- paste('5\'', out_seq2, '3\'')
      }
      solutions[[p]] <- list( seq1_ali = out_seq1, alignment = ali_seq, seq2_ali = out_seq2, score = max(H), location = end_pts[p,],
                              seq1_start_end = c(start_s1, end_s1), seq2_start_end = c(start_s2, end_s2))
    } else {
      solutions[[p]] <- list(score = max(H), location = end_pts[p,])
    }
    
    
    
  } # end larger for loop
  return(solutions)
}

print_solutions <- function(solutions) {
  print(paste('MAX SCORE:', as.character(solutions[[1]]$score)))
  for (i in 1:length(solutions)) {
    solution_i <- solutions[[i]]
    print(paste('SOLUTION:', as.character(i)))
    print(paste('seq1_ali', solution_i$seq1_ali))
    print(paste('align   ', solution_i$alignment))
    print(paste('seq2_ali', solution_i$seq2_ali))
    print('seq1 position')
    print(solution_i$seq1_start_end)
    print('seq2 position')
    print(solution_i$seq2_start_end)
  }
}