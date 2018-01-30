#### SMITH WATERMAN FUNCTION #####

# trace back function recursively constructs sequences and adds them to a list of solutions
# basically, the base case case is a situation where no branches are occured in traceback, and the termination condition of
# having a '0' neighbor occurs
# if a branch is encountered, for each branch direction, this function is called within itself, starting
# a new frame for trace bcak

# function is meant to be used within smith_wat function

trace_back <- function(inp1, # character vector, used to seed the constructed sequence for seq1 alignment
                       inp2, # same for sequence 2 alignment
                       seq1, # actual sequence for seq1
                       seq2, # same for seq2
                       H, # matrix constructed in first part of smith waterman algorithm
                       BLOSUM, # BLOSUM matrix (pairwise scoring matrix)
                       i, # starting position in H
                       j, # end position in H
                       gap_penalty) {
  
  cont <- TRUE
  gap_count <- 0
  inp1_append <- ''
  inp2_append <- ''
  while(cont) {
    # check neighbors
    base_ind1 <- seq1[i-1]
    base_ind2 <- seq2[j-1]
    option1 <- H[i-1, j - 1] + BLOSUM[base_ind1, base_ind2]
    option2 <- H[i -1, j] - gap_penalty #here we have a gap in seq 2
    option3 <- H[i , j - 1] - gap_penalty # gap in seq 1
    choices <- c(option1, option2, option3, 0)
    choice <- which(choices %in% max(choices))
      if (length(choice) == 1) { # if no branchpoints, append appropriate sequence characters
          if (choice == 1) {
            inp1 <- paste0 (inp1, base_ind1)
            inp2 <- paste0(inp2, base_ind2)
            i <- i - 1
            j <- j - 1
          } else if (choice == 2) {
            inp1_append <- paste0(inp1_append, base_ind1)
            inp2_append <- paste0(inp2_append, '-')
            i <- i - 1
            j <- j
            gap_count <- gap_count + 1
          } else if (choice == 3) {
            inp1_append <- paste0(inp1, '-')
            inp2_append <- paste0(inp2, base_ind2)
            i <- i
            j <- j - 1
            gap_count <- gap_count + 1
          } else if ((choice == 4) || (H[i,j] == 0)) {
            cont <- FALSE
            break
          } else {
            stop()
          }
        if (H[i,j] == 0) {
          cont <- FALSE
        }
        }
        else if (length(choice) > 1) { # if we have a branch point, then we recurse back into our traceback function and return a character vector of solutions
          solution <- character() # character vector to contain all solutions found in branches
          # first, append sequence, or make something to be appended
          for (n in 1:length(choice)) {
            if (choice[n] == 1) {
              inp1 <- paste0 (inp1, base_ind1)
              inp2 <- paste0(inp2, base_ind2)
              ri <- i - 1
              rj <- j - 1
            } else if (choice[n] == 2) {
              inp1 <- paste0(inp1, base_ind1)
              inp2 <- paste0(inp2, '-')
              ri <- i - 1
              rj <- j
            } else if (choice[n] == 3) {
              inp1 <- paste0(inp1, '-')
              inp2 <- paste0(inp2, base_ind2)
              ri <- i
              rj <- j -1
            } else if ((choice[n] == 4) || (H[i,j] == 0)) {
              cont <- FALSE
            } else {
              stop()
            }
            solution_n <- trace_back(inp1, inp2, seq1, seq2, H, BLOSUM, i = ri, j = rj, gap_penalty = gap_penalty)
            solution <- c(solution, solution_n)
          }
        return(solution)
          break
        }
      } # end while loop
  start_pt <- c(as.character(i), as.character(j))
  inp1 <- paste0(inp1, start_pt[1])
  inp2 <- paste0(inp2, start_pt[2])
  solution <- c(inp1, inp2)
  return(solution)
  } # end function

##### SMITH WATERMAN FUNCTION ##########

# implements smith waterman algorithm with max gap penalty of 1, considering gaps of length 1 for penalty
# by default, second sequence is reversed in function
smith_wat <- function(seq1, #character vector of length 1
                      seq2, # character vector of length 2
                      BLOSUM, # pairwise scorirng matrix
                      gap_penalty = 1.0,
                      rev_seq2 = TRUE ) { # reverse sequence? default true
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
  # Find all maxima in H, and use recursive traceback function to find all solutions
  # Locate the max value in matrix
  end_pts <- which(H == max(H), arr.ind = TRUE)
  # in case there are multiple solutions
  solutions <- list()
  
  solution_list <- list()
  solution_set <- character()
  for (m in 1:nrow(end_pts)) {
    i <- end_pts[m,1]
    j <- end_pts [m,2]
    inp1 <- paste0('seq1_', as.character(i - 1))
    inp2 <- paste0('seq2_', as.character(j - 1))
    solution_set_i <- trace_back(inp1, inp2, seq1, seq2, H, BLOSUM, i, j, gap_penalty)
    solution_set <- c(solution_set, solution_set_i)
  }
  for (m in seq(2, length(solution_set), 2)) {
    # slice and dice things with regular expressions to extract aligned sequences and their positions
    # note that every two sequences in solution_set is an aligned pair
    solution_list_i <- list(seq1_ali = character(), seq2_ali = character(), alignment = character(), ends1 = numeric(), ends2 = numeric())
    seq1_ali <- solution_set[m-1]
    ends1 <- regmatches(seq1_ali, regexpr('_[0-9]*', seq1_ali))
    ends1 <- as.numeric(sub('_', '', ends1))
    seq1_ali <- sub('[a-z]*[0-9]*_[0-9]*', '',  seq1_ali)
    ends1 <- c(ends1, as.numeric(regmatches(seq1_ali, regexpr('[0-9]*$', seq1_ali))))
    seq1_ali <- gsub('[0-9]*', '', seq1_ali)
    seq1_ali <- gsub('[0-9]', '', seq1_ali)
    seq2_ali <- solution_set[m]
    ends2 <- regmatches(seq2_ali, regexpr('_[0-9]*', seq2_ali))
    ends2 <- as.numeric(sub('_', '', ends2))
    seq2_ali <- sub('[a-z]*[0-9]*_[0-9]*', '',  seq2_ali)
    ends2 <- c(ends2, as.numeric(regmatches(seq2_ali, regexpr('[0-9]*$', seq2_ali))))
    seq2_ali <- gsub('[0-9]*', '', seq2_ali)
    seq2_ali <- gsub('[0-9]', '', seq2_ali)
    if (rev_seq2) { # reverse positions of sequence to reflect forward sequence position
      ends2[1] <- length(seq2) - ends2[1] + 1
      ends2[2] <- length(seq2) - ends2[2] + 1
    }
    seq2_ali <- sub('[0-9]*', '', seq2_ali)
    if (length(seq1_ali) == length(seq2_ali)) {
      seq1_ali_split <- strsplit(seq1_ali, split = '')[[1]]
      seq2_ali_split <- strsplit(seq2_ali, split = '')[[1]]
      alignment <- rep('|', length(seq1_ali_split))
      alignment[seq1_ali_split == '-'] <- ' '
      alignment[seq2_ali_split == '-'] <- ' '
      seq1_ali <- paste(rev(seq1_ali_split), collapse = '')
      seq2_ali <- paste(rev(seq2_ali_split), collapse = '')
      alignment <- paste(rev(alignment), collapse = '')
    } else {
      stop()
    }
    solution_list_i$seq1_ali <- seq1_ali
    solution_list_i$seq2_ali <- seq2_ali
    solution_list_i$alignment <- alignment
    solution_list_i$ends1 <- rev(ends1)
    solution_list_i$ends2 <- ends2
    solution_list_i$score <- max(H)
    solution_list[[m/2]] <- solution_list_i
  }
  
  
  return(solution_list)
}

######### PRINT SOLUTIONS ###########
# take solutions object (list) from smith_wat function, and print summary of all alignments
print_solutions <- function(solutions) {
  print(paste('MAX SCORE:', as.character(solutions[[1]]$score)))
  for (i in 1:length(solutions)) {
    solution_i <- solutions[[i]]
    print(paste('SOLUTION:', as.character(i)))
    print(paste('seq1_ali', solution_i$seq1_ali))
    print(paste('align   ', solution_i$alignment))
    print(paste('seq2_ali', solution_i$seq2_ali))
    print('seq1 position')
    print(solution_i$ends1)
    print('seq2 position')
    print(solution_i$ends2)
  }
}