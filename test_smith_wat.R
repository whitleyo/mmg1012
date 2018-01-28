#### test smith waterman ####
source('smith_waterman.R')

# setup pairwise scoring matrix BLOSUM_mat
BLOSUM_mat <- matrix(nrow = 4, ncol = 4, data = rep(0, 16))
rownames(BLOSUM_mat) <- c('A', 'U', 'G', 'C')
colnames(BLOSUM_mat) <- rownames(BLOSUM_mat)
BLOSUM_mat['G','C']  <- 3
BLOSUM_mat['A', 'U'] <- 2
BLOSUM_mat['U', 'G'] <- 1
BLOSUM_mat <- BLOSUM_mat + t(BLOSUM_mat)
BLOSUM_mat[BLOSUM_mat == 0] <- -1

## TEST CASES ##

#Test case handout
print('Test case from Handout')
seq1 <- 'CUUUUUGCGGUCUGGGCUUGC'
seq2 <- 'GUACAGAGCACAACCCCAAGUGCAAAAAUGUAUUUUAUUCUAACACUGAUGCAGAAAGUUGGCAGUGAACCAUUU'
solutions_hdout <- smith_wat(seq1, seq2, BLOSUM_mat, gap_penalty = 1.0)
print_solutions(solutions_hdout)

print('======================')
print('Extra Test Cases')
#Test case 1
seq1 <- 'GGGGGAAAAAUUUUUCCCCCGGGG'
seq2 <- 'AAAAGGGGGAAAAUUUUUUUAAA'
print('Test case 1')
solutions1 <- smith_wat(seq1, seq2, BLOSUM_mat, gap_penalty = 1.0)
print_solutions(solutions1)

#Test case 2
print('Test case 2')
seq1 <- 'GGGGGAAAAAUUUUCCCCCGGGG'
seq2 <- 'AAAAGGGGGAAUUUUUUUAAA'
solutions2 <- smith_wat(seq1, seq2, BLOSUM_mat, gap_penalty = 1.0)
print_solutions(solutions2)
#Test case 3
print('Test case 3')
seq1 <- 'CAAGUGACGGUUGAAA'
seq2 <- 'UUUCAACCGCACUUG'
solutions3 <- smith_wat(seq1, seq2, BLOSUM_mat, gap_penalty = 1.0)
print_solutions(solutions3)
#Test case 4
print('Test case 4')
seq1 <- 'CAAGUGCAGUUGAAA'
seq2 <- 'UUUCAAGCACUUG'
solutions4 <- smith_wat(seq1, seq2, BLOSUM_mat, gap_penalty = 1.0)
print_solutions(solutions4)
#Test case 5
print('Test case 5')
seq1 <- 'CCGAUAGG'
seq2 <- 'CCUUUCGGCCAAUCGG'
solutions5 <- smith_wat(seq1, seq2, BLOSUM_mat, gap_penalty = 1.0)
print_solutions(solutions5)
#Test case 6
print('Test case 6')
seq1 <- 'UGAGGUAGUAGGUUGUAUA'
seq2 <- 'GCAAUGAUGCCUACCAAACAUUUCCAGACUUAACAUUUUGGUCUCUG'
solutions6 <- smith_wat(seq1, seq2, BLOSUM_mat, gap_penalty = 1.0)
print_solutions(solutions6)
#Test case 7
print('Test case 7')
seq1 <- 'GUGAAAUGUU'
seq2 <- 'AUUUCCAGGAAUUUAUUCCCCUUCAUAAUUUGUCUCAUUUCAUUUUAUUUCAUCCACUUGGUAGAUGAAGUCACG'
solutions7 <- smith_wat(seq1, seq2, BLOSUM_mat, gap_penalty = 1.0)
print_solutions(solutions7)
#Test case 8
print('Test case 8')
seq1 <- 'AAAGAAUUC'
seq2 <- 'CAUGAAUGAAGAUAGGUUGUAAACUGAAUGCUGUGAUAAUACUCUGUAUUCUUUAUGGAAAAUGUUGUCCUGU'
solutions8 <- smith_wat(seq1, seq2, BLOSUM_mat, gap_penalty = 1.0)
print_solutions(solutions8)