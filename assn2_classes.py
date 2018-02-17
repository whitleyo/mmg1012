#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 12:51:00 2018

@author: owen
"""

import numpy as np
import random as rand
################ Define Classes #############

class nuc_seq:
    id_list = []
    id_index = 0
    ''' contains a list of nucleotides. initialize based off of input string'''
    def __init__(self, inpstr = str):
        inpstr_list = list(inpstr)
        bases = ['A','C','G','U']
        in_bases = np.in1d(inpstr_list, bases)
        in_bases.astype(dtype = int)
        if (sum(in_bases) != len(in_bases)):
            raise ValueError('inpstr must contain only bases acgu')
        # setup the seq list
        # setup the sequence id
        self.seq_list = inpstr_list
        self.seq_id = nuc_seq.id_index
        # increment sequence id by 1
        nuc_seq.id_index += 1
        
    def get_seq(self):
        return(self.seq_list)
    
    def get_subseq(self, k = int, offset = int):
        # k: length of k-mer window for extracting subsequence
        # offset: positive integer indicating offset for sequence extracted
        
        # eg given seq.seq_list = ['a','b','c','d'],
        # seq.get_subseq(k = 2, offset = 3)
        # -> ['c','d']
        try:
            start_ind = offset - 1
            end_ind = offset + k - 1
            if start_ind < -1:
                raise ValueError
            elif (end_ind > len(self.seq_list)):
                raise ValueError
        except ValueError:
            raise ValueError('offset must be between 0 and len(seq) - k + 1')
            
        subseq = ''
        
        for i in range(start_ind, end_ind):
            subseq = subseq + self.seq_list[i]
        
        return(nuc_seq(inpstr = subseq))
        
    def get_seq_id(self):
        return(self.seq_id)
    
class seq_set:
    
    # seqs : set of sequences
    def __init__(self):
        self.seqs = []
        self.id_mapping = {}
    def add_seq(self, nuc_seq = nuc_seq):
        # add sequence and map its id to its index
        current_id = nuc_seq.get_seq_id()
        if (current_id in self.id_mapping.keys()):
            raise ValueError('seq_id must be unique')
        self.seqs.append(nuc_seq)
        self.id_mapping[current_id] = len(self.seqs) - 1
            
    def pop_seq(self, seq_id = int):
        ind = self.id_mapping[seq_id] # retrieve index of sequence
        x = self.seqs.pop(ind) # pop the sequence
        del(self.id_mapping[seq_id]) # remove the id to index mapping
        return(x)
         
    def get_seq(self, seq_id = int):
        ind = self.id_mapping[seq_id]
        x = self.seqs[ind]
        return(x)
        
    def get_mult_seqs(self, seq_ids = list):
        seqs_return = []
        
        for i in range(0, len(seq_ids)):
            seq_id_i = seq_id[i]
            if not (seq_id_i == int):
                raise TypeError('seq_ids must be list of sequence ids')
            ind = self.id_mapping[seq_id_i]
            x = self.seqs[ind]
            seqs_return.append(x)
            
        return(seqs_return)
    
    def get_all_seqs(self):
        return(self.seqs)
        
class offset_list:
    
    def __init__(self):
        self.offsets = []
        self.seq_2_offset = {}
        
    def add_offset(self, seq_id = int, offset = int):
        
        if (seq_id in self.seq_2_offset.keys()):
            raise ValueError('sequence id already in offset list')
        self.offsets.append(offset)
        self.seq_2_offset[seq_id] = len(self.offsets) - 1
        
    def change_offset(self, seq_id = int, offset = int):
        
        if (seq_id in self.seq_2_offset.keys()):
            ind = self.seq_2_offset[seq_id]
            self.offsets[ind] = offset
        else:
            raise ValueError
            
    def get_offset(self, seq_id):
        
        if (seq_id in self.seq_2_offset.keys()):
            ind = self.seq_2_offset[seq_id]
            return(self.offsets[ind])
        else:
            raise ValueError
            
    def pop_offset(self, seq_id):
        
        if (seq_id in self.seq_2_offset.keys()):
            ind = self.seq_2_offset[seq_id]
            x = self.offsets.pop(ind)
            del self.seq_2_offset[seq_id]
            return(x)
        else:
            raise ValueError
##########################################################################
nucleotides = {}
nucleotides['A'] = 1
nucleotides['C'] = 2
nucleotides['G'] = 3
nucleotides['U'] = 4
### FUNCTIONS ####

def calc_M(S, o, b, k, nucleotides, mode = 'hardOOPs'):
    # Calculate numpy matrix M given a set of sequences, an offset_list, b (a numpy array), and k (an integer)
    # output should be a k x len(nucleotides) matrix containing estimated frequencies
    # of nucleotides at poition kk in k
    try:
        if type(S) != seq_set:
            raise TypeError
        if type(o) != offset_list:
            raise TypeError
        if type(b) != np.ndarray:
            raise TypeError
        if type(k) != int:
            raise TypeError
    except TypeError:
        raise TypeError('S must be of class seq_set, o of type offset_list, b of type np.ndarray, k of type int')
        
    if not (mode in ['hardOOPs', 'hardZOOPs', 'gibbZOOPs', 'sim_anneal']):
        raise ValueError('Mode must be one of following: hardOOPs, hardZOOPs, gibbZOOPs, sim_anneal')
    
    seqs = S.get_all_seqs()
    # make k x nucleotides matrix
    M_mat = np.zeros((k, len(nucleotides)), dtype = float)
    M_mat = np.matrix(M_mat, dtype = float)
    for kk in range(0,k,1):
        W_mat_kk = np.zeros(shape = (len(seqs), len(b)), dtype = int)
        for i in range(0,len(seqs)):
            # get sequence, offset, and subsequence
            seq_i = seqs[i]
            seq_id_i = seq_i.get_seq_id()
            offset_i = o.get_offset(seq_id_i)
            # if offset is 0, either skip finding matches for seq_i position kk,
            # or raise an error if running hardoops. The error should not be raised
            # in case of hardoops if the offsets are set correctly in the E-step.
            if offset_i == 0:
                if mode == 'hardOOPs':
                    raise ValueError('hardOOPs mode chosen, cannot have offsets of 0')
                else:
                    continue
                
            subseq_i = seqs[i].get_subseq(k = k, offset = offset_i)
            subseq_i_list = subseq_i.get_seq()
            nuc = subseq_i_list[kk]
            W_mat_kk[i,:] = np.in1d(nucleotides, nuc).astype(dtype = int) # set row i column n to 1 if matches our nucleotide
            
            if np.sum(W_mat_kk[i:]) == 0:
                if mode == 'hardOOPs':
                    print(nuc)
                    print('nuc')
                    print('nucleotide not found in nucleotide set provided')
                    raise ValueError
        # calculate given row of M
        M_mat[kk,:] = W_mat_kk.sum(axis = 0) + b 
        M_mat[kk,:] = M_mat[kk,:]/np.sum(M_mat[kk,:])
        del(W_mat_kk)
        
    
    # check for any rows with     
    rowsums = M_mat.sum(axis = 1)
    if (sum(np.in1d(rowsums, np.ndarray(0, dtype = int))) < np.ndarray(1, dtype = int)):
        print(rowsums)
        print('following indices have < 1 propabilities')
        print(range(k)[np.in1d(rowsums, np.ndarray(0, dtype = int))])
        raise ValueError
    M = M_mat    
    return(M)

### END FUNCTION ###
    
def calc_log_LR(nuc_seq, M, b, nucleotides):
    ## Calculate log likelihood ratio for sequence given matrix M and background frequenceies b
    if type(b) != np.ndarray:
        raise TypeError
    if not (type(M) == np.matrix):
        raise TypeError
    seq = nuc_seq.get_seq()
    dim_M = np.shape(M)
    equal_length = dim_M[0] == len(seq)
    nucs4 = (dim_M[1] == len(nucleotides)) & (len(b) == len(nucleotides))
    if not equal_length:
        raise ValueError('M must have number of rows same as length of provided seq')
    if not nucs4:
        raise ValueError('M must have 4 columns, b must have length of 4')
    # initialize log_LR
    log_LR = 0
    
    for k in range(0, len(seq)):
        
        current_nuc = seq[k] # set current nucleotide
        nuc_ind = np.in1d(nucleotides, current_nuc)
        num_matches = sum(nuc_ind.astype(dtype = bool))
        
        if num_matches != 1:
            print('num_matches in nucleotides at position ' + str(k))
            print(str(num_matches))
            raise ValueError('Most have one match in provided nucleotides per position in sequence')
        # get log ratio of Pc/Qc, add to log_LR
        Pc = M[k, nuc_ind]
        Qc = b[nuc_ind]
        log_LR += np.log(Pc/Qc)
    
    return(log_LR)
    

        
        
    
def calc_E(S, M, b, k, nucleotides, mode = 'hardOOPs', T = None):
    # ACGU -> 1,2,3,4 for nucleotides
    if type(S) != seq_set:
        raise TypeError
    else:
        seqs = S.get_all_seqs()
    if type(b) != np.ndarray:
        raise TypeError
    if type(k) != int:
        raise TypeError
    ## Check M argument
    is_matrix = (type(M) == np.matrix)
    M_shape = (np.shape(M) == (k, len(nucleotides)))
    if not (is_matrix and M_shape):
        raise TypeError
    if not (mode in ['hardOOPs', 'hardZOOPs', 'gibbZOOPs', 'sim_anneal']):
        raise ValueError('Mode must be one of following: hardOOPs, hardZOOPs, gibbZOOPs, sim_anneal')
    
    if mode == 'sim_anneal':
        if type(T) != float:
            raise TypeError('sim_anneal mode selected: T must be floating point number!')
    
    # Stores max likelihood ratios of sequences
    Max_LR_List = []
    offsets = offset_list()
    
    # set start of range for offsets to use.
    if (mode == 'hardOOPs') or (mode == 'hardZOOPs'):
        q = 1 
    elif (mode == 'gibbZOOPs') or (mode == 'sim_anneal'):
        q = 0
    else:
        print('Expect 1 of 4 modes to be satisfied')
        assert(False)
    
    for i in range(0,len(seqs)):
        
        seq_i = seqs[i]
        LR_list_i = []
        
        offset_range = np.arange(q,(len(seq_i.get_seq()) - k + 2)) # want q:len-k + 1, python requires an extra 1 to be added
        ## Calculate log likelihood ratios
        for j in offset_range:
            offset_i_j = j
            if offset_i_j > 0:
                subseq_i_j = seq_i.get_subseq(k = k, offset = offset_i_j)
                LR_list_i.append(calc_log_LR(subseq_i_j, M, b, nucleotides))
            elif offset_i_j == 0:
                LR_list_i.append(int(1))
            else:
                raise ValueError('offset_i_j should be >= 0')
        
        if mode == 'hardOOPs':
            ## take max likelihood ratio, pick the offset corresponding to it
            Max_LR_i = max(LR_list_i)
            Max_LR_ind = np.in1d(LR_list_i, Max_LR_i)
            num_max = sum(Max_LR_ind.astype(dtype = int)) # number of maximum LR satisfying positions
            if (num_max < 1):
                raise ValueError('num_max should be greater than or equal to 1')
            if Max_LR_i < 1:
                offset_i = int(0)
            else:
                offset_i = offset_range[Max_LR_ind]
                if len(offset_i) > 1:
                    print('WARNING: Offset for sequence ' + str(i) + ' seq_id ' + str(seq_i.get_seq_id()) +
                          ' has more than one optimal offset. Taking lowest of '+  str(len(offset_i)) + 
                          ' offsets' )
                    print('offsets')
                    print(offset_i)
                    offset_i = int(offset_i[0])
                
        elif mode == 'hardZOOPs':
            ## Calculate max log ratio, if < 1, set offset to 0
            # get max Log likelihood ratio
            Max_LR_i = max(LR_list_i)
            Max_LR_ind = np.in1d(LR_list_i, Max_LR_i)
            num_max = sum(Max_LR_ind.astype(dtype = int))
            if (num_max < 1):
                raise ValueError('num_max should be greater than or equal to 1')
            if Max_LR_i < 1: # if max likelihood ratio less than 
                Max_LR_i = 1
                offset_i = 0
            else:
                offset_i = offset_range[Max_LR_ind]
                if len(offset_i) > 1:
                    print('WARNING: Offset for sequence ' + str(i) + ' seq_id ' + str(seq_i.get_seq_id()) +
                          ' has more than one optimal offset. Taking lowest of '+  str(len(offset_i)) + 
                          ' offsets' )
                    print('offsets')
                    print(offset_i)
                    offset_i = int(offset_i[0])
                    
        elif (mode == 'gibbZOOPs') or (mode == 'sim_anneal'):
            ## Calculate g(p) as outlined in handout
            # LR_list_i holds the log likelihood ratios for each
            offset_probs = np.array(LR_list_i, dtype = float)
            if mode == 'sim_anneal': #if simulated annealing, moderate probabilit by temperature factor
                offset_probs = np.power(offset_probs, T)
            gp_vect = offset_probs/np.sum(offset_probs) # 
            ##sample offset_i from g[p] distribution
            # basic procecdure: pick a random num 0<x<1.
            # sum gp_vect[0:j] (cumulative probability up to j) > x, index j set to true
            # pick lowest index j, and offset j is selected
            pick_num = rand.randrange(0,1000)/1000
            
            for j in range(0,len(gp_vect)):
                
                if np.sum(gp_vect[0:j]) >= pick_num:
                    offset_i = j
                    break
                #below if statement a sanity check: if you haven't picked an index,
                # something's wrong because sum(gp_vect[0:len(gp_vect)]) should be 1
                if j == (len(gp_vect) - 1):
                    raise ValueError('j hit end of gp_vect, offset_i still not assigned')
        
            
            
        # add offset
        offset_i = int(offset_i)
        offsets.add_offset(seq_id = seq_i.get_seq_id(), offset = offset_i)
        Max_LR_List.append(Max_LR_i)
        
    ## end for loop
    EMO = sum(Max_LR_List)
    return(EMO, offsets)
    
##### END FUNCTION ###############
    
    
def hardOOPs(S, o, b, k, nucleotides, delta = 0.01, max_iter = 10^5, mode = 'hardOOPs'):
    
    if type(S) != seq_set:
        raise TypeError
    if type(o) != offset_list:
        raise TypeError
    if type(b) != np.ndarray:
        raise TypeError
    if type(k) != int:
        raise TypeError
    if not (mode in ['hardOOPs', 'hardZOOPs', 'gibbZOOPs', 'sim_anneal']):
        raise ValueError('Mode must be one of following: hardOOPs, hardZOOPs, gibbZOOPs, sim_anneal')

    # initialize empty containers for log EMO calculated values and determined offsets after each E step
    EMO_list = [] 
    o_list = []
    cont = True
    N = 0
    while cont == True:
        if mode == 'sim_anneal':
            T = N
        M = calc_M(S, o, b, k, nucleotides, mode)
        (EMO, o_new) = calc_E(S, M, b, k, nucleotides, mode, T) # calculate log EMO, determine new offsets
        o = o_new
        print('iteration')
        print(N)
        if N > 0:
            delta_N = EMO - last_EMO
            if delta_N < abs(delta):
                cont = False
        elif N >= max_iter:
            cont = False
        last_EMO = EMO
        EMO_list.append(EMO)
        o_list.append(o) 
        N += 1
    
    EMO_list = np.array(EMO_list)
    o_list = np.ndarray(o_list)
    if mode == 'gibbZOOPs':
        
        max_EMO = np.max(EMO_list)
        max_EMO_ind = np.in1d(EMO_list, max_EMO)
        max_offset = o_list[max_EMO_ind]
        if np.sum(max_EMO_ind.astype(dtype = int)) > 1:
            print('WARNING, more than 1 maximum reached, taking location of first maximum reached')
            max_offset = max_offset[0]
            
        o = max_offset
        M = calc_M(S, o, b, k, nucleotides, mode)
    return(EMO, M, o, EMO_list, o_list)
        
### END FUNCTION
    
    
    
        
            
            
            
            
                
    
        
    
    