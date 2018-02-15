#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 12:51:00 2018

@author: owen
"""

import numpy as np
import re

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
            if start_ind < 0:
                raise ValueError
            elif (end_ind > len(self.seq_list)):
                raise ValueError
        except ValueError:
            raise ValueError('offset must be between 1 and len(seq) - k + 1')
            
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
        if (current_id in self.ids):
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

### FUNCTIONS ####

def calc_M(S, o, b, k, nucleotides):
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
    
    seqs = S.get_all_seqs()
    M_mat = np.zeros((k, len(nucleotides)), dtype = int)
    for kk in range(0,k,1):
        W_mat_kk = np.zeros(shape = (i, len(b)), dtype = int)
        for i in range(0,len(seqs)):
            # get sequence, offset, and subsequence
            seq_i = seqs[i]
            seq_id_i = seq_i.get_seq_id()
            offset_i = o.get_offset(seq_id_i)
            subseq_i = seqs[i].get_subseq(k = k, offset = offset_i)
            nuc = subseq_i[kk]
            W_mat_kk[i,:] = np.in1d(nucleotides, nuc).astype(dtype = int) # set row i column n to 1 if matches our nucleotide
            
            if sum(W_mat_kk[i:]) == 0:
                print(nuc)
                print('nuc')
                print('nucleotide not found in nucleotide set provided')
                raise ValueError
        # calculate given row of M
        M_mat[kk] = W_mat_kk.sum(axis = 0) + b 
        M_mat[kk] = M_mat[kk]/sum(M_mat[kk])
    
    # check for any rows with     
    rowsums = M_mat.sum(axis = 1)
    if (sum(np.in1d(rowsums, np.ndarray(0, dtype = int))) > np.ndarray(0, dtype = int)):
        print(rowsums)
        print('following indices have 0 propabilities for anything')
        print(range(k)[np.in1d(rowsums, np.ndarray(0, dtype = int))])
        raise ValueError
        
    return(M_mat)
        
            
            
            
            
                
    
        
    
    