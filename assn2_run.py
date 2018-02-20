#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 14:51:12 2018

@author: owenwhitley
"""

import assn2_classes as as2
import re
import numpy as np
import sys
from matplotlib import pyplot as plt
import os

def test_assignment2(fname, input_dir, mode, max_iter):
    
    
    
    fileobj = open(os.path.join(input_dir, fname))
    seq_list = []
    for line in fileobj:
        
        if '#' not in line:
            seq_string = re.findall('[ACGU]*', line)[0]
            seq_list.append(seq_string)
            
    
    k = 6
    o = as2.offset_list()
    S = as2.seq_set()
    for i in range(len(seq_list)):    
        seq_i = as2.nuc_seq(seq_list[i])
        S.add_seq(nuc_seq = seq_i)
        o.add_offset(seq_id = seq_i.get_seq_id(), offset = 1)
    
    # be is an array holding probabilities of seeing A,C,G,U, respectively, by random chance    
    b = np.array([0.25, 0.25, 0.25, 0.25])
    nucleotides = ('A','C','G','U')
    
    #### Get results for test case using one of 4 algorithms ####
    
    if mode == 'hardOOPs':
        print('getting results for first test case')
        ## hardOOPs
        print('hardOOPs')
        (EMO, M, o_output, EMO_list, o_list) = as2.hardOOPs(S, o, b, k, nucleotides, 
                                                            delta = 0.01, max_iter = max_iter)
        print('========= RESULTS =======')
        print('log E(M,o)')
        print(EMO)
        print('Matrix M: rows (A,C,G,U), columns k = 1:6')
        print(np.transpose(M))
        print('offsets')
        print(o_output.offsets)
    elif mode == 'hardZOOPs':
        ## hardZOOPs
        print('getting results for first test case')
        print('hardZOOPs')
        (EMO, M, o_output, EMO_list, o_list) = as2.hardOOPs(S, o, b, k, nucleotides, 
                                                        delta = 0.01, max_iter = max_iter, mode = 'hardZOOPs')
        print('========= RESULTS =======')
        print('log E(M,o)')
        print(EMO)
        print('Matrix M: rows (A,C,G,U), columns k = 1:6')
        print(np.transpose(M))
        print('offsets')
        print(o_output.offsets)
    elif mode == 'gibbZOOPs':
    ## gibbZOOPs (Part c)
        print('getting results for first test case')
        print('gibbZOOPs')
        (EMO, M, o_output, EMO_list, o_list) = as2.hardOOPs(S, o, b, k, nucleotides, 
                                                        delta = 0.01, max_iter = max_iter, mode = 'gibbZOOPs')
        ## Get location of highest log EMO
        EMO_ind = np.in1d(EMO_list, EMO)
        loc_EMO = np.array(range(0,len(EMO_list)))[EMO_ind]
        if len(loc_EMO) > 1:
            loc_EMO = loc_EMO[0]
        print('========= RESULTS =======')
        print('log E(M,o)')
        print(EMO)
        print('Matrix M: rows (A,C,G,U), columns k = 1:6')
        print(np.transpose(M))
        print('offsets')
        print(o_output.offsets)
        print('iteration location of highest log(E(M,o))')
        print(loc_EMO)
        
    elif mode == 'sim_anneal':
        ## Simulated Annealing
        print('getting results for first test case')
        print('simulated annealing')
        (EMO, M, o_output, EMO_list, o_list) = as2.hardOOPs(S, o, b, k, nucleotides, 
                                                        delta = 0.01, max_iter = max_iter, mode = 'sim_anneal')
        print('========= RESULTS =======')
        print('log E(M,o)')
        print(EMO)
        print('Matrix M: rows (A,C,G,U), columns k = 1:6')
        print(np.transpose(M))
        print('offsets')
        print(o_output.offsets)
    else:
        raise ValueError('mode must be one of following: hardOOPs, hardZOOPs, gibbZOOPs, or sim_anneal')
        
    fileobj.close()
    plot_file_prefix = fname
    plot_fname = plot_file_prefix + '_' +  mode +  '_log_E_v_iteration.png'
    log_E = EMO_list
    iteration = range(0,len(EMO_list))
    p = plt.figure()
    ax1 = p.add_subplot(111)
    ax1.plot( iteration, log_E, alpha=0.5)
    ax1.set_xlabel("iteration", fontsize=15)
    ax1.set_ylabel("log E(M, o)", fontsize=15)
    ax1.set_title(plot_fname)
    plt.savefig(plot_fname)
    del p
    
    return((EMO, M, o_output, EMO_list, o_list))
        
    
################################################################################

if __name__ == '__main__':
    test_assignment2(fname = sys.argv[1], mode = sys.argv[3], max_iter = int(sys.argv[4]), input_dir = sys.argv[2])
    
    
