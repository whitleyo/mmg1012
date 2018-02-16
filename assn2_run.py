#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 14:51:12 2018

@author: owenwhitley
"""

import assn2_classes as as2
import re
import numpy as np

fname = 'example_input.txt'
fileobj = open(fname)

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
    
b = np.array([0.25, 0.25, 0.25, 0.25])
nucleotides = ('A','C','G','U')

(EMO, M, o, EMO_list) = as2.hardOOPs(S, o, b, k, nucleotides, delta = 0.01, max_iter = 10^3)
