#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 12:23:21 2018

@author: owen
"""
import subprocess
import os
from matplotlib import pyplot as plt
import assn2_classes as AC
### RUN MULTIPLE TEST CASES ###
system_call = True
input_dir = './input_dir'
input_files = os.listdir(input_dir)

for f in input_files: 
    
    fname = f
    output_file_prefix = f # set output file prefix to input filename. Okay to do since no '.txt' or '.tsv'
    
    print('running for ' + f)
    for mode in ['hardOOPs', 'hardZOOPs', 'gibbZOOPs', 'sim_anneal']:
        print(mode)
        if system_call == True:
            if mode == 'gibbZOOPs':
                max_iter = 150
            elif (mode == 'sim_anneal') and (f == 'sample_input_hard_1'):
                max_iter = 1000
            else:
                max_iter = 500
            outfile =  output_file_prefix + mode + '_out.txt'
            if os.path.exists(outfile): # remove output file if it already exists
                os.remove(outfile)
            command = 'python3 assn2_run.py ' + fname + ' ' + input_dir + ' ' + mode + ' ' + str(max_iter)
            # run code caputuring stdout
            fileobj = open(outfile, 'a')
            subprocess.run(command.split(), stdout = fileobj)
            fileobj.close()
            del fileobj
                
    
