#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 12:23:21 2018

@author: owen
"""
import subprocess
import assn2_run as AR
import os
### RUN FIRST TEST CASE ###
system_call = True
fname = 'example_input'
input_dir = 'input_dir'
output_file_prefix = 'example_input_'
if system_call == True:
    
    for mode in ['hardOOPs', 'hardZOOPs', 'gibbZOOPs', 'sim_anneal']:
        if mode == 'gibbZOOPs':
            max_iter = 150
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
    

#AR.test_assignment2(fname = 'example_input.txt', mode = 'hardOOPs')
#AR.test_assignment2(fname = 'example_input.txt', mode = 'hardZOOPs')
#AR.test_assignment2(fname = 'example_input.txt', mode = 'gibbZOOPs')
#AR.test_assignment2(fname = 'example_input.txt', mode = 'sim_anneal')