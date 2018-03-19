##!/usr/bin/env python3
## -*- coding: utf-8 -*-
#"""
#Created on Sat Mar  3 19:31:58 2018
#
#@author: owen
#"""
import os
import re
import numpy as np
from scipy import stats
import pandas as pd
import warnings
assn3_files = os.listdir()

#
#### Check that all lines below header are QC+
##for file in assn3_files:
##    
##    file_match = re.search('\.hmap$', file)
##    print(file)
##    if not file_match:
##        
##        continue
##    print('match found for hmap file')
##    my_file = file
##    my_fileobj = open(my_file)
##    
##    N = 0
##    QC_count = 0
##    my_regex = 'QC\+'
##    for line in my_fileobj:
##        N += 1
##        
##        QC_match = re.search(my_regex, line)
##        
##        if QC_match:
##            QC_count += 1
##            
##    print('QC_count')            
##    print(QC_count)
##    print('number of lines in file')
##    print(N)
##    
##    if (N - QC_count) != 1:
##        assert False
##    
##    my_fileobj.close()
#    
#    
###############################################################
#
### OUTLINE ##
#
### Functions for part 1
### PART1:
###   Loading in data for parts 1a and 1b
###   Part 1a: 
#
################################################################################
### DEFINE FUNCTIONS TO RETRIEVE ALLELES< CALCULATE PROPORTIONS, Fst ###########
#    
#def retrieve_header(file):
#    
#    ## In first line, find index i with QC code, select index i + 1
#    ## as index
#    
#    N = 0
#    
#    fileobj = open(file)
#    for line in fileobj:
#    
#        N += 1
#        line_split = re.split('\s', line)
#        line_split = line_split[0:(len(line_split) - 1)]
#        
#        for i in range(0,len(line_split)):
#            
#            if re.match('QCcode', line_split[i]):
#                samples_start = i + 1
#                break
#            elif i == (len(line_split) - 1):
#                ## Stop program if cannot find QCcode
#                assert False
#                
#        if N > 0:
#            break
#        
#    sample_names = line_split[samples_start:len(line_split)]
#    ## check that all sample names have regex 'NA[0-9]*'
#    for i in range(0,len(sample_names)):
#        samp_name = sample_names[i]
#        if not (re.match('NA[0-9]*', samp_name)):
#            print(samp_name + ' does not fit regex NA[0-9]*. Occurs at element ' + str(samples_start + i) + ' in header')
#
#    return((sample_names, samples_start))
#    
#def retrieve_allele_pairs(file, allele, samples_start):
#    
#    ## Search each file line for the appropriate allele. When found,
#    ## return a list of allele pairs
#    
#    fileobj = open(file)
#    N = 0
#    
#    for line in fileobj:
#        N +=1
#        if re.match(allele, line):
#            
#            print('found match at line ' + str(N))
#            line_split = re.split('\s', line)
#            line_split = line_split[0:(len(line_split) - 1)]
#            
#            
#            
#            sample_alleles = line_split[samples_start:len(line_split)]
#            
#            removal_inds = []
#        
#            for n in range(0, len(sample_alleles)):
#                
#                samp = sample_alleles[n]
#                
#                    
##                if re.search('^[ATCGN]', samp):
##                    print('element' + str(n))
##                    print(samp)
##                    raise ValueError
#                if re.search('N', samp):
#                    print('element' + str(n))
#                    print(samp)
#                    print('excluding from allele pairs used in FST')
#                    removal_inds.append(n)
#                elif not re.match('[ATCG]{2}', samp):
#                    print('element' + str(n))
#                    print(samp)
#                    raise ValueError
#                    
#            sample_alleles = np.array(sample_alleles)
#            sample_alleles = sample_alleles[(np.logical_not(np.in1d(sample_alleles, sample_alleles[removal_inds])))]
#        
#            return(sample_alleles)
#    
#def allele_set_2_counts(allele_pairs):
#    
#    ## Find all unique alleles in the set of allele pairs, and return counts
#    
#    allele_list = []
#    
#    for i in range(0,len(allele_pairs)):
#        
#        allele_list = allele_list + (list(allele_pairs[i]))
#        
#    
#    a_count = int(sum(np.in1d(allele_list, np.array('A'))))
#    c_count = int(sum(np.in1d(allele_list, np.array('C'))))
#    g_count = int(sum(np.in1d(allele_list, np.array('G'))))
#    t_count = int(sum(np.in1d(allele_list, np.array('T'))))
#    allele_counts = np.array([a_count, c_count, g_count, t_count])
#    
#    return(allele_counts)
#    
#
#def calc_pairwise_diff(alleles1, alleles2):
#    
#    allele1_arr = np.array(list(alleles1))
#    allele2_arr = np.array(list(alleles2))
#    pairwise_diff = 2 - int(sum(np.in1d(allele1_arr, allele2_arr)))
#    
#    return(pairwise_diff)
#    
#def calc_pi(allele_pairs1, allele_pairs2 = None, within = True):
#    
#    num_pwds = 0 #number of pairwise differences
#    sum_pwd = 0 # sum of pairwise differences
#    
#    ## If doing within sample pi calculation, exclude self pairwise difference
#    ## calculations
#    if within == True:
#        
#        for i in range(0, len(allele_pairs1) - 1):
#            
#            alleles1 = str(allele_pairs1[i])
#            
#            for j in range(i + 1, len(allele_pairs1)):
#                
#                
#                alleles2 = str(allele_pairs1[j])
#                num_pwds += 1
#                sum_pwd += calc_pairwise_diff(alleles1, alleles2)
#    ## Otherwise, since dealing with two different sets of samples, calculate all
#    ## pairwise differences        
#    else:
#        
#        if not isinstance (allele_pairs2, np.ndarray):
#            raise TypeError('Must be numpy array containing alleles')
#            
#        for i in range(0, len(allele_pairs1) ):
#            
#            alleles1 = str(allele_pairs1[i])
#            
#            for j in range(0, len(allele_pairs2)):
#                
#                
#                alleles2 = str(allele_pairs2[j])
#                num_pwds += 1
#                sum_pwd += calc_pairwise_diff(alleles1, alleles2)
#                
#    pi = sum_pwd/num_pwds # pi = average pairwise difference
#    return((pi, num_pwds))
#    
#
#def calc_fst(allele_pairs_list):
#    
#    ## Run all within and non-within pi calculations, and take
#    ## weighted average of pi between and pi within
#    
#    # List of within results (pi value and number of pairwise comparisons in each entry)
#    pi_vals_within = []
#    pairwise_counts_within = []
#    # List of between results (pi value and pairwise comparisons in each entry)
#    pi_vals_between = []
#    pairwise_counts_between = []
#    
#    ## Get the pi within results and their weights
#    for i in range(0,len(allele_pairs_list)):
#        allele_pairs_1 = allele_pairs_list[i]
#        (pi_w_i, pw_ct_w_i) =  calc_pi(allele_pairs_1, within = True)
#        pi_vals_within.append(pi_w_i)
#        pairwise_counts_within.append(pw_ct_w_i)
#        del pi_w_i
#        del pw_ct_w_i
#        
#    
#    ## get the pi between results    
#    for i in range(0, len(allele_pairs_list) - 1):
#        
#        allele_pairs_1 = allele_pairs_list[i]
#        
#        for j in range(i + 1, len(allele_pairs_list)):
#            
#            allele_pairs_2 = allele_pairs_list[j]
#            (pi_btw_i_j, pw_ct_btw_i_j) = calc_pi(allele_pairs_1, allele_pairs_2, within = False)
#            pi_vals_between.append(pi_btw_i_j)
#            pairwise_counts_between.append(pw_ct_btw_i_j)
#            del pi_btw_i_j
#            del pw_ct_btw_i_j
#            
#            
#    ## calculate average pi for within group differences
#    pi_within_avg = np.dot(np.array(pi_vals_within, dtype = 'f'), 
#                           np.array(pairwise_counts_within, dtype = 'f'))/np.sum(pairwise_counts_within, dtype = 'f')
#    ## calculate average pi for between group differences
#    pi_between_avg = np.dot(np.array(pi_vals_between, dtype = 'f'),
#                            np.array(pairwise_counts_between, dtype = 'f'))/np.sum(pairwise_counts_between, dtype = 'f')
#    
#    fst = (pi_between_avg - pi_within_avg)/pi_between_avg
#    
#    return(fst)
#    
#
################################################################################
#
### Part 1 ###
#    
#my_files = os.listdir('./')
#rs_683_allele_sets = []
#rs_910_allele_sets = []
#group_names = []
#
#for file in my_files:
#    
#    if re.search('.hmap', file):
#        print('hapmap file found')
#        print(file)
#        ## population group will have name preceding .hmap
#        group_name = re.search('[A-Z]*(?=.hmap)', file).group(0)
#        ## only take in hapmapfiles from desired groups
#        if group_name in ['ASW', 'GIH', 'LWK', 'MEX', 'MKK', 'TSI']:
#            print('group name')
#            print(group_name)
#            group_names.append(group_name)
#            
#            ## find the column wheree samples are, then retrieve the rows corresponding to each SNV position
#            ## and get all the allele pairs
#            (header, start_ind) = retrieve_header(file)
#            rs_683_allele_sets.append(retrieve_allele_pairs(file, 'rs683', start_ind))
#            rs_910_allele_sets.append(retrieve_allele_pairs(file, 'rs910', start_ind))
#        
#### Part 1a
##        
### calculate proportions of alleles within populations and
### across all populations
#        
#all_gr_array_rs_683 = np.zeros(shape = (len(group_names), 4))
#all_gr_array_rs_910 = np.zeros(shape = (len(group_names), 4))
#
#for i in range(0, len(group_names)):
#    
#    print(group_names[i])
#    print('calculating proportions of alleles in population')
#    rs_683_alleles_i = rs_683_allele_sets[i]
#    rs_910_alleles_i = rs_910_allele_sets[i]
#    rs_683_counts_i = allele_set_2_counts(rs_683_alleles_i)
#    rs_910_counts_i = allele_set_2_counts(rs_910_alleles_i)
#    rs_683_proportions_i = np.divide(rs_683_counts_i, np.sum(rs_683_counts_i), dtype = 'f')
#    rs_910_proportions_i = np.divide(rs_910_counts_i, np.sum(rs_910_counts_i), dtype = 'f')
#    fname1 = group_names[i] + '_rs_683_allele_proportions.tsv'
#    fname2 = group_names[i] + '_rs_910_allele_proprotions.tsv'
#    np.savetxt(fname1, rs_683_proportions_i, delimiter = '\t', header = 'A\tC\tG\tT')
#    np.savetxt(fname2, rs_910_proportions_i, delimiter = '\t', header = 'A\tC\tG\tT')
#    
#    all_gr_array_rs_683[i, :] = rs_683_counts_i[:]
#    all_gr_array_rs_910[i, :] = rs_910_counts_i[:]
#    
### get counts and proportions for each allele across all populations
#all_gr_counts_rs_683 = np.sum(all_gr_array_rs_683, axis = 0)
#all_gr_counts_rs_910 = np.sum(all_gr_array_rs_910, axis = 0)
#all_gr_proportions_rs_683 = np.divide(all_gr_counts_rs_683, np.sum(all_gr_counts_rs_683), dtype = 'f')
#all_gr_proportions_rs_910 = np.divide(all_gr_counts_rs_683, np.sum(all_gr_counts_rs_910), dtype = 'f')
#fname1 = 'all_groups_rs_683_allele_proportions.tsv'
#fname2 = 'all_groups_rs_910_allele_proprotions.tsv'
#np.savetxt(fname1, all_gr_proportions_rs_683, delimiter = '\t', header = 'A\tC\tG\tT')
#np.savetxt(fname2, all_gr_proportions_rs_910, delimiter = '\t', header = 'A\tC\tG\tT')
#
#
#### Part 1b
#### Calculate FST
#
#rs_683_fst = calc_fst(rs_683_allele_sets)
#rs_910_fst = calc_fst(rs_910_allele_sets)
#fst_vals = np.array([rs_683_fst, rs_910_fst])
#fname3 = 'fst_values_rs_683_rs_910'
#np.savetxt(fname3, fst_vals, header = 'rs_683_fst rs_910_fst')



###############################################################################
######################   PART 2 ###############################################
###############################################################################

## 2a: YFGWAS

MEF2_var_predict = np.loadtxt('MEF2_variant_effect_prediction.txt', dtype = str)
shape = np.shape(MEF2_var_predict)
MEF2_var_predict = np.array(MEF2_var_predict[1:shape[0], :], dtype = str)
d = {}
d['variant'] = np.array(MEF2_var_predict[:,0], dtype = str)
d['predictor'] = np.array(MEF2_var_predict[:,1], dtype = float)
MEF2_var_predict = pd.DataFrame(d)

TNG_var_predict = np.loadtxt('TNG1_variant_effect_prediction.txt', dtype = str)
shape = np.shape(TNG_var_predict)
TNG_var_predict = np.array(TNG_var_predict[1:shape[0], :], dtype = str)
d = {}
d['variant'] = np.array(TNG_var_predict[:,0], dtype = str)
d['predictor'] = np.array(TNG_var_predict[:,1], dtype = float)
TNG_var_predict = pd.DataFrame(d)

TOS_var_predict = np.loadtxt('TOS1_variant_effect_prediction.txt', dtype = str)
shape = np.shape(TOS_var_predict)
TOS_var_predict = np.array(TOS_var_predict[1:shape[0], :], dtype = str)
d = {}
d['variant'] = np.array(TOS_var_predict[:,0], dtype = str)
d['predictor'] = np.array(TOS_var_predict[:,1], dtype = float)
TOS_var_predict = pd.DataFrame(d)

## Make the dataframes into a pandas series, each index of series corresponding to a gene
## and mapping to appropriate set of variants/predictions
d = {}
d['TOS1'] = TOS_var_predict
d['TNG1'] = TNG_var_predict
d['MEF2'] = MEF2_var_predict
df_series = pd.Series(d)

genotype_results = np.loadtxt('genotyping_results.txt', dtype = str)
shape = np.shape(genotype_results)
d = {}
d['KepralSyndrome'] = genotype_results[1:,0]
d['BendiiSyndrome'] = genotype_results[1:,1]
d['Xenopolycythemia'] = genotype_results[1:,2]
d['TOS1'] = genotype_results[1:,3]
d['TNG1'] = genotype_results[1:,4]
d['MEF2'] = genotype_results[1:,5]
genotype_results_df = pd.DataFrame(d)

def count_variants(df_series, genotype_results, disease):
    
    ## Input: pandas  dataframe
    
    ## For a given disease, quantify co-occurence of variant and 
    ## control/case status
    ## Overivew: First construct a numpy array to contain variant identifiers
    ## and the number of times the variant appears with case and control outcomes.
    ## Keep a dictionary mapping variant names to row index to avoid unnecessary
    ## looping
    ## Then, go through the genotyping results, and map outcomes to variants,
    ## adding occurences of case or control outcomes to the respective counts
    ## for a given variant
    ## It is expected that all variants found within the genotyping
    ## results be found within the set of variants identified in the
    ## variant effect 
    
    if type(genotype_results) != pd.DataFrame:
        raise TypeError('genotype_results must be pandas dataframe')
        
    if disease not in genotype_results.columns:
        raise ValueError(disease + ' not in genotype results dataframe')
    
    

    
    index_dict = {} # store indices for each variant. 
    ind = 0
    
    ## for each dataframe containing variants for a particular gene
    for gene in df_series.index:
        
        ## set a row with the wild type variant for the gene
        wt_variant_name = gene + '_WT_variant'
        
        if ind == 0:
            var_name_array = np.array([wt_variant_name], dtype = str)
            
        else:
            var_name_array = np.concatenate((var_name_array, np.array([wt_variant_name], dtype = str)), axis = 0)
            
        index_dict[wt_variant_name] = ind
        ind += 1
        current_df = df_series[gene]
        current_vars = current_df['variant'].as_matrix()
        
        ## for each of the variants that differ from the wildtype (i.e. pulled from the prediction file)
        for var in current_vars:
            
            var_name = gene + '_' + var
            var_name_array = np.concatenate((var_name_array, np.array([var_name], dtype = str)), axis = 0)
            index_dict[var_name] = ind
            ind += 1
            
    ## Create dataframe with variant names and columns to hold counts for cases
    ## and controls for each variant
    
    case_count = np.zeros(shape = (len(var_name_array),), dtype = int)
    ctrl_count = np.zeros(shape = (len(var_name_array),), dtype = int)
    
    
    
    
    ## now that we have our empty count arrays we go row by row in the genotype results
    ## df, and match variant status to disease status for our chosen disease
    
    for i in range(0, len(genotype_results[disease])):
        
        disease_status_i = genotype_results[disease][i]
        
        for gene in df_series.index:
            ## variant name for variant i for given gene
            var_gene_i = genotype_results[gene][i]
            var_gene_name_i_list = []
            
            ## make name for variant i for given gene, with name for given gene
            ## appended to name. Used to map to appropriate row in our count
            ## table and add count to either case or control
            
            if var_gene_i == 'p.=':
                # wildtype case
                var_gene_name_i = gene + '_WT_variant'
                var_gene_name_i_list.append(var_gene_name_i)
            else:
                # non-wildtype
                # check if there's 1 variant or >1.
                # make a list that contains all identified variants for this
                # row/gene
                # loop through list and map variants to rows of our count
                # dataframe for variant/case/control status
                rex1 = 'p.[A-Z][a-z]{2}[0-9]*[A-Z][a-z]{2}' # one variant
                rex2 = 'p.\[([A-Z][a-z]{2}[0-9]*[A-Z][a-z]{2}[;]*)*\]' # 2 or more repetitions
                match_rex1 = re.match(rex1, var_gene_i)
                match_rex2 = re.match(rex2, var_gene_i)
                
                
                if match_rex1:
                    # one variant
                    var_gene_name_i = gene + '_' + var_gene_i
                    var_gene_name_i_list.append(var_gene_name_i)
                    
                elif match_rex2:
                    rex3 = '[A-Z][a-z]{2}[0-9]*[A-Z][a-z]{2}'
                    var_gene_i_mult = re.findall(rex3, match_rex2.group(0))
                    ## take individual variant names, append gene name, and append to lsit
                    for var_gene_i_mult_j in var_gene_i_mult:
                        var_gene_name_i_list.append(gene + '_p.' + var_gene_i_mult_j)
                        
                else:
                    warnings.warn(var_gene_i + ' does not match regex1 or regex 2 used to distinguish single and multiple variant cases')
                    print('expect variant entry such as p.=, p.Gly437Tyr, or p.[Gly457Tyr;Ala123Trp]')
                    continue
            ## iterate through list, do mapping to row (variant) and column (disease status)
            ## and add instance to count of variant/status instances
            
            for var_gene_name_i_j in var_gene_name_i_list:
                
                if var_gene_name_i_j not in var_name_array:
                    warnings.warn(var_gene_name_i_j + ' not found in count array rows')
                elif var_gene_name_i not in index_dict:
                   warnings.warn(var_gene_name_i_j + ' not found in the index dictionary (index_dict) for count array')
                else:
                    # add 1 to the entry at ctrl_counts[row of variant] or case_counts[row_of_variant]
                    row_use = index_dict[var_gene_name_i_j]
                    
                    if disease_status_i == 'ctrl':
                        ctrl_count[row_use] += 1
                    elif disease_status_i == 'case':
                        case_count[row_use] += 1
                    else:
                        warnings.warn('At row ' + str(i) + ' gene ' + gene + ' disease status listed as ' + disease_status_i)
                    
    d = {}
    d['variant'] = var_name_array
    d['case'] = case_count
    d['ctrl'] = ctrl_count
    count_df = pd.DataFrame(d)
    return(count_df)
                
    
bendii_syndrome_counts = count_variants(df_series, genotype_results_df, disease = 'BendiiSyndrome')
kepral_syndrome_counts = count_variants(df_series, genotype_results_df, disease = 'KepralSyndrome')
Xenopolycythemia_counts = count_variants(df_series, genotype_results_df, disease = 'Xenopolycythemia')
