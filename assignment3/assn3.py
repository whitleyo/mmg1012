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
import random as rand
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

## Functions for Part 2

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


def fisher_p_table(count_df):
    ## get variants with at least one identification
    ## run fisher's exact test
    ## return dataframe with fisher p values and the contingency table counts
    ## in each row
    
    sorted_df = count_df.sort_values(['case', 'ctrl'], ascending = False)
    variants = sorted_df['variant']
    counts_only_df = sorted_df[['case', 'ctrl']]
    counts_only_array = counts_only_df.as_matrix()
    variants_keep = 0
    
    for i in range(0,sorted_df.shape[0]):
        
        if (np.sum(counts_only_array[i,0:2] > 0)):
            variants_keep += 1
        else:
            print('breaking, found ' + str(i - 1) + ' variants with at least one instance')
            break
        
    if variants_keep > 0:
        kept_variants = variants[0:variants_keep]
    else:
        print('did not find any variants with > 1 instance')
        return(None)
    
    counts_keep = counts_only_array[0:variants_keep,:]
    
    not_var_case = []
    not_var_ctrl = []
    p_vals = []
    ## run fisher's exact test for every row
    for i in range(0, variants_keep):
        
        # get counts for contingency table
        case_count_var_i = counts_keep[i, 0]
        ctrl_count_var_i = counts_keep[i, 1]
        case_count_not_var_i = np.sum(counts_keep[:,1], dtype = 'int64') - case_count_var_i
        ctrl_count_not_var_i = np.sum(counts_keep[:,1], dtype = 'int64') - ctrl_count_var_i
        
        cont_table = np.array([[case_count_var_i, ctrl_count_var_i], 
                               [case_count_not_var_i, ctrl_count_not_var_i]])
        ## get p value for 2 sided fisher's exact test
        odds_ratio, fisher_p_val = stats.fisher_exact(cont_table)
        not_var_case.append(case_count_not_var_i)
        not_var_ctrl.append(ctrl_count_not_var_i)
        p_vals.append(fisher_p_val)
    
    # make dataframe
    new_df ={}
    new_df['case_var'] = counts_keep[0:variants_keep,0]
    new_df['ctrl_var'] = counts_keep[0:variants_keep,1]
    new_df['case_not_var'] = not_var_case
    new_df['ctrl_not_var'] = not_var_ctrl
    new_df['p_val'] = p_vals
    new_df['variant'] = kept_variants
    new_df = pd.DataFrame(new_df)
    return(new_df)
    
    
def BH_filter(fisher_p_df, FDR_cutoff = 0.01):
    ## Filter using Benjamini Hochberg FDR control
    fisher_p_sorted = fisher_p_df.sort_values('p_val', ascending = True)
    p_vals = fisher_p_sorted['p_val'].as_matrix()
    
    num_kept_variants = 0
    for i in range(0, len(p_vals)):
        q_val = p_vals[i]*float(i)/float(len(p_vals))
        
        if q_val >= FDR_cutoff:
            break
        
        else:
            num_kept_variants += 1
    
    if num_kept_variants == 0:
        return(None)
    else:
        fisher_p_sorted = fisher_p_sorted[0:num_kept_variants]
        return(fisher_p_sorted)
    
    
    
#def sep_by_gene(count_df):
#    
#    ## Variant burden: number of 
#    
#    sorted_df = count_df.sort_values(['variant'], ascending = False)
#    gene_names = [] # will hold names of genes found
#    gene_indices = [] # will hold indices of dataframe corresponding to gene
#    
#    
#    last_gene_name = '' # reset within loop. used to check if current gene name matches last one
#    current_index_list = [] # list of indices for current gene. Reset if new gene encountered
#    rex = '[A-Z]{1,3}[0-9]'
#
#    for i in range(0, sorted_df.shape[0]):
#        gene_name_i = re.match(rex, sorted_df['variant'][i])
#        
#        try:
#            if last_gene_name != gene_name_i.group(0):
#                
#                gene_names.append(gene_name_i.group(0))
#                last_gene_name = gene_name_i.group(0)
#                if i > 0:
#                    current_index_list.append(i-1)
#                    gene_indices.append(current_index_list)
#                else:
#                    current_index_list.append(i)
#                    gene_indices.append(current_index_list)
#                
#            elif i == (sorted_df.shape[0] - 1):
#                current_index_list.append(i)
#                gene_indices.append(current_index_list)
#                
#                
#        except:
#            if not gene_name_i:
#                raise ValueError(sorted_df['variant'][i] + ' not a gene name')
#            
#    ## make pandas panel containing dataframe for each gene
#    d = {}
#    for i in range(0, len(gene_names)):
#        gene = gene_names[i]
#        indices = gene_indices[i]
#        start = indices[0]
#        end = indices[1] + 1
#        df_subs = sorted_df[start:end]
#        d[gene] = df_subs
#        
#    gene_disease_dfs = pd.Series(d)
#    return(gene_disease_dfs)
    
##
    
def calculate_variant_burden(genotype_results_df, # table of case/control status for multiple disease and genes
                             disease, # string, should be disease name
                             gene, # gene name
                             resample = False, # if set to true, will randomly resample 
                             resample_pct = 0.80,
                             seed = 2,
                             stat = 'case',
                             mode = 'Normal', 
                             polyphen2_scores = None # pandas dataframe containing variants and predicted scores
                             ):
    # calculate variant burden for a given gene for a given disease status for a given disease
    # calculated as number of variants in poulation with disease status divided by
    # number of individuals of population with disease status
    
    ## provide a list of indices if subsampling desired
    if disease not in genotype_results_df.columns:
        raise ValueError('disease name must be in columns of genotype_results_df')
        
    disease_status = genotype_results_df[disease].as_matrix()
    
    if not (np.sum(np.in1d(disease_status, ['case', 'ctrl'])) == len(disease_status)):
        raise ValueError('check disease name, as not all elements of column are case or ctrl')
        
    if stat not in ['case', 'ctrl']:
        raise ValueError('stat must be case or control')
    
    num_subpop = np.sum(np.in1d(disease_status, stat), dtype = 'int64') # number of individuals in subpop (case or control)
    # get variants corresponding to our set of case/control population
    variant = genotype_results_df[gene].as_matrix()
    variant_sub = variant[np.in1d(disease_status, stat)]
    
    # resample variants within our subpop
    if resample:
        
        num_resample = int(num_subpop*resample_pct)
        rand.seed(seed)
        indices_use = resample(range(0,num_subpop, num_resample))
        sample_size = num_resample
        variant_sub = variant_sub[indices_use]
    
    else:
        
        sample_size = num_subpop
        
        
    # this score will be incremented each time we encounter a variant
    # within our subpopulation
    variant_score = 0
    
    # regexes to find variants
    rex1 = 'p.[A-Z][a-z]{2}[0-9]*[A-Z][a-z]{2}' # one variant
    rex2 = 'p.\[([A-Z][a-z]{2}[0-9]*[A-Z][a-z]{2}[;]*)*\]' # 2 or more variants
    rex3 = '[A-Z][a-z]{2}[0-9]*[A-Z][a-z]{2}' # regex for all variants if 2 or more variants     
    
    if mode == 'PolyPhen2':
    
        ppn2_variant = polyphen2_scores['variant']
        ppn2_predictor = polyphen2_scores['predictor']
        
        
    for var in variant_sub:
        # find out whether wiltype, variant, or multiple variants, add appropriate
        # amount to variant_score
        match_rex1 = re.match(rex1, var) # check for there being one variant
        match_rex2 = re.match(rex2, var) # check for there being 2 or more variants

        if match_rex1:
            
            match_list = [match_rex1.group(0)]
        
        elif match_rex2:
            
            match_list = re.findall(rex3 ,match_rex2.group(0))
            
        elif var == 'p.=':
            continue
            
        else:
            
            raise ValueError(var + ' does not match up to variant finding regexes')
            
        
        if mode == 'Normal':
            
            variant_score += len(match_list)
            
        elif mode == 'PolyPhen2':
            
            for var_i in match_list:
                
                var_ind = np.in1d(ppn2_variant, var_i)
                polyphen2_score = ppn2_predictor[var_ind]
                variant_score += polyphen2_score
                
        else:
            
            raise ValueError('mode must be Normal or PolyPhen2')
                
    norm_score = float(variant_score)/float(sample_size)
    
    return(norm_score)
    
def score_and_p_val( genotype_results_df, # table of case/control status for multiple disease and genes
                     disease, # string, should be disease name
                     gene, # gene namee 
                     resample_pct = 0.80,
                     base_seed = 2,
                     mode = 'Normal', 
                     polyphen2_scores = None ,
                     iterations = 100):
    
    if mode == 'PolyPhen2':        
        
        if type(polyphen2_scores) != pd.DataFrame:
            
            raise TypeError('polyphen2_scores must be a pandas dataframe')
            
    elif not mode == 'Normal':
        
        raise ValueError('mode argument must be Normal or PolyPhen2')
        
    # First calculate variant score for case condition
    var_score_case = calculate_variant_burden(genotype_results_df, # table of case/control status for multiple disease and genes
                                             disease, # string, should be disease name
                                             gene, # gene name
                                             resample = False, # if set to true, will randomly resample 
                                             resample_pct = 0.80,
                                             seed = 2,
                                             stat = 'case',
                                             mode = mode, 
                                             polyphen2_scores = None # pandas dataframe containing variants and predicted scores
                                             )
    ## calcluate null distribution by resampling control group for n iteratios
    var_scores_null = []
    
    for i in range(0, iterations):
        
        current_seed = base_seed + i
        var_score_null_i = calculate_variant_burden(genotype_results_df, # table of case/control status for multiple disease and genes
                                                     disease, # string, should be disease name
                                                     gene, # gene name
                                                     resample = False, # if set to true, will randomly resample 
                                                     resample_pct = 0.80,
                                                     seed = current_seed,
                                                     stat = 'ctrl',
                                                     mode = mode, 
                                                     polyphen2_scores = polyphen2_scores # pandas dataframe containing variants and predicted scores
                                                                 )
        var_scores_null.append(var_score_null_i)
        
    ## calculate one sided p value
    var_scores_null_arr = np.array(var_scores_null)
    
    num_gt = np.sum(var_scores_null_arr >= var_score_case)
    num_lt = np.sum(var_scores_null_arr > var_score_case)
    
    if num_gt <= num_lt:
        p_val = float(num_gt)/len(var_scores_null_arr)
    else:
        p_val = float(num_lt)/len(var_scores_null_arr)
        
    return(var_score_case, p_val, var_scores_null)
    

############################################
######### 2A: YFGWAS #######################

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

bendii_syndrome_counts = count_variants(df_series, genotype_results_df, disease = 'BendiiSyndrome')
kepral_syndrome_counts = count_variants(df_series, genotype_results_df, disease = 'KepralSyndrome')
Xenopolycythemia_counts = count_variants(df_series, genotype_results_df, disease = 'Xenopolycythemia')

bendii_syndrome_fisher = fisher_p_table(bendii_syndrome_counts)
kepral_syndrome_fisher = fisher_p_table(kepral_syndrome_counts)
Xenopolycythemia_fisher = fisher_p_table(Xenopolycythemia_counts)


########################################################################################
###### 2B: Scoring by mutational burden #####################################
#######################################################################


BS_burden_MEF1_score, BS_burden_MEF1_p_val, null  = score_and_p_val(genotype_results_df, 'BendiiSyndrome', 'MEF2')
BS_burden_TNG1_score, BS_burden_TNG1_p_val, null  = score_and_p_val(genotype_results_df, 'BendiiSyndrome', 'TNG1')
BS_burden_TOS1_score, BS_burden_TOS1_p_val, null  = score_and_p_val(genotype_results_df, 'BendiiSyndrome', 'TOS1')

KS_burden_MEF1_score, KS_burden_MEF1_p_val, null  = score_and_p_val(genotype_results_df, 'KepralSyndrome', 'MEF2')
KS_burden_TNG1_score, KS_burden_TNG1_p_val, null  = score_and_p_val(genotype_results_df, 'KepralSyndrome', 'TNG1')
KS_burden_TOS1_score, KS_burden_TOS1_p_val, null  = score_and_p_val(genotype_results_df, 'KepralSyndrome', 'TOS1')

X_burden_MEF1_score, XS_burden_MEF1_p_val, null  = score_and_p_val(genotype_results_df, 'Xenopolycythemia', 'MEF2')
X_burden_TNG1_score, XS_burden_TNG1_p_val, null  = score_and_p_val(genotype_results_df, 'Xenopolycythemia', 'TNG1')
X_burden_TOS1_score, XS_burden_TOS1_p_val, null  = score_and_p_val(genotype_results_df, 'Xenopolycythemia', 'TOS1')

#####################################################################################
########## 2C: Scoring by polyphen2 prediction (aggregate) ##########################
#####################################################################################

BS_burden_MEF1_score_ppn2, BS_burden_MEF1_p_val_ppn2, null  = score_and_p_val(genotype_results_df, 'BendiiSyndrome', 'MEF2', polyphen2_scores = MEF2_var_predict)
BS_burden_TNG1_score_ppn2, BS_burden_TNG1_p_val_ppn2, null  = score_and_p_val(genotype_results_df, 'BendiiSyndrome', 'TNG1', polyphen2_scores = TNG_var_predict)
BS_burden_TOS1_score_ppn2, BS_burden_TOS1_p_val_ppn2, null  = score_and_p_val(genotype_results_df, 'BendiiSyndrome', 'TOS1', polyphen2_scores = TOS_var_predict)

KS_burden_MEF1_score_ppn2, KS_burden_MEF1_p_val_ppn2, null  = score_and_p_val(genotype_results_df, 'KepralSyndrome', 'MEF2', polyphen2_scores = MEF2_var_predict)
KS_burden_TNG1_score_ppn2, KS_burden_TNG1_p_val_ppn2, null  = score_and_p_val(genotype_results_df, 'KepralSyndrome', 'TNG1', polyphen2_scores = TNG_var_predict)
KS_burden_TOS1_score_ppn2, KS_burden_TOS1_p_val_ppn2, null  = score_and_p_val(genotype_results_df, 'KepralSyndrome', 'TOS1', polyphen2_scores = TOS_var_predict)

X_burden_MEF1_score_ppn2, XS_burden_MEF1_p_val_ppn2, null  = score_and_p_val(genotype_results_df, 'Xenopolycythemia', 'MEF2', polyphen2_scores = MEF2_var_predict)
X_burden_TNG1_score_ppn2, XS_burden_TNG1_p_val_ppn2, null  = score_and_p_val(genotype_results_df, 'Xenopolycythemia', 'TNG1', polyphen2_scores = TNG_var_predict)
X_burden_TOS1_score_ppn2, XS_burden_TOS1_p_val_ppn2, null  = score_and_p_val(genotype_results_df, 'Xenopolycythemia', 'TOS1', polyphen2_scores = TOS_var_predict)
