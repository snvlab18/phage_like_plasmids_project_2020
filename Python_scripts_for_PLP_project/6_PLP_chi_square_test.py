#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:22:27 2020

@author: lutra
"""
import  sqlite3
from scipy.stats import chi2_contingency

def chi_square_test(table):
    '''Perform a chi-square test with Yates correction on the observed data'''
    #example
    #table = [[90, 605], [178, 1588]]
    #chi, p = chi_square_test(table)
    
    chi_square, p_value, dof, f_expected = chi2_contingency(table, correction=True)
    print(f_expected, '', sep = '\n')
    print (f'Scipy solution:\nChi-squared Statistic: {chi_square}\np-value: {p_value}\nDegrees of Freedom: {dof}')
    return chi_square, p_value


PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
else:
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'

table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'
database = main + table + '.sqlite3'

conn = sqlite3.connect(database)
cur = conn.cursor()


#ARGs are accosiated with additional replicons
if 1:
    task = 'SELECT organism, project_ID, Major_replicon, enterobacteriaceae, ResFinder_N FROM ' + table
    data = []
    for row in cur.execute(task):
        host, proj_ID, major, reps, resf = [str(r) for r in row]
        if proj_ID != '0':
            if resf == '0':
                resf = '_noARG'
            else:
                resf = '_ARGpos'
            
            minor = 'no_minor_reps'
            if reps != '0':
                reps = reps.split('\n')
                reps = [r for r in reps if major not in r]
                if reps:
                    minor = 'minor_reps'
            
            data.append(minor + resf)
                
        
    print('N =', len(data), set(data))

    var = ['minor_reps_ARGpos', 'minor_reps_noARG', 'no_minor_reps_ARGpos', 'no_minor_reps_noARG']
    print(*var, '', sep = '\n')
    
    chi_table = [[data.count(var[0]), data.count(var[1])], [data.count(var[2]), data.count(var[3])]]
    print(*chi_table, '', sep = '\n')

    chi, p = chi_square_test(chi_table)
    print(f'\nchi = {chi}, p = {p}')


#Y_pestis_is_associated_with_more_vir_genes_than_all_other_species_in_all_three_PLP_groups
if 0:
    task = 'SELECT organism, project_ID, VFDB_setB_nt FROM ' + table
    data = []
    for row in cur.execute(task):
        host, proj_ID, vfdb = [str(r) for r in row]
        if proj_ID != '0':
            virulence = '_2no'
            if vfdb != '0':
                virulence = '_1yes'
            if 'Yersinia pestis' in host:
                data.append('Yersinia' + virulence)
            else:
                data.append('Z_other' + virulence)
        
    print('N =', len(data))
    var = ['Yersinia_1yes', 'Yersinia_2no', 'Z_other_1yes', 'Z_other_2no']
    print(*var, '', sep = '\n')    
    
    chi_table = [[data.count(var[0]), data.count(var[1])], [data.count(var[2]), data.count(var[3])]]
    print(*chi_table, '', sep = '\n')

    chi, p = chi_square_test(chi_table)
    print(f'\nchi = {chi}, p = {p}')