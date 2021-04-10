#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 10:50:32 2020

@author: lutra
"""
import sqlite3
import altair as alt
import pandas as pd
from PLP_main_functions import find_overlaps


PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
else:   
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'


print(main)


table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'

database = main + table + '.sqlite3'
conn = sqlite3.connect(database)
cur = conn.cursor()

arg_db = []
overlap_reports = []
count = 0

task = f'SELECT nucleotide, completeness, slen, PLP_status, ResFinder FROM {table}'
for row in cur.execute(task):
    gid, complete, slen, status, resfinder = [str(r) for r in row]
    if complete == 'complete' and int(slen) < 2000000 and len(status) > 4 and resfinder != '0':
        count +=1
        resfinder = resfinder.strip().split('\n')
        resf = [aa.split()[1:] for aa in resfinder]
        resf = [[aa[0], int(aa[1]), int(aa[2]), aa[3], float(aa[4]), float(aa[5])] for aa in resf]
        
        resf_classes = {}
        for r in resfinder:
            arg_cl, arg = r.split()[:2]
            arg = arg.replace('**', '').replace('#', '')
            if arg not in resf_classes:
                resf_classes[arg] = []
            resf_classes[arg].append(arg_cl)
            
        
        upd_resf, report = find_overlaps(resf)
        for r in upd_resf:
            arg = r[0].replace('**', '').replace('#', '')
            arg_cl = resf_classes[arg]
            arg_cl.sort()
            arg_cl = ','.join(arg_cl)
            arg_db.append([status, arg_cl, arg])
            
        if report:
            report = f'{status} {gid}\n{report}'
            overlap_reports.append(report)


print('ARG_db', count, len(arg_db))

args = [a[2] for a in arg_db]
arg_list = list(set(args))
arg_list.sort(key = lambda x: args.count(x), reverse = True)

print(len(arg_list))
for a in arg_list[:10]:
    print(a, args.count(a))

arg_list.sort()
check = 0
with open(main + 'ARG_maps/ARG_lists_upd.tsv', 'w') as fh:
    header = ['ARG class', 'ARG', 'Total N', 'D6', 'P1', 'SSU5']
    fh.write('\t'.join(header) + '\n')
    
    collect = []
    for a_gene in arg_list:
        a_extract = [a for a in arg_db if a[2] == a_gene]
        a_class = a_extract[0][1]
        total = len(a_extract)
        check += total
        entry = [a_class, a_gene, str(total)]
        check2 = 0
        for phage in ('D6', 'P1', 'SSU5',):
            phage_count = len([a for a in arg_db if a[2] == a_gene and phage in a[0]])
            check2 += phage_count
            entry.append(str(phage_count))
            
        if check2 != total:
            print('local count problem', a_gene, total, check2)
        collect.append(entry)
    collect.sort(key = lambda x: (x[0], int(x[2])*-1, x[1]))
    collect = ['\t'.join(c) for c in collect]
    fh.write('\n'.join(collect) + '\n')


if overlap_reports:
    with open(main + 'ARG_maps/ARG_lists_overlap_report.txt', 'w') as fh:
        fh.write('\n'.join(overlap_reports))


if check != len(arg_db):
    print('GRAND COUNT PROBLEM', len(arg_db), check)
        
    

if 1:
    d = {'PLP_group': [a[0] for a in arg_db], 'ARG_class': [a[1] for a in arg_db], 'ARGs': [a[2] for a in arg_db]}
    df = pd.DataFrame(data=d)
      
    chart = alt.Chart(df).mark_bar().encode(
        column='PLP_group',
        x='ARG_class',
        y='count(ARG_class):Q',
        color='ARG_class').properties(width=220)
    
    chart.show()
