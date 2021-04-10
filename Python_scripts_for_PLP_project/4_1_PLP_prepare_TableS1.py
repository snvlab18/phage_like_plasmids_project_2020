#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 22:49:52 2020

@author: lutra

Phage_like_plasmids_SSU5_P1_D6_12Nov20

project_ID INTEGER, 
project_ID_number INTEGER, 
nucleotide TEXT, 
biosample TEXT, 
organism TEXT, 
completeness TEXT, 
genome TEXT, 
slen INTEGER, 
shape INTEGER, 
PLP_status TEXT, 
MF356679_D6_ref_cov INTEGER, 
AF234172_P1_ref_cov INTEGER, 
JQ965645_SSU5_ref_cov INTEGER, 
Major_replicon INTEGER, 
enterobacteriaceae_N INTEGER, 
enterobacteriaceae TEXT, 
ResFinder_N INTEGER, 
ResFinder TEXT, 
AF234172_P1_ref_CDS_N INTEGER, 
AF234172_P1_ref_CDS TEXT, 
JQ965645_SSU5_ref_CDS_N INTEGER, 
JQ965645_SSU5_ref_CDS INTEGER, 
ISel_db_N INTEGER, 
ISel_db INTEGER, 
VFDB_setB_nt_N INTEGER, 
VFDB_setB_nt INTEGER, 
oriT_all_N INTEGER, 
oriT_all INTEGER, 
MF356679_D6_ref_CDS_N INTEGER, 
MF356679_D6_ref_CDS INTEGER, 
D6_putative_replicon_orf42_N INTEGER, 
D6_putative_replicon_orf42 INTEGER, 
t4cp_all_blastx_N INTEGER, 
t4cp_all_blastx INTEGER, 
relaxase_all_blastx_N INTEGER, 
relaxase_all_blastx INTEGER, 
auxiliary_all_blastx_N INTEGER, 
auxiliary_all_blastx INTEGER, 
strain TEXT, 
host TEXT, 
plasmid TEXT, 
country TEXT, 
isolation_source TEXT, 
lat_lon TEXT, 
collection_date TEXT, 
collected_by TEXT, 
host_disease TEXT, 
latitude_and_longitude TEXT, 
geographic_location TEXT, 
Other_features TEXT
"""
from csv_in_tables_to_sqlite3 import csv_to_sqlite3
import sqlite3
import os

import pandas as pd 
import numpy as np 
  
  

PC_lab = True

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
else:
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'


columns = '''
project_ID INTEGER,
nucleotide TEXT, 
biosample TEXT, 
organism TEXT, 
completeness TEXT, 
genome TEXT, 
slen INTEGER, 
shape INTEGER, 
PLP_status TEXT, 
MF356679_D6_ref_cov INTEGER, 
AF234172_P1_ref_cov INTEGER, 
JQ965645_SSU5_ref_cov INTEGER, 
Major_replicon INTEGER, 
enterobacteriaceae_N INTEGER, 
enterobacteriaceae TEXT, 
ResFinder_N INTEGER, 
ResFinder TEXT,
ISel_db_N INTEGER, 
ISel_db INTEGER, 
VFDB_setB_nt_N INTEGER, 
VFDB_setB_nt INTEGER,
strain TEXT, 
host TEXT, 
plasmid TEXT, 
country TEXT, 
isolation_source TEXT, 
lat_lon TEXT, 
collection_date TEXT, 
collected_by TEXT, 
host_disease TEXT, 
latitude_and_longitude TEXT, 
geographic_location TEXT
'''
columns = columns.strip().split('\n')
columns = [c.replace(',', '').split() for c in columns]

task_columns = [c[0] for c in columns]
print(*task_columns, len(task_columns), '--\n\n', sep = '\n')


table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'
database = main + table + '.sqlite3'
print(database)

conn = sqlite3.connect(database)
cur = conn.cursor()

data = []
task = 'SELECT ' + ', '.join(task_columns) + ' FROM ' + table
for row in cur.execute(task):
    row = [str(r) for r in row]
    if row[0] != '0':
        shorten_annot = row[task_columns.index('Major_replicon')+1:task_columns.index('strain')]
        shorten = []
        for s in shorten_annot:
            if row.index(s) == task_columns.index('ResFinder'):
                s = [ss.split()[1] for ss in s.split('\n')]
            else:
                s = [ss.split()[0] for ss in s.split('\n')]
            shorten.append('\n'.join(s))
        row = row[:task_columns.index('Major_replicon')+1] + shorten + row[task_columns.index('strain'):]
        data.append(row)

print(len(data))
for phage in ['D6', 'P1', 'SSU5']:
    print(phage, len([d for d in data if phage in d[0]]))
    

data.sort(key = lambda x: (x[0].split('_')[0], int(x[0].split('_')[-1])))
print([d[0] for d in data][:10], [d[0] for d in data][45:55], [d[0] for d in data][-10:], '--\n\n', sep = '\n')



name = 'STable1.Phage_like_plasmids_SSU5_P1_D6_26Aug20'
table1_csv = main + name + '.csv'
with(open(table1_csv, 'w')) as fh:
    header = '\t'.join(task_columns)
    fh.write(header + '\n')
    
    data = ['\t'.join(d) for d in data]
    data = [d.replace('\n', '<<AND>>').replace('**', '#').replace('--', '_').replace('__', '_').replace('__', '_') for d in data]
    fh.write('\n'.join(data) + '\n')
    
csv_to_sqlite3(table1_csv, main + name + '.sqlite3')


# Reading the csv file 
df_new = pd.read_csv(table1_csv, sep='\t')
df_new = df_new.replace('<<AND>>', ' ', regex=True)
  
# saving xlsx file 
GFG = pd.ExcelWriter(main + name + '.xlsx') 
df_new.to_excel(GFG, index = False) 
  
GFG.save() 
os.remove(table1_csv)