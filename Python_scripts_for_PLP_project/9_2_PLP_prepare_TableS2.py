#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 22:49:52 2020

@author: lutra

Phage_like_plasmids_SSU5_P1_D6_12Nov20
CREATE TABLE "Phage_like_plasmids_SSU5_P1_D6_12Nov20" (
	"id"	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	"project_ID"	INTEGER,
	"project_ID_number"	INTEGER,
	"nucleotide"	TEXT,
	"biosample"	TEXT,
	"organism"	TEXT,
	"completeness"	TEXT,
	"genome"	TEXT,
	"slen"	INTEGER,
	"shape"	TEXT,
	"PLP_status"	TEXT,
	"MF356679_D6_ref_cov"	INTEGER,
	"AF234172_P1_ref_cov"	INTEGER,
	"JQ965645_SSU5_ref_cov"	INTEGER,
	"Major_replicon"	TEXT,
	"Major_replicon_variant"	TEXT,
	"Major_replicon_sequence"	TEXT,
	"enterobacteriaceae_N"	INTEGER,
	"enterobacteriaceae"	TEXT,
	"ResFinder_N"	INTEGER,
	"ResFinder"	TEXT,
	"AF234172_P1_ref_CDS_N"	INTEGER,
	"AF234172_P1_ref_CDS"	TEXT,
	"JQ965645_SSU5_ref_CDS_N"	INTEGER,
	"JQ965645_SSU5_ref_CDS"	TEXT,
	"ISel_db_N"	INTEGER,
	"ISel_db"	TEXT,
	"VFDB_setB_nt_N"	INTEGER,
	"VFDB_setB_nt"	TEXT,
	"oriT_all_N"	INTEGER,
	"oriT_all"	TEXT,
	"MF356679_D6_ref_CDS_N"	INTEGER,
	"MF356679_D6_ref_CDS"	TEXT,
	"D6_putative_replicon_orf42_N"	INTEGER,
	"D6_putative_replicon_orf42"	TEXT,
	"t4cp_all_blastx_N"	INTEGER,
	"t4cp_all_blastx"	TEXT,
	"relaxase_all_blastx_N"	INTEGER,
	"relaxase_all_blastx"	TEXT,
	"auxiliary_all_blastx_N"	INTEGER,
	"auxiliary_all_blastx"	TEXT,
	"strain"	TEXT,
	"host"	TEXT,
	"plasmid"	TEXT,
	"country"	TEXT,
	"isolation_source"	TEXT,
	"lat_lon"	TEXT,
	"collection_date"	TEXT,
	"collected_by"	TEXT,
	"host_disease"	TEXT,
	"latitude_and_longitude"	TEXT,
	"geographic_location"	TEXT,
	"Other_features"	TEXT
)
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
    
main = '/data/Current_work/Phage_like_plasmids/PLP_final/scripts_to_Github/'


columns = '''
project_ID INTEGER,
nucleotide TEXT, 
biosample TEXT, 
organism TEXT, 
completeness TEXT, 
genome TEXT, 
slen INTEGER, 
shape TEXT, 
PLP_status TEXT, 
MF356679_D6_ref_cov INTEGER, 
AF234172_P1_ref_cov INTEGER, 
JQ965645_SSU5_ref_cov INTEGER, 
Major_replicon INTEGER, 
enterobacteriaceae_N INTEGER, 
enterobacteriaceae TEXT, 
D6_putative_replicon_orf42_N INTEGER,
D6_putative_replicon_orf42 TEXT,
ResFinder_N INTEGER, 
ResFinder TEXT,
ISel_db_N INTEGER, 
ISel_db TEXT, 
VFDB_setB_nt_N INTEGER, 
VFDB_setB_nt TEXT,
MF356679_D6_ref_CDS_N INTEGER,
MF356679_D6_ref_CDS TEXT,
AF234172_P1_ref_CDS_N INTEGER,
AF234172_P1_ref_CDS TEXT,
JQ965645_SSU5_ref_CDS_N INTEGER,
JQ965645_SSU5_ref_CDS TEXT,
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


replicons = ['D6_putative_replicon_orf42_c45699_46596', 'IncY_1__K02380', 'p0111_1__AP010962', 'IncFIB_pHCM2_1__AL513384', 'IncFIB_pKPHS1_1__CP003223', 'IncFIB_H89-PhagePlasmid_1__HG530657', 'IncFIB_pLF82-PhagePlasmid_1__CU638872']
ref_cov = {'D6': 'MF356679_D6_ref_cov', 'IncY': 'AF234172_P1_ref_cov', 'p0111': 'AF234172_P1_ref_cov', 'IncFIB': 'JQ965645_SSU5_ref_cov'}
ref_cds_n = {'D6': 'MF356679_D6_ref_CDS_N', 'IncY': 'AF234172_P1_ref_CDS_N', 'p0111': 'AF234172_P1_ref_CDS_N', 'IncFIB': 'JQ965645_SSU5_ref_CDS_N'}

table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'
database = main + table + '.sqlite3'
print(database)

conn = sqlite3.connect(database)
cur = conn.cursor()

data = []

for rep in replicons[:]:
    print(rep)
    if rep[:2] == 'D6':
        rep_column = 'D6_putative_replicon_orf42'
    else:
        rep_column = 'enterobacteriaceae'
    
    rep_data = []
    task = 'SELECT ' + ', '.join(task_columns) + ' FROM ' + table
    task += f" WHERE {rep_column} LIKE '%{rep}%'"
    for row in cur.execute(task):
        row = [str(r) for r in row]
        
        rep_match = row[task_columns.index(rep_column)].split('\n')
        rep_match = [r for r in rep_match if rep in r]
        rep_match = '\n'.join(rep_match)
        # print(rep_match)
        
        ref = rep.split('_')[0]
        complete =row[task_columns.index('completeness')]
        slen = row[task_columns.index('slen')]
        cov = row[task_columns.index(ref_cov[ref])]
        cds_n = row[task_columns.index(ref_cds_n[ref])]
        
        cov = float(cov)
        
        reason = ''
        if complete == 'complete' and int(slen) < 2000000 and cov >= 40:
            reason = 'grouped'
        elif complete != 'complete':
            reason = 'not complete'
        elif int(slen) >= 2000000:
            reason = 'chromosome'
        elif cov < 40:
            reason = f'lower_coverage qcovs:{cov}% CDSs:{cds_n}'
            
        if not reason:
            print('REASON PROBLEM', row[1])
        
        shorten_annot = row[task_columns.index('Major_replicon')+1:task_columns.index('strain')]
        shorten = []
        for s in shorten_annot:
            if row.index(s) == task_columns.index('ResFinder'):
                s = [ss.split()[1] for ss in s.split('\n')]
            else:
                s = [ss.split()[0] for ss in s.split('\n')]
            shorten.append('\n'.join(s))
        row = row[:task_columns.index('Major_replicon')+1] + shorten + row[task_columns.index('strain'):]
        
        new_row = [rep, rep_match, reason] + row
        rep_data.append(new_row)
        
    rep_data.sort(key = lambda x: (x[2], x[3].split('_')[0], int(x[3].split('_')[-1])))
    possible_reasons = ['grouped', 'not complete', 'chromosome', 'lower_coverage']
    for pr in possible_reasons:
        pr_filt = [r for r in rep_data if pr in r[2]]
        print(pr, len(pr_filt))
    
    print(len(rep_data))
    data += rep_data
    
    print('--\n\n\n')


print(len(data))

name = 'STable2.Distribution_of_PLP_related_replicons_26Aug20'
table2_csv = main + name + '.csv'
with(open(table2_csv, 'w')) as fh:
    header = '\t'.join(['Replicon', 'Replicon_match', 'Grouped_or_why_not_grouped'] + task_columns)
    fh.write(header + '\n')
    
    data = ['\t'.join(d) for d in data]
    data = [d.replace('\n', '<<AND>>').replace('**', '#').replace('--', '_').replace('__', '_').replace('__', '_') for d in data]
    fh.write('\n'.join(data) + '\n')
    
# csv_to_sqlite3(table2_csv, main + name + '.sqlite3')


# Reading the csv file 
df_new = pd.read_csv(table2_csv, sep='\t')
df_new = df_new.replace('<<AND>>', ' ', regex=True)
  
# saving xlsx file 
GFG = pd.ExcelWriter(main + name + '.xlsx') 
df_new.to_excel(GFG, index = False) 
  
GFG.save() 
os.remove(table2_csv)