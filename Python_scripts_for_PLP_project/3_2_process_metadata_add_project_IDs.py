#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 15:14:22 2020

@author: lutra
"""
import sqlite3


PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
    dropbox = '/home/shiri/Dropbox/'
else:
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'
    dropbox = '/home/lutra/Dropbox/'
    
    

table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'
database = main + table + '.sqlite3'
print(database)

conn = sqlite3.connect(database)
cur = conn.cursor()

priority = ['MF356679', 'AF234172', 'JQ965645']
ref_names = {'MF356679': 'D6', 'AF234172': 'P1', 'JQ965645': 'SSU5'}
Inc_types = {'D6': ['D6_putative_replicon'], 'P1': ['IncY_1__K02380', 'p0111_1__AP010962'], 'SSU5': ['IncFIB_pHCM2_1__AL513384', 'IncFIB_pKPHS1_1__CP003223', 'IncFIB_H89-PhagePlasmid_1__HG530657', 'IncFIB_pLF82-PhagePlasmid_1__CU638872']}


data = []
task = 'SELECT id, nucleotide, organism, completeness, slen, PLP_status, enterobacteriaceae, D6_putative_replicon_orf42 FROM ' + table
upd_rep = {}
for row in cur.execute(task):
    gid, nucleotide, host, complete, slen, status, replicon, D6_rep = [str(r) for r in row]
    replicon += D6_rep
    if complete == 'complete' and int(slen) < 2000000 and '_' in status:
        proj_rep = Inc_types[status.split('_')[1]]
        define_replicons = [p for p in proj_rep if p in replicon]
        if not define_replicons:
            replicon_number = len(proj_rep)
            upd_rep[gid] = 'other'
        else:
            if len(define_replicons) > 1:
                print('replicon problem!!!', nucleotide, define_replicons)
            replicon_number = proj_rep.index(define_replicons[0])
            if 'D6' in status:
                print(define_replicons[0])
            upd_rep[gid] = define_replicons[0]
        
        num = -1
        if nucleotide in priority:
            num = 0
            ref = nucleotide
        else:
            ref = [p for p in priority if p in status][0]
            num = priority.index(ref) + 1
        
        data.append([gid, ref, num, replicon_number, ' '.join(host.split()[:2]), int(slen)])

upd = {}
for ref in priority:
    data_ref = [d for d in data if d[1] == ref]
    print(ref, len(data_ref))
    data_ref.sort(key = lambda x: (x[2], x[3], x[4], x[5]))
    print(*data_ref[:4], '...', *data_ref[-4:], '--\n\n', sep = '\n')
    
    if ref == 'JQ965645':
        data_ref = data_ref[:-4] + data_ref[-2:] + data_ref[-4:-2]
    
    count = 1
    for d in data_ref:
        upd[d[0]] = f'{ref_names[ref]}_{count}'
        count += 1

print(len(upd))
for u in upd:
    upd_val = upd[u]
    upd_num = upd_val.split('_')[1]
    task = f'UPDATE {table} SET project_ID="{upd_val}", project_ID_number="{upd_num}", Major_replicon="{upd_rep[u]}" WHERE id="{u}"'
    cur.execute(task)

conn.commit()