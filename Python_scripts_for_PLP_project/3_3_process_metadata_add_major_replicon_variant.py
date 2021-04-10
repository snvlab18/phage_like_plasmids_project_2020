#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 04:03:51 2020

@author: lutra
"""
import sqlite3

def dna_to_fasta_format(dna):
    dna = [dna[i*60:(i+1)*60] for i in range(int(len(dna)/60)+1)]
    dna=[d for d in dna if d]
    fasta = '\n'.join(dna)
    return fasta



PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
    dropbox = '/home/shiri/Dropbox/'
else:
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'
    dropbox = '/home/lutra/Dropbox/'
    

annotate_fold = main + 'annotate/'    

table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'
database = main + table + '.sqlite3'
print(database)

conn = sqlite3.connect(database)
cur = conn.cursor()

data = []
task = 'SELECT id, nucleotide, Major_replicon FROM ' + table
upd = {}
collect_vars = {}
for row in cur.execute(task):
    gid, nucleotide, major = [str(r) for r in row]
    if major != '0' and major != 'other':                   #currently exclude others1
        annotate = annotate_fold + nucleotide + '.txt'
        annot = open(annotate).read().strip().split('\n')
        rep_seq = [a for a in annot if major in a.replace(')','').replace('(','_')]
        if len(rep_seq) == 2 and nucleotide != 'MK356558':
            rep_seq = [r for r in rep_seq if 'decircle' in r]
        if len(rep_seq) != 1 and nucleotide != 'MK356558':
            print(nucleotide)
            print(rep_seq)
            
        rep_seq = [r.split('\t')[7] for r in rep_seq]
        upd[gid] = rep_seq
        
        if major not in collect_vars:
            collect_vars[major] = {}
        
        for seq in rep_seq:
            if seq not in collect_vars[major]:
                collect_vars[major][seq] = 0
            collect_vars[major][seq] += 1

collect_names = {}
for rep in collect_vars:
    order_and_name = []
    print(rep, len(collect_vars[rep]))
    for seq in collect_vars[rep]:
        order_and_name.append(seq)
    order_and_name.sort(key = lambda x: collect_vars[rep][x], reverse = True)
    for i in range(len(order_and_name)):
        seq = order_and_name[i]
        seq_n = collect_vars[rep][seq]
        collect_names[seq] = f'{rep}_variant{i+1}_n{seq_n}'
    print('--\n\n')
    
for u in upd:
    upd_seq = []
    upd_name = []
    for seq in upd[u]:
        seq_name = collect_names[seq]
        # seq_fasta = dna_to_fasta_format(seq)
        # fasta = f'>{seq_name}\n{seq_fasta}\n'
        # upd_seq.append(fasta)
        upd_seq.append(seq)
        upd_name.append(seq_name)
    upd_seq = '\n'.join(upd_seq)
    upd_name = '\n'.join(upd_name)
    task = f'UPDATE {table} SET Major_replicon_variant="{upd_name}", Major_replicon_sequence="{upd_seq}" WHERE id="{u}"'
    cur.execute(task)

conn.commit()