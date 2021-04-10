#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 17:24:51 2020

@author: lutra
"""
import subprocess
import sqlite3
import glob


def run(cmd):
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    
    
def mkdir(fold):
    if fold[-1:] == '/':
        fold = fold[:-1]
    loc = '/'.join(fold.split('/')[:-1]) +'/'
    check = glob.glob(loc + '*')
#    print(*check, '', '', sep = '\n')
    if fold not in check:
        run(f'mkdir {fold}')
    # else:
    #     print(fold, 'exists')
    return
    
    
def rev_comp(dna):
    trantab = "".maketrans("ACTGactg","TGACtgac")
    return dna[::-1].translate(trantab)


def dna_to_fasta_format(dna):
    dna = [dna[i*60:(i+1)*60] for i in range(int(len(dna)/60)+1)]
    dna=[d for d in dna if d]
    fasta = '\n'.join(dna)
    return fasta


def read_fasta_file(fasta_file):
    fasta = open(fasta_file).read().strip().split('>')[1:]
    fasta_dict = {}
    for f in fasta:
        f = f.split('\n')
        fasta_dict[f[0]] = ''.join(f[1:])
    return fasta_dict


def table_to_sqlite(db_name, table_name, work_table_as_list_of_tab_sep_raws):
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    
    ex1='DROP TABLE IF EXISTS ' + table_name
    cur.execute(ex1)
    conn.commit()
    
    work = work_table_as_list_of_tab_sep_raws
    work = [w.replace('|','_') for w in work]
    work = [w.replace('"','_') for w in work]
    work = [w.replace("'",'_') for w in work]
    
    work = [w.split('\t') for w in work if w]
    
    check_size = max([max([len(ww) for ww in w]) for w in work])
#    print (check_size)
    field_limit = check_size+1
    
    header = work[0]
    header = [h.replace(' ','_') for h in header]
    header = [h.replace('-','_') for h in header]
#    print (*header,'\n')

    while len(work[-1]) != len(header):
        work[-1].append('')
    
    column = 'id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE'
    insert = []
    for i in range(len(header)):
        check_type = work[1]
        insert.append(header[i])
        try:
            check = int(check_type[i])
            column += (', ' + header[i] + ' INTEGER')
        except:
            column += (', ' + header[i] + ' TEXT')
            if 'sequence' in header[i]:
                column += ('({})'.format(field_limit))
    
    ex2='CREATE TABLE {} ({})'.format(table_name, column)
    cur.execute(ex2)
    conn.commit()
    
#    print (len(work[1:]), end = ', ')
    count = 0
    
    insert = ', '.join(insert)
    for w in work[1:]:
        w = ['"'+ww+'"' for ww in w]        
        w = ', '.join(w)
        task = 'INSERT INTO {} ({}) VALUES ({})'.format(table_name, insert, w)
#        print (task)
        cur.execute(task)
        count +=1
#    print (count)
    conn.commit()
    conn.close()
#    print ('\n\n')
    return


def find_overlaps(feats):
    '''needs a prepared list of lists:
        [[aa[0], int(aa[1]), int(aa[2]), aa[3], float(aa[4]), float(aa[5])] for aa in a]
        name, c1, c2, nap, pident, qcovs'''
    feats.sort(key = lambda x: (x[4], x[5], x[1]*-1), reverse=True)
    # print(*feats, '--\n', sep = '\n')
    
    save = []
    report = []
    for feat in feats:
        name, c_st, c_end, cnap, pident, qcovs = feat
        add = True
        
        for s in save:
            if ('+++red' not in name) or ('+++red' in name and '+++red' in s[0]): #ARGs are replaces only with ARGs!
                s_st, s_end, snap = s[1:4]
                check = max(c_st, s_st) - min(c_end, s_end)
                if check <= 0:
                    overlap = abs(check) / (c_end - c_st)
                    rep_block = f'{name} overlap = {round(overlap,1)} => '
                    
                    if  overlap >= 0.6:
                        add = False
                        rep_block += ('removed!\nrm ' + ' '.join(map(str, feat)) +'\n++ ' + ' '.join(map(str, s)) + '\n\n')
                        report.append(rep_block)
                    else:
                        if add:
                            rep_block += ('saved\n++ ' + ' '.join(map(str, feat)) +'\n++ ' + ' '.join(map(str, s)) + '\n\n')
                            report.append(rep_block)

        if add:
            save.append(feat)
    # print(*save, '* * *\n\n', sep ='\n')
    if len(feats) - len(save) != len([r for r in report if 'removed!' in r]):
        print('REPORT PROBLEM!!!', len(feats), len(save), len(report))
        print(*feats, '--\n', sep = '\n')
        print('reported', report, '---', sep = '\n')
        print(*save, '* * *\n\n', sep ='\n')
        
    
    if report:
        report = ''.join(report).strip() + '\n----------------\n\n'
    else:
        report = ''
    return save, report


def ANI(seq1, seq2):
    if len(seq1) != len(seq2):
        print('ERROR! Wrong length!!!', len(seq1), len(seq2))
        return 'error'
    
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    size = len(seq1)
    ident = [1 if seq1[i] == seq2[i] and seq1[i] in 'ACTG' else 0 for i in range(size)]
    if len(ident) != size:
        print('WEIRD SIZE PROBLEM', size, len(ident))
        
    ident = sum(ident)
    return ident/size
    

names = ['run(cmd) -> None', 'mkdir(fold) -> None', 'rev_comp(dna) -> dna_rev_comp']
names += ['dna_to_fasta_format(dna) -> dna_60_nt_per_row', 'read_fasta_file(fasta_file) -> fasta_dict']
names += ['table_to_sqlite(db_name, table_name, work_table_as_list_of_tab_sep_raws) -> None']
names += ['find_overlaps(spec list of feats) -> upd_feats(spec list), report(txt)']
names += ['ANI(seq1, seq2) -> identi/length(float, <=1)']
print('From PLP_main_functions load', *names, '--\n', sep = '\n')