#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 00:42:18 2020

@author: lutra
"""
import glob


main = '/data/Current_work/Phage_like_plasmids/PLP_final/GenBank/'

cd_clstr = main + 'PLP_genbanks.clstr'
all_prots = main + 'PLP_genbanks.faa'

all_prots = [a[1:].split(' ') for a in open(all_prots).read().split('\n') if a[:1] == '>']
prot_names = [a[0] for a in all_prots]
prot_descr = [' '.join(a[1:]) for a in all_prots]

all_prots = dict(zip(prot_names, prot_descr))
print('all_prots', len(all_prots))


prot_db = {}
gb_files = [g for g in glob.glob(main + '*') if '.gb' in g]
for gb in gb_files:
    gb_name = gb.split('/')[-1].split('_')[0]
    gb = open(gb).read().split('\nLOCUS   ')
    print(gb_name, len(gb))
    
    for plp in gb:
        plp = plp.split('\n')
        accession = [p for p in plp if 'ACCESSION   ' in p]
        source = [p for p in plp if 'SOURCE      ' in p]
        if len(accession) != 1 or len(source) != 1:
            print('GB problems', plp[:2], accession, source, '--\n\n', sep = '\n')
        else:
            accession = accession[0].split('ACCESSION   ')[-1].split(' ')[0]
            source = source[0].split('SOURCE      ')[-1]
            source = ' '.join(source.split(' ')[:2])
            protein_ids = [p.split('protein_id="')[-1].split('"')[0] for p in plp if '/protein_id=' in p]
            for p in protein_ids:
                prot_db[p] = [gb_name, accession, source]

print('prot_db', len(prot_db))


cd_clstr = open(cd_clstr).read().strip().split('\n>')
improve_names = []

compare_nums = []
fut_replace = {}
clstr_n = 0
for clstr in cd_clstr[:]:
    clstr_n += 1
    clstr = clstr.split('\n')
    
    add = {}
    func = []
    main_name = ''
    
    for c in clstr:
        if '>' in c and '...' in c:
            nid = c.split('>')[-1].split('...')[0]
            if nid not in all_prots:
                print('NOT FOUND', nid)
            else:
                add[nid] = [all_prots[nid], c.split('...')[-1].strip()] + prot_db[nid]
                func.append(all_prots[nid])
            
            if '*' in c and 'at' not in c:
                if main_name:
                    print('TWO NAMES?!!!', main_name, nid)
                main_name = nid
            
            compare_nums.append(nid)
            if nid not in prot_db:
                print('LOST IN THE PROT DB', nid)
                
    if not main_name:
        print('NO NAME?!', clstr[0])
    
    func = list(set(func))
    
    tail = [t for t in func if 'tail' in t.lower()]
    if tail:
        tail = 'tail'
    else:
        tail = '-'
        
    conjug = [t for t in func if 'conjug' in t.lower() or ('plasmid' in t.lower() and 'transfer' in t.lower())]
    if conjug:
        conjug = 'conjug'
    else:
        conjug = '-'
        
    phage = [t for t in func if 'phage' in t.lower()]
    if phage:
        phage = 'phage'
    else:
        phage = '-'
        
    func = [f.split('~~~') for f in func]
    trusted = [f for f in func if f[0]]
    genes = [f for f in func if f[1] and f not in trusted]
    others = [f for f in func if f not in trusted and f not in genes]
    
    main_func = all_prots[main_name].split('~~~')
    
    name_it = []
    if trusted:
        status = 'trusted'
        name_func = trusted[0]
        options = len(trusted)
    elif genes:
        status = 'genes'
        name_func = genes[0]
        options = len(genes)
    else:
        status = 'no_name'
        name_func = main_func
        options = len(others)
        
    variants = 'best_option'
    if options > 1:
        variants = f'{options} options'
        
    flag = 'cluster_name'
    if name_func != main_func:
        flag = 'from_inside_cluster'
        
    full_name = '~~~'.join(name_func)
    short_name = name_func[1]
    if not short_name:
        short_name = [s for s in name_func[-1].split(' ') if s != 'DNA' and (len(s) == 3 or len(s) == 4) and len(set(list(s))) != 1]
        if name_func[-1] == 'death on curing protein':
            short_name = 'doc'
            
        if short_name:
            short_name = short_name[-1]
            short_name = '~~~' + short_name + '~~~' + name_func[-1]
        else:
            short_name = name_func[-1]
            short_name = '~~~~~~' + short_name
    else:
        short_name = name_func[0] + '~~~' + name_func[1] + '~~~' + name_func[2]
        
    short_name = short_name[:1].lower() + short_name[1:]
    
    
    fut_replace_name = short_name
    for add_class in (tail, conjug, phage):
        if add_class != '-':
            fut_replace_name += ('___' + add_class)
    fut_replace[main_name] = fut_replace_name 
    
        
    tail_trust = '-'
    if tail == 'tail':
        if 'tail' in full_name.lower():
            tail_trust = 'tail_pos'
        else:
            tail_trust = 'problems_with_tail'
    
        
    name_it = [str(clstr_n), main_name, status, short_name, full_name, variants, flag, phage, conjug, tail, tail_trust, str(len(add))]
    
    add_group = [[main_name] + add[main_name]]
    for a in add:
        if a != main_name:
            add_group.append([a] + add[a])
            
    add_group = ['\t'.join(a) for a in add_group]
    # print(*add_group, sep = '\n')
    
    name_it.append(add_group)
    improve_names.append(name_it)
    
    

compare_nums = list(set(compare_nums))
print('compare_nums', len(compare_nums))

'''
#only hypothetical proteins wer not added by prokka-genbank_to_fasta_db
count = 0
for p in prot_db:
    if  count < 10 and p not in compare_nums:
        count += 1
        print(p, *prot_db[p])
        '''

print('improve_names', len(improve_names))

with open(main + 'cd_hit_annotation_table.tsv', 'w') as fh:
    header = ['N', 'ID_in_the_database', 'Gene_name_trust', 'Gene', 'Gene_descr', 'Consider_option', 'Belong_to_the_ID', 'Phage_related', 'Conjugation', 'Tail_protein', 'Tail_trusted', 'Cluster_size', 'ID', 'description', 'similarity', 'PLP_group', 'Host_accession' ,'Host']
    filler = '.\t' * (len(header) - 6)
    fh.write('\t'.join(header) + '\n')
    for imp in improve_names:
        add_later = imp[-1]
        
        fh.write('\t'.join(imp[:-1]) + '\n')
        for a in add_later:
            fh.write(filler + a + '\n')
    fh.write('\n')
    

overwrite = main + 'PLP_genbanks'
overwrite = open(overwrite).read().strip().split('\n')
with open(main + 'PLP_genbanks_upd', 'w') as fh:
    for line in overwrite:
        if line[:1] == '>':
            pid = line[1:].split(' ')[0]
            new_line = f'>{pid} {fut_replace[pid]}\n'
            fh.write(new_line)
        else:
            fh.write(line + '\n')

    