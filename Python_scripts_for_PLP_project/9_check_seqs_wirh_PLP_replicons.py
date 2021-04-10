#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 17:35:32 2020

@author: lutra
"""
import glob
import sqlite3
from PLP_main_functions import read_fasta_file


PC_lab = True

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
else:   
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'


print(main)

proj = 'SSU5'
plas_fold = f'{main}prophages_groups/{proj}/'


plas_reps = {}
db_orig = []

table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'

database = main + table + '.sqlite3'
conn = sqlite3.connect(database)
cur = conn.cursor()


if proj != 'D6':
    rep_column = 'enterobacteriaceae'
else:
    rep_column = 'D6_putative_replicon_orf42'
    
ref_cov = {'D6': 'MF356679_D6_ref', 'P1': 'AF234172_P1_ref', 'SSU5': 'JQ965645_SSU5_ref'}

task = f'SELECT nucleotide, completeness, slen, PLP_status, {rep_column}, {ref_cov[proj]}_cov, {ref_cov[proj]}_CDS_N FROM {table}'
for row in cur.execute(task):
    gid, complete, slen, status, reps, ref_cov, cds_n = [str(r) for r in row]
    label = 'included'
    if complete == 'complete' and int(slen) < 2000000 and f'_{proj}' in status:
        label = 'grouped'
        
    reason = 'grouped'
    if complete != 'complete':
        reason = 'not complete'
    elif int(slen) >= 2000000:
        reason = 'chromosome'
    elif float(ref_cov) < 40:
        reason = float(ref_cov)
    db_orig.append([gid, reason, int(cds_n)])
    
    if reps != '0':
        reps = [r.split()[0].replace('#','').replace('**','') for r in reps.split('\n')]
    else:
        reps = []
    plas_reps[gid] = [label, reps]


check_PLP_size = [p for p in plas_reps if plas_reps[p][0] == 'grouped']
print(proj, 'size', len(check_PLP_size), '\n')

rep_sizes = {}
orig_file = main + 'db/enterobacteriaceae.fsa'
if proj == 'D6':
    orig_file = main + 'db/D6_putative_replicon_orf42.fsa'

rep_dict = read_fasta_file(orig_file)
for r in rep_dict:
    rep_dict[r] = len(rep_dict[r])
    

work_file = main + '0_blastn_26Aug20/'
work_file = glob.glob(work_file + '*')
if proj != 'D6':
    work_file = [w for w in work_file if f'{proj}_all_blastn' in w][0]
else:
    work_file = [w for w in work_file if f'{proj}_orf42_blastn' in w][0]
print(work_file)

work = open(work_file).read().strip().split('\n')
work = [w.split() for w in work if w[:1]!='#']


good = {}
for w in work:
    rep, gid, pident, length = w[:4]
    qlen = rep_dict[rep]
    flag = True
    
    if proj == 'SSU5':
        flag = False
        my_reps = ['pHCM2', 'pLF82-PhagePlasmid', 'H89-PhagePlasmid', 'pKPHS1']        
        if [m for m in my_reps if m in rep]:
            flag = True
            
    if flag and float(pident) >= 90 and float(length)/qlen > 0.6:
        if rep not in good:
            good[rep] = []
        good[rep].append(gid)

answer = []

rep_class = []
reason = []
collect_all = []
for rep in good:
    print(rep)
    contain_replicons = [c.split('.')[0] for c in good[rep]]
    duplicates = [c for c in contain_replicons if contain_replicons.count(c)>1]
    
    contain_replicons = list(set(contain_replicons))
    rep = rep.replace('(', '_').replace(')', '')
    
    print(rep, '\n')
    
    found_grouped = [c for c in contain_replicons if c in plas_reps and plas_reps[c][0] == 'grouped']
    found_included = [c for c in contain_replicons if c in plas_reps and plas_reps[c][0] == 'included']
    missed = [c for c in contain_replicons if c not in plas_reps]
    
    db = [d for d in db_orig if d[0] in found_included]
    db_range = [d[1:] for d in db if d[1] not in ('not complete', 'chromosome')]
    db_range_cov = [d[0] for d in db_range]
    db_range_cds = [d[1] for d in db_range]
    if db_range_cds:
        if len(db_range_cds) > 1:
            db_range = f'qcovs: {int(min(db_range_cov))}-{int(max(db_range_cov))}, CDSs: {min(db_range_cds)}-{max(db_range_cds)}, n={len(db_range_cov)}'
        else:
            db_range = f'qcovs: {int(min(db_range_cov))}, CDSs: {min(db_range_cds)}, n=1'
        
        print('db_range', len(db_range_cds), db_range, *db, '\n* * *\n\n')
        
    db = [d[1] if d[1] in ('not complete', 'chromosome') else db_range for d in db]    
    db = [d for d in db if d == db_range] + [d for d in db if d == 'chromosome'] + [d for d in db if d == 'not complete']
    
    rep_class += [rep[:26] for d in db]
    reason += db
    
    rep_name = rep.split('_1')[0].replace('IncFIB_', 'FIB_').replace('-PhagePlasmid', '').replace('H89', 'pH89')
    if rep_name[:2] == 'D6':
        rep_name = 'D6_orf42'
    collect_all += [rep_name + '__' + d for d in db]
    collect_all += [rep_name + '__' + proj + '-PLP group' for g in found_grouped]
    
    answer += missed
    
    check = sum([len(found_grouped), len(found_included), len(missed)])
    if check != len(contain_replicons):
        print('SIZE CHECK PROBLEM', check)
        
    existed_replicons = [p for p in plas_reps if rep in plas_reps[p][1]]
    existed_not_blasted = [e for e in existed_replicons if e not in contain_replicons]
    
    missed_in_exist = [c for c in contain_replicons if c not in missed and c not in existed_replicons]
    
    if duplicates:
        duplicates = list(set(duplicates))
        print('Duplicates exist!', len(duplicates), *duplicates, '\n')
    
    print(len(contain_replicons), *contain_replicons[:5])
    print('Found: grouped -', len(found_grouped), ', not grouped -',  len(found_included))
    print('Missed - ', len(missed), '\n')
    print('Exist in DB:', len(existed_replicons))
    print('Exist&Grouped -', len([e for e in existed_replicons if e in found_grouped]), ', Exist&Not grouped -', len([e for e in existed_replicons if e in found_included]), ', not found by blast -', len(existed_not_blasted))
    if missed_in_exist:
        print('\nMissed in exist!', len(missed_in_exist), *missed_in_exist)
    print('--\n\n')
    
if answer:
    with open(f'{main}Entrez_add_replicons_for_{proj}.txt', 'w') as fh:
        fh.write('\n'.join(answer) + '\n')


print('collect_all', len(collect_all))
collect_all_clust = list(set(collect_all))
collect_all_clust.sort()
collect_all_clust.sort(key = lambda x: (x.split('__')[-1]))

names = []
colors = []
sizes = []

color_lib = {'PLP group': 'green', 'chromosome': 'blue', 'not complete': 'grey', 'CDSs': 'orange'}

print('N of clusters', len(collect_all_clust))
for clust in collect_all_clust:
    print(clust, )
    
    size = collect_all.count(clust)
    color = ''
    for col in color_lib:
        if col in clust:
            if not color:
                color = color_lib[col]
            else:
                print('COLOR EXISTS', color, col, clust)
    if not color:
        print('COLOR LOST', clust)
    
    if size > 1:
        name = clust.split('__')[0] + '\nn=' + str(size)
    else:
        name = clust.split('__')[0] + ', ' + str(size)
    
    if 'CDSs' in clust:
        name = clust.replace('__', '\n')
        name = name.replace(', CDSs', '%\nCDSs')
        # name = name.replace(', CDSs', '%, CDSs')
        name = name.replace(', n=', '\nn=')
        
        name, qcovs, cdss, n = name.split('\n')
        name = name + ', ' + n + '\n' + qcovs + '\n' + cdss
        
        
    # elif proj == 'SSU5':
    #     name = name.replace(', n=', '\nn=')
        
        
    names.append(name)
    colors.append(color)
    sizes.append(size)
    
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import squarify # pip install squarify (algorithm for treemap)
from pylab import rcParams

# if proj == 'P1':
#     rcParams['figure.figsize'] = 8,6
#     SMALL_SIZE = 14
# elif proj == 'SSU5':
rcParams['figure.figsize'] = 6,7
SMALL_SIZE = 14
    

plt.rcParams["font.family"] = "Arial"
MEDIUM_SIZE = SMALL_SIZE + 0
BIGGER_SIZE = SMALL_SIZE + 1

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
 
# Create tree plot
alpha_level = 0.75
squarify.plot(sizes=sizes, label=names, color=colors, alpha=alpha_level, pad=True)
plt.axis('off')

#add user legend
legend_dict = { 'PLP_group': 'green', 'Chromosome': 'blue', 'Complete, lower coverage': 'orange', 'Not complete sequences': 'grey' }
patchList = []
for key in legend_dict:
    data_key = mpatches.Patch(color=legend_dict[key], label=key)
    patchList.append(data_key)
leg = plt.legend(handles=patchList, bbox_to_anchor=(1, 0.24), loc='upper left', frameon=False)

#add opacity to colors in the legend
for lh in leg.legendHandles: 
    lh.set_alpha(alpha_level)


fig_name = {'SSU5': 'SSU5 four FIB replicons*', 'P1': 'P1-group IncY and p0111 replicons', 'D6': 'D6_repB'}

plt.title(f'{fig_name[proj]}, n = {len(collect_all)}', y=-0.08)

plt.savefig(f'{main}Charts/{proj}_replicons.svg', format='svg', bbox_inches='tight')

plt.show()
            

if 0:       
    import pandas as pd
    import altair as alt
    
    
    if len(rep_class) != len(reason):
        print('HORRIBLE LEN PROBLEM', len(rep_class), len(reason))
    else: 
        d = {'Replicon': rep_class, 'Reason of exclusion': reason}
        df = pd.DataFrame(data=d)
        
        
        chart = alt.Chart(df).mark_bar(
            cornerRadiusTopLeft=3,
            cornerRadiusTopRight=3
        ).encode(
            x='Replicon',
            y='count(Replicon):Q',
            color='Reason of exclusion'
        )
        chart.show()
