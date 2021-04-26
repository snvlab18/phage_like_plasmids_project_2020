#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 17:35:32 2020

@author: lutra
"""
import sqlite3


PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
else:   
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'


print(main)

proj = 'SSU5'

table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'

database = main + table + '.sqlite3'
conn = sqlite3.connect(database)
cur = conn.cursor()

if proj != 'D6':
    rep_column = 'enterobacteriaceae'
else:
    rep_column = 'D6_putative_replicon_orf42'
    
ref_cov = {'D6': 'MF356679_D6_ref', 'P1': 'AF234172_P1_ref', 'SSU5': 'JQ965645_SSU5_ref'}
replicons = {'D6': ['D6_putative_replicon_orf42_c45699_46596'], 'P1':['IncY_1__K02380', 'p0111_1__AP010962'], 'SSU5': ['IncFIB_pHCM2_1__AL513384', 'IncFIB_pKPHS1_1__CP003223', 'IncFIB_H89-PhagePlasmid_1__HG530657', 'IncFIB_pLF82-PhagePlasmid_1__CU638872']}

collect_all = []

for rep in replicons[proj]:
    task = f'SELECT nucleotide, completeness, slen, PLP_status, {rep_column}, {ref_cov[proj]}_cov, {ref_cov[proj]}_CDS_N FROM {table}'
    task += f" WHERE {rep_column} LIKE '%{rep}%'"
    
    qcovs_ranges = []
    CDSs_ranges = []
    
    rep_collect = []
    
    for row in cur.execute(task):
        gid, complete, slen, status, reps, cov, cds_n = [str(r) for r in row]
        cov = float(cov)
        if complete == 'complete' and int(slen) < 2000000 and cov:
            label = 'grouped'
        
        reason = ''
        if complete == 'complete' and int(slen) < 2000000 and cov >= 40:
            reason = 'grouped'
        elif complete != 'complete':
            reason = 'not complete'
        elif int(slen) >= 2000000:
            reason = 'chromosome'
        elif cov < 40:
            reason = 'lower_coverage'
            qcovs_ranges.append(cov)
            CDSs_ranges.append(int(cds_n))
            
        if not reason:
            print('REASON PROBLEM', gid)
        
            
        upd_rep = ''
        if proj == 'D6':
            upd_rep = 'D6_orf42'
        elif proj == 'P1':
            upd_rep = rep.split('_')[0]
        else:
            upd_rep = '_'.join(rep.split('_')[:2])
            
        if reason == 'grouped':
            reason = f'{proj}-PLP group'
            
        rep_collect.append(f'{upd_rep}__{reason}')
     
    db_range = ''
    if qcovs_ranges:
        if len(qcovs_ranges) > 1:
            db_range = f'qcovs: {int(min(qcovs_ranges))}-{int(max(qcovs_ranges))}, CDSs: {min(CDSs_ranges)}-{max(CDSs_ranges)}, n={len(CDSs_ranges)}'
        else:
            db_range = f'qcovs: {int(qcovs_ranges[0])}, CDSs: {CDSs_ranges[0]}, n=1'
        
    rep_collect = [r.replace('lower_coverage', db_range) for r in rep_collect]
    collect_all += rep_collect
        
            

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


fig_name = {'SSU5': 'SSU5 four FIB replicons*', 'P1': 'P1-group IncY and p0111 replicons', 'D6': 'D6_orf42'}

plt.title(f'{fig_name[proj]}, n = {len(collect_all)}', y=-0.08)

plt.savefig(f'{main}Charts/{proj}_replicons.svg', format='svg', bbox_inches='tight')

plt.show()

