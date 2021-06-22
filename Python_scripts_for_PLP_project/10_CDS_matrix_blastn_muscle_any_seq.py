#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 17:17:57 2020

@author: lutra

https://stackoverflow.com/questions/62894563/add-3-or-more-legends-to-a-seaborn-clustermap
"""
import os
import glob
import sqlite3

import pandas as pd
import seaborn as sns
import scipy

from pylab import rcParams
import matplotlib
import matplotlib.pylab as plt
import matplotlib.patches as mpatches

from matplotlib.pyplot import gcf

import PLP_main_functions as myf




def blastn_cds(ref, plasmid, outp):
    cmd = f'blastn -task blastn -query {ref} -subject {plasmid} -out {outp} -outfmt "6 qaccver saccver qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"'
    myf.run(cmd)
    return


def ResFinder(txt):
    if txt != '0':
        txt = txt.strip().split('\n')
        save_size = len(txt)
        txt = [t.split(' ') for t in txt]
        arg_classes = list(set([t[0] for t in txt]))
        order = [o for o in ('colistin', 'beta-lactam', 'aminoglycoside', 'quinolone', ) if o in arg_classes  ]
        arg_classes = [a for a in arg_classes if a not in order]
        arg_classes.sort()
        arg_classes = order + arg_classes
        
        args = []
        priority = {'NDM': 1, 'KPC': 1, 'CTX-M':2, 'SHV':3, 'TEM':3}
        for arg_cl in arg_classes:
            curr_class = list(set([t[1] for t in txt if t[0] == arg_cl and t[1] not in args]))
            curr_class = ['_'.join(a.split('_')[:-2]) for a in curr_class]
            curr_class = [tt.strip().replace('_)', ')').replace('_)', ')').replace('(', '').replace(')', '').replace('_', '').split('-Hangzhou')[0] for tt in curr_class]
            
            curr_class = [tt + '@' if '#' in tt or '**' in tt else tt for tt in curr_class]
            curr_class = [tt.replace('**', '').replace('#', '').replace('@', '#') for tt in curr_class]
            curr_class = [tt.replace('bla', '') if tt[:3] == 'bla' else tt for tt in curr_class]
            
            prior = []
            if arg_cl == 'beta-lactam':
                for curr in curr_class:
                    check = [p for p in priority if p in curr]
                    if check:
                        prior.append([curr, priority[check[0]]])
            prior.sort(key = lambda x: x[1])
            args += [p[0] for p in prior]
            add_args = [c for c in curr_class if c not in args]
            add_args.sort()
            args += add_args
    else:
        args=['-']
        
    # print(args)
    if False and args!=['-'] and len(args) != save_size:
        print('\n---')
        print('SIZE PROBLEM!!!', save_size, len(args))
        print(*args)
        check = [t[1] for t in txt]
        check.sort()
        print(*check, sep = '\n')
        print('---\n\n')
    
    args = '_'.join(args)
    # if args != '-':
    #     print(args)
    return args


main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
main = '/data/Current_work/Phage_like_plasmids/PLP_final/'

proj = 'SSU5'
plas_fold = f'{main}prophages_groups/{proj}/'

proj_params = {'P1': (30, 40), 'D6': (30, 15), 'SSU5': (35, 70)} #'P1': (30, 5)


        
rcParams['xtick.bottom'] = rcParams['xtick.labelbottom'] = False
rcParams['xtick.top'] = rcParams['xtick.labeltop'] = True
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
plt.rcParams["figure.figsize"] = proj_params[proj]
sns.set_context("paper")
sns.set_style('whitegrid')
plt.rc('legend', fontsize='24', title_fontsize='32')

fasta = plas_fold + 'fasta/'
fasta = glob.glob(fasta + '*')

table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'
db = main + table + '.sqlite3'

conn = sqlite3.connect(db)
cur = conn.cursor()

Inc_types = {'D6': [''], 'P1': ['IncY_1__K02380', 'p0111_1__AP010962'], 'SSU5': ['IncFIB_pKPHS1_1__CP003223', 'IncFIB_pLF82-PhagePlasmid_1__CU638872', 'IncFIB_H89-PhagePlasmid_1__HG530657', 'IncFIB_pHCM2_1__AL513384', '']}

new_reps = open(main + 'New_replicons.txt').read().strip().split('\n')
new_reps = [n.split('\t') for n in new_reps]

data = []
task = 'SELECT nucleotide, organism, completeness, genome, slen, PLP_status, enterobacteriaceae, D6_putative_replicon_orf42, ResFinder, Major_replicon_variant FROM ' + table
for row in cur.execute(task):
    gid, host, complete, genome, slen, status, enterobacteriaceae, d6_orf42, resfinder, major_var = [str(r) for r in row]
    major_var = major_var.replace('\n', '_')
    if major_var == 'None':
        major_var = 'other'
    if complete == 'complete' and int(slen) < 2000000 and f'_{proj}' in status:
        if genome != 'plasmid': genome = 'phage'
        
        inc = [i for i in Inc_types[proj] if i in enterobacteriaceae][0]
        inc_for_num = inc
        if not inc:
            if 'D6_putative_replicon_orf42' in d6_orf42:
                inc = 'D6_orf42'
            else:
                found = [n[1] for n in new_reps if gid == n[0]]
                if found:
                    inc = found[0]
                    print(gid, inc)
                else:
                    print('Replicon is not found!!!', gid)
        
        num = Inc_types[proj].index(inc_for_num)*2
        if genome == 'phage': num += 1
        
        res = ResFinder(resfinder)
        if res != '-':
            # print(res)
            res = ' ' + res.replace('_', ' ')
        else:
            res = ''
        
        data.append([num, inc, genome, res, gid, host, major_var])
        
data.sort(key = lambda x: (x[0], x[5]))
print(*data[-10:], sep = '\n')



heatmap_fold = plas_fold + 'Heatmap_CDSs/'
list_cdss = [g.split('/')[-1] for g in glob.glob(heatmap_fold + '*') if '.' not in g]

list_cdss = [l for l in list_cdss if 'AF234172' in l or 'MF356679' in l or 'JQ965645' in l]

for ref in list_cdss[:]:
    print(ref)
    
    heatmap = f'{heatmap_fold}{proj}_{ref}_CDSs_heatmap_clustermap.tsv'
    
    cds = heatmap_fold + ref + '/'
    ref_cds = [c for c in glob.glob(cds + '*') if 'CDS' in c and ('fasta' in c or 'fsa' in c)][0]
    print(ref_cds)
    
    blastn_fold = cds + 'blastn/'
    print(blastn_fold)
    myf.mkdir(blastn_fold)
    
    
    if 0:
        print('Do blastn')
        count = 0
        for plasmid in fasta[:]:
            count += 1
            name = plasmid.split('/')[-1].split('.')[0]
            
            if count % 5 == 0:
                print(str(round(100*count/len(fasta), 1)) + '%', name)
            
            prolonged_query = plas_fold + 'long_query.fasta'
            
            query_seq = open(plasmid).read().strip().split('\n')
            dna = ''.join(query_seq[1:])
            
            dna += dna[:10000]
            with open(prolonged_query, 'w') as fh:
                fh.write(f'{query_seq[0]}\n{dna}')
            
            output = blastn_fold + name + '.txt'
            blastn_cds(ref_cds, prolonged_query, output)
            os.remove(prolonged_query)
    
    
    res_cds = glob.glob(blastn_fold + '*')
    
    cds_list = open(ref_cds).read().strip().split('\n')
    
    cds_list = [c[1:].replace(']', '').split(' [') for c in cds_list if '>' in c]
    
    cds_names = []
    name_dict = {}
    for c in cds_list:
        gene = [cc for cc in c if 'protein=' in cc]
        if gene: 
            gene = gene[0].split('=')[1]
        else:
            gene = ''
            
        loc = [cc for cc in c if 'location=' in cc]
        if loc:
            loc = loc[0].split('=')[1]
        else:
            loc = ''
        c_name = '_'.join([c[0].split('_')[-1], gene.upper(), loc])
        cds_names.append(c_name)
        name_dict[c_name] = c[0]
    
    if 1:
        print('collect blastn results into a heatmap')
        
        print(len(cds_names), *cds_names[:5])
        
        header = ['Index', 'Host_species', 'Major_replicons', 'Major_replicon_variant'] + cds_names
        collect = [header]
        
        for entry in data[:]:
            num, inc, genome, resistance, gid, host, major = entry
            res = blastn_fold + gid + '.txt'
            
            if genome != 'phage' or True:
                host = ' '.join(host.split()[:2])
                
            gid_index = gid + '   ' + host
            if genome == 'phage':
                gid_index += ('   ' + ' '.join(host.split()[2:]))
            
            while len(gid_index) < 44:
                gid_index += ' '
            
                
            host = host.replace('virus', 'phage')
            add = [gid_index + resistance, host, inc, major]
            
            if res not in res_cds:
                print(gid, 'res file is missing!')
            else:
                res = open(res).read().strip().split('\n')
                res = [r.split('\t') for r in res]
                for c in cds_list[:]:
                    res_c = [r for r in res if r[0] == c[0]]
                    similarity = 0
                    if res_c:
                        # print(*res_c, sep = '\n')
                        res_c.sort(key = lambda x: float(x[-2]), reverse = True)
                        res_c = res_c[0]
                        # print(res_c)
                        qlen, slen, pident, length = [float(r) for r in res_c[2:6]]
                        # print(qlen, pident, length)
                        
                        if length > qlen:
                            length = qlen
                        
                        similarity = length/qlen * pident
                        # if genome == 'phage':
                        #     similarity *= -1
                    add.append(similarity)
                collect.append(add)
        
        almost_removed = []
        for i in range(4, len(collect[0])):
            curr_cds = collect[0][i]
            curr_sim = [abs(c[i]) for c in collect[1:]]
            curr_sim = [c for c in curr_sim if c >= 90*0.6]
            presence_perc = len(curr_sim)/len(collect[1:])
            if presence_perc < 0.5:
                almost_removed.append([curr_cds, len(curr_sim), round(presence_perc,3)])
         
            
        print(len(almost_removed))
        for a in almost_removed:
            print(*a)

            
        with open(heatmap, 'w') as fh:
            collect = ['\t'.join(map(str, c)) for c in collect]
            fh.write('\n'.join(collect) + '\n')
            
        print('--\n\n\n')
    
    
    if 0:
        print('draw a heatmap')
        df = pd.read_csv(heatmap, index_col = 'Index', sep = '\t')
        print(df.head())
        print(df.shape) 
        df = df.drop(columns = ['Host_species', 'Major_replicons', 'Major_replicon_variant'])       
    
        sns.heatmap(df, cbar=True, cbar_kws = dict(shrink=0.25, ), annot=False, cmap='RdBu_r')
    #    plt.xlabel(x_name, fontname="Arial Bold", fontweight="bold", fontsize = "16")
    #    plt.ylabel(y_name, fontname="Arial Bold", fontweight="bold", fontsize = "16")
        plt.xticks(fontname="Arial", fontsize = "14")
        plt.yticks(fontname="Arial", fontsize = "14")
        plt.savefig(heatmap[:-4] + '.png', bbox_inches='tight', dpi=300)
        plt.show()
        
        
    if 0:
        print('draw a clustermap')
        df = pd.read_csv(heatmap, index_col = 'Index', sep = '\t')
        print(df.head())
        print(df.shape)
        
        if proj == 'P1' and False:
            my_index = list(df.index.values)
            print('orig index length', len(my_index))
            to_drop = [i for i in my_index if 'mcr-' not in i and 'AF234172' not in i]
            df = df.drop(to_drop)
        
        species = df.Host_species
        replicons = df.Major_replicons
        major_vars = df.Major_replicon_variant
        df = df.drop(columns = ['Host_species', 'Major_replicons', 'Major_replicon_variant'])
        print(df.head())
        # print(df.min(), df.max())
        
        
        #add colors to rows and columns
        lut = {}
        lut['Escherichia coli'] = 'red'
        lut['Escherichia marmotae'] = 'red'
        lut['Escherichia fergusonii'] = 'red'
        lut['Escherichia albertii'] = 'red'
        lut['Shigella sonnei'] = 'red'
        
        lut['Kluyvera cryocrescens'] = 'pink'
        
        lut['uncultured bacterium'] = 'grey'
        
        lut['Salmonella enterica'] = 'green'
        lut['Salmonella sp.'] = 'green'
        
        lut['Klebsiella pneumoniae'] = 'blue'
        lut['Klebsiella grimontii'] = 'blue'
        lut['Klebsiella variicola'] = 'blue'
        lut['Klebsiella oxytoca'] = 'blue'
        
        lut['Citrobacter freundii'] = 'aqua'
        lut['Citrobacter sp.'] = 'aqua'
        
        lut['Enterobacter cloacae'] = 'maroon'
        lut['Enterobacter hormaechei'] = 'maroon'
        
        
        lut['Yersinia pestis'] = 'yellow'
        
        lut['Cronobacter sakazakii'] = 'purple'
        
        lut['Escherichia phage'] = 'red'
        lut['Klebsiella phage'] = 'blue'
        lut['Enterobacteria phage'] = 'orange'
        lut['Salmonella phage'] = 'green'
        
        species_colors = species.map(lut)
        
        
        #concise list of species for the legend
        species_lut = {'Escherichia sp.': ['red', ['Escherichia coli','Escherichia marmotae', 'Escherichia fergusonii', 'Escherichia albertii']]} 
        species_lut['Shigella sonnei'] = ['red', []]
        species_lut['Escherichia phage'] = ['red', []]
        
        species_lut['Klebsiella sp.'] = ['blue', ['Klebsiella pneumoniae','Klebsiella grimontii','Klebsiella variicola','Klebsiella oxytoca']]
        species_lut['Klebsiella phage'] = ['blue', []]
        
        species_lut['Salmonella sp.'] = ['green', ['Salmonella enterica']]
        species_lut['Salmonella phage'] = ['green', []]
        
        species_lut['Yersinia pestis'] = ['yellow', []]
        species_lut['Citrobacter sp.'] = ['aqua', ['Citrobacter freundii']]
        species_lut['Enterobacter sp.'] = ['maroon', ['Enterobacter cloacae','Enterobacter hormaechei']]
        species_lut['Cronobacter sakazakii'] = ['purple', []]
        species_lut['Kluyvera cryocrescens'] = ['pink', []]
        species_lut['uncultured bacterium'] = ['grey', []]
        species_lut['Enterobacteria phage'] = ['orange', []]
        
        
        
        palette_size = len(replicons.unique())
        colors = ['orange', 'maroon', 'yellow', 'green', 'purple', 'pink', 'aqua', 'blue'][:palette_size]
        print(colors)
        
        replicons_lut = dict(zip(replicons.unique(), colors))
        print(replicons_lut)
        rep_colors = replicons.map(replicons_lut)
        
        colors = list(matplotlib.colors.CSS4_COLORS)
        colors = [colors[i] for i in range(len(major_vars.unique()))]
        major_vars_lut = dict(zip(major_vars.unique(), colors))
        major_vars_colors = major_vars.map(major_vars_lut)
        
        node_colors = pd.DataFrame(species_colors).join(pd.DataFrame(rep_colors)).join(pd.DataFrame(major_vars_colors))
        node_colors = pd.DataFrame(rep_colors).join(pd.DataFrame(species_colors))
        
        #draw
        #col_colors=rep_colors, 
        #row_linkage=
        my_method = ['ward', 'average', 'single', 'complete'][0]
        my_robust = True
        plot = sns.clustermap(df, figsize=proj_params[proj], method=my_method, metric='euclidean', 
                              robust=my_robust, col_cluster = False, row_colors=node_colors, cmap="coolwarm", 
                              cbar_pos=(0, .8, .03, .2), vmin=0, vmax=100)
        
        
        #specify linewidth in Seaborn's clustermap dendrograms
        for a in plot.ax_row_dendrogram.collections:
            a.set_linewidth(4)
        
        #add legend for replicons
        for label in replicons_lut:
            plot.ax_col_dendrogram.bar(0, 0, color=replicons_lut[label], label=label, linewidth=0);
        l1 = plot.ax_col_dendrogram.legend(title='Major replicons', loc="upper left", 
                                           bbox_to_anchor=(0.85, 1), bbox_transform=gcf().transFigure, frameon = False)
        
        #add legend for species
        shorten_list = set(list(species))
        print(shorten_list)
        save = []
        for label in species_lut:
            check = [label] + species_lut[label][1]
            check = [c for c in check if c in shorten_list]
            save += check
            if check:
                plot.ax_row_dendrogram.bar(0, 0, color=species_lut[label][0], label=label, linewidth=0);
        l2 = plot.ax_row_dendrogram.legend(title='Host species', loc="upper left", ncol=2, 
                                           bbox_to_anchor=(0.85, 0.9), bbox_transform=gcf().transFigure, frameon = False)
        
        
        if [s for s in shorten_list if s not in save]:
            print('LOST SPECIES!!!', [s for s in shorten_list if s not in save])
        else:
            print('ALL SPECIES WERE FOUND))')
        
        plot.savefig(heatmap[:-4] + '_rep1.png', dpi=300, format='png')
