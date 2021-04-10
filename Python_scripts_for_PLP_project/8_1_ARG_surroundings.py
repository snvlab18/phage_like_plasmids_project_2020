#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 13:08:49 2020

@author: lutra

annotate - /home/lutra/Dropbox/Kira/scripts/phage_like_plasmids/PLP_Easyfig_master_v2.py
"""
import glob
import sqlite3
from collections import OrderedDict as orddic
from dna_features_viewer import GraphicFeature, GraphicRecord
import PLP_main_functions as myf


PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
else:   
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'


print(main)
    

proj = 'SSU5'
plas_fold = f'{main}prophages_groups/{proj}/'


with_args = main + 'all_ARG_encoding.txt'
with_args = open(with_args).read().strip().split('\n')

not_only_ARG_encoding = False
ARG_extra_Inc = ['CP021720', 'CP055609', 'CP017632', 'CP021203', 'CP033632', 'CP042617', 'KP453775', 'KX928752', 'LT985263', 'MK419152', 'CP017387', 'CP050276', 'CP052423', 'CP052452', 'CP052672', 'CP056516', 'CP058165']
extract_region = True


plas_lengths = {}
table = 'Phage_like_plasmids_SSU5_P1_D6_8Sep20'

database = main + table + '.sqlite3'
conn = sqlite3.connect(database)
cur = conn.cursor()

task = 'SELECT nucleotide, slen FROM ' + table
for row in cur.execute(task):
    gid, slen = [str(r) for r in row]
    plas_lengths[gid] = int(slen)


clusters = plas_fold + 'IDs_grouped_by_Rep.txt'
clusters = open(clusters).read().strip().split('\n')

cluster_fold = plas_fold + 'HomBlocks_per_cluser/'

answer = []
types = []
feat_db = {}
db = orddic()

arg_classes = ['beta-lactam', 'aminoglycoside', 'trimethoprim', 'sulphonamide', 'phenicol', 'tetracycline', 'fosfomycin', 'macrolide', 'colistin', 'quinolone', 'rifampicin']
oriTDB_classes = ['relaxase_all_blastx', 'auxiliary_all_blastx', 'oriT_all', 't4cp_all_blastx']
remove_feats = ['phenicol cml_1_M22614', 'tetracycline tet_A_6_AF534183', 'quinolone qepA_1_FJ167861', 'quinolone aac_6_-Ib-cr_2_EF636461', 'IncFII(pSFO)_1__AF401292', 'IncFII_pSFO_1__AF401292', 'IS15DIV_IS6_820bp', 'IS100X_IS21_1924bp', 'IS100X_IS21_1950bp', 'IS1541A_IS200/IS605_IS200_710bp', 'ISEc9_IS1380_1656bp', 'quinolone aac_6_-Ib-cr_1_DQ303918']

#improve vfdb and oriTDB feat names!  
expand_feats = {}
for oritdb in glob.glob(main + 'db/oriTDB/*') + [main + 'db/VFDB_setB_nt.fsa']:
    oritdb_name = oritdb.split('/')[-1].split('.')[0]
    oritdb_db = [v[1:].strip().replace(' ', '--') for v in open(oritdb).read().strip().split('\n') if v[:1] == '>']
    genes = [v.split('--')[0] for v in oritdb_db]
    gene_names = [v.split('--')[1].replace('[','').replace(']','').replace('(','').replace(')','') if len(v.split('--')) > 1 else v for v in oritdb_db]
    # print(oritdb_name, genes[:10], gene_names[:10], '--\n\n', sep = '\n')
    
    oritdb_dict = dict(zip(genes, gene_names))
    expand_feats[oritdb_name] = oritdb_dict
    
for e in expand_feats:    
    print(e, len(expand_feats[e]))
            
print('--\n\n')

clusters = [c for c in clusters if 'ARG_extra' not in c]


visualize = {}
count = 0
save_dna = {}
for cluster in clusters[:]:
    cluster = cluster.split(' ')
    plas_group = cluster[0]
    IDs = cluster[1:]
    print(proj, plas_group, len(IDs))
    
    annotate = cluster_fold + plas_group + '/annotate/'
    annotate = glob.glob(annotate + '*')
    print('Len ANNOTATE', len(annotate))
    
    backbone = cluster_fold + plas_group + '/backbone/'
    backbone = glob.glob(backbone + '*')
    
    fasta_fold = cluster_fold + plas_group + '/fasta/'
    print(fasta_fold)
    fasta_fold = glob.glob(fasta_fold + '*')
    
    for plas in IDs:
        if (plas in ARG_extra_Inc) and (plas in with_args or not_only_ARG_encoding) and count < 500:
            print('PLASMID!!!', plas)
            count += 1
            # print(plas)
            annot = [a for a in annotate if plas in a.split('/')[-1]]
            backb = [b for b in backbone if plas in b.split('/')[-1]]
            fasta_file = [f for f in fasta_fold if plas in f.split('/')[-1]][0]
            
            if len(annot) != 1 or len(backb) != 1:
                print('FILE PROBLEM', plas, '\nAnnotate:', *annot, '\nBackbone:', *backb)
            else:
                # print(fasta_file)
                this_dna = myf.read_fasta_file(fasta_file)
                upd_dna = {}
                for d in this_dna:
                    print(d, len(this_dna[d]))
                    new_d = d.split('_')[0]
                    upd_dna[new_d] = this_dna[d]
                save_dna[plas] = upd_dna[plas]
                
                annot = [a.split('\t')[:7] for a in open(annot[0]).read().strip().split('\n')]
                
                backb = [b.split('\t')[:5] + [100, 100] for b in open(backb[0]).read().strip().split('\n')]
                annot += backb
                
                vals = [a for a in annot if a[0] in arg_classes + oriTDB_classes + ['enterobacteriaceae', 'D6_putative_replicon_orf42', 'ISel_db', 'VFDB_setB_nt', 'ARG'] or 'backbone' in a[0]]
                # print(*vals,'--\n', sep = '\n')
                
                new_vals = []
                for v in vals:
                    check = False
                    if v[0] in expand_feats:
                        check = True
                        suff = [s for s in ('**', '#',) if s in v[1]]
                        if suff:
                            suff = suff[0]
                            v[1] = v[1].replace(suff, '')
                        else:
                            suff = ''
                        
                        v[1] = suff + expand_feats[v[0]][v[1]]
                        # print(v[0], v[1])
                        
                    if 'decircle' not in v:
                        v = ['_SPACE_'.join(v[:-5])] + v[-5:]
                    else:
                        dec_index = v.index('decircle')
                        last = '_'.join(v[dec_index:])
                        v = v[:dec_index]
                            
                        v = ['_SPACE_'.join(v[:-5]) + '_SPACE_' + last] + v[-5:]
                        
                    # if check:
                    #     print(v, len(v))
                    new_vals.append(v)
                
                vals = new_vals[::]
                try:
                    vals = [[aa[0], int(aa[1]), int(aa[2]), aa[3], float(aa[4]), float(aa[5])] for aa in vals]                       
                except:
                    # print(vals)
                    print('UNK PROBLEMS')
                    vals = []
                vals_upd, report = myf.find_overlaps(vals)
                # print(report, '\n- - -\n')
                
                annot = [v[0].split('_SPACE_') + list(map(str, v[1:])) for v in vals_upd]
                annot = [v[:5] for v in annot]
                
                
                if len(annot[0]) > 1:
                    # print(*annot,'--\n', sep = '\n')
                    annot.sort(key = lambda x: int(x[2]))
                    # print(*annot,'--\n', sep = '\n')
                    collect = []
                    for a in annot:
                        if 'backbone' in a[0]:
                            collect.append(' '.join(a))
                        # elif a[0] in ('oriT_all', 'relaxase_all_blastx', 't4cp_all_blastx', 'auxiliary_all_blastx',):
                            #not interesting gene groups
                            # pass
                        else:
                            if not collect or 'backbone' in collect[-1].split()[0]:
                                collect.append(' '.join(a))
                            else:
                                collect[-1] = collect[-1] + ';' + ' '.join(a)
                    
                    # if  len(collect)>1 and 'backbone' not in collect[0].split()[0] and 'backbone' not in collect[-1].split()[0]:
                    #     print('NEED TO ROTATE', plas + ' ' + plas_group, collect[0], collect[-1], sep = '\n\n')
                    #     collect = collect[1:-1] + [collect[-1] + ';' + collect[0]]
                    #     print('--\n\n')
                    
                    interesting = [c for c in collect if 'backbone' not in c.split()[0]]
                    # print(len(annot), len(collect), len(interesting))
                    size = len(collect)
                    for i in interesting:
                        feat = i.split(';')
                        feat = [f.split() for f in feat]
                        
                        feat = [f for f in feat if f[1].replace('#', '').replace('**', '').replace('(', '_').replace(')', '') not in remove_feats and ' '.join(f[:2]).replace('#', '').replace('**', '') not in remove_feats]
                        
                        
                        arg = [f[0] for f in feat if f[0]]
                        arg = list(set(arg))
                        types += [a for a in arg if a not in types]
                        
                        arg = [a for a in arg if a in arg_classes]                
                        if arg:
                            arg = '"' + '\n'.join(arg) + '"'
                        else:
                            arg = '-'
                        
                        feat_descr = [f[1].split('_')[0] for f in feat]
                        feat_descr = [f if f not in ('D6', '#D6', '**D6') else f + '_orf42' for f in feat_descr]
                        feat_descr = [' '.join(feat_descr[i*5:(i+1)*5]) for i in range(int(len(feat_descr)/5)+1)]
                        feat_descr = [f for f in feat_descr if f]
                        feat_descr = '"' + '\n'.join(feat_descr) + '"'
                        
                        feat_size = int(feat[-1][3]) - int(feat[0][2])
                        # if feat_size > 0:
                        #     print('this feat is', plas, feat_descr)
                        
                        
                        feat_db[feat_descr] = arg
                        
                        if arg != '-':
                            new_name = feat_descr.replace('#', '').replace('**', '').replace('"', '').replace('\n', ' ')
                            reverse_name = ' '.join(new_name.split(' ')[::-1])
                            if reverse_name in visualize:
                                new_name = reverse_name
                            
                            if new_name not in visualize:
                                visualize[new_name] = [[], [], [], []]
                            visualize[new_name][0].append(feat)
                            visualize[new_name][1].append(feat_descr.replace('"', '').replace('\n', ' '))
                            visualize[new_name][2].append(feat_size)
                            visualize[new_name][3].append(plas)
                            
                            # print('***\nto viz', *visualize[new_name], '% % %', *visualize[new_name][0], '\n* * *\n\n', sep='\n')
                            
                        
                        feat_coords = f'feat_coords_{feat[0][2]}{feat[0][4]}..{feat[-1][3]}{feat[-1][4]}'
                        
                        index_i = collect.index(i)
                        
                        prev_i, next_i = index_i - 1, index_i + 1
                        if prev_i == -1:
                            prev_i = size - 1
                        if next_i == size:
                            next_i = 0
                        
                        prev_cds, next_cds = collect[prev_i].split(), collect[next_i].split()
                        
                        prev_name, next_name = prev_cds[1], next_cds[1]
                        prev_name = prev_name.split('lcl')[0] + 'cds' + prev_name.split('_')[-1]
                        next_name = next_name.split('lcl')[0] + 'cds' + next_name.split('_')[-1]
                        
                        cds_pair = prev_name + '_' + next_name
                        cds_coords = f'ref_CDSs_ends_{prev_cds[3]}{prev_cds[4]}..{next_cds[2]}{next_cds[4]}'
                        
                        if cds_pair not in db:
                            db[cds_pair] = orddic()
                        
                        if feat_descr not in db[cds_pair]:
                            db[cds_pair][feat_descr] = []
                            
                        entry = [plas_group.split('_')[0], plas, str(feat_size), feat_coords, cds_coords]
                        db[cds_pair][feat_descr].append(entry)
                    print('* * *\n\n')
# print(types)    


def pick_color(db):
    global arg_classes
    global oriTDB_classes
    
    colors_db = {'PlasmidFinder': '#0000FF', 'ISel_db' : '#FFA500', 'ARG' : '#FF0000', 'oriTDB': '#87CEFA', 'VFDB_setB_nt': '#006600'}
    color = '#000000'
    if db in arg_classes:
        color = 'ARG'
    elif db in oriTDB_classes:
        color = 'oriTDB'
    elif db == 'ISel_db':
        color = 'ISel_db'
    elif db in ['enterobacteriaceae', 'D6_putative_replicon_orf42']:
        color = 'PlasmidFinder'
    elif db == 'VFDB_setB_nt':
        color = 'VFDB_setB_nt'
    
    if color in colors_db:
        color = colors_db[color]
    else:
        print('COLOR PROBLEM', db, color)
        
    return color
        
nap_db = {'+': +1, '-': -1}

if 1:
    print('\n\nVISUALIZE\n')
    count = 0
    for viz in visualize:
        if 'NDM' in viz or True:
            # print(viz, '\n+ + +\n', len(visualize[viz]))
            
            viz_vars, viz_names, viz_lengths, viz_ids = visualize[viz]
            # if len(viz_lengths) > 1:
            #     print(viz, *viz_lengths)
            viz_long = viz_lengths.index(max(viz_lengths))
            
            viz_use = viz_vars[viz_long]
            viz_len = viz_lengths[viz_long]
            viz_title = f'{viz_names[viz_long]} {viz_len} bp'
            
            N = len(viz_ids)
            viz_ids = ' '.join(viz_ids)
            viz_descr = f'N = {N}, {proj}-PLP: {viz_ids}, {viz_len} bp'
            # print(*viz_use, sep = '\n')
            
            min_coord = min([min(map(int, v[2:4])) for v in viz_use])
            max_coord = max([max(map(int, v[2:4])) for v in viz_use])
            print(viz_ids, min_coord, max_coord, 'size', max_coord - min_coord + 1)
            
            
            
            
            
            fut_name = []
            features = []
            for use in viz_use:
                db, name, c1, c2, nap = use
                my_color = pick_color(db)
                
                shorten_name = name.split('_')[0]
                if db in arg_classes:
                    fut_name.append(shorten_name.replace('#', '').replace('**', ''))
                
                feat = GraphicFeature(start=int(c1) - min_coord, end=int(c2) - min_coord, strand=nap_db[nap], color=my_color, label=shorten_name, thickness=14, linewidth = 0.5, box_linewidth = 0.5, fontdict={'family':'Arial'})
                features.append(feat)
                
            fut_name = '_'.join(fut_name)
            
            if extract_region:
                region = save_dna[viz_ids][min_coord-1 : max_coord]
                region = myf.dna_to_fasta_format(region)
                reg_name = f'{proj}_{viz_ids}_{fut_name}_{min_coord}_{max_coord}_size{max_coord - min_coord +1}' 
                with open(f'{main}ARG_maps/cut_ARG_extra_Incs/{reg_name}.fasta', 'w') as fh:
                    fh.write(f'>{reg_name}\n{region}\n')
            
            fut_name += ('_n' + str(N) + '_' + str(count))
            count += 1
            
            
            width = 5
            if len(viz_use) > 5:
                width = len(viz_use)/2
                
            add_width = int(viz_len/1000)
            print(viz_len, width, add_width)
            
            
            record = GraphicRecord(sequence_length=viz_len, features=features)
            ax, _ = record.plot(figure_width=width)
            print(viz_title)
            print(viz_descr)
            print(width, add_width)
            
            # ax.set_xlabel('\n' + viz_title + '\n' + viz_descr, fontdict=dict(weight='bold'))
            ax.set_xlabel('\n' + viz_descr, fontdict=dict(weight='bold'))
            
            # ax.figure.savefig(f'{main}ARG_maps/{proj}/{proj}_{fut_name}.png', bbox_inches='tight', dpi = 300)
            ax.figure.savefig(f'{main}ARG_maps/cut_ARG_extra_Incs/images/{reg_name}.png', bbox_inches='tight', dpi = 300)
            
            print('\n\n')
            



if 0:
    with open(f'{main}ARG_surroundings_{proj}.tsv', 'w') as fh:
        header = ['Ref_CDS_borders', 'Feature', 'ARG_classes', 'Feature_N', 'Plasmid_group', 'Nucleotide', 'Feature_size', 'Feature_coords', 'Ref_CDSs_ends']
        fh.write('\t'.join(header) + '\n')
        
        fut_ans = []
        for ref_cds in db:
            for feat in db[ref_cds]:
                feats = db[ref_cds][feat]
                entry = [ref_cds, feat, feat_db[feat], str(len(feats))]
                
                if len(feats) == 1:
                    entry += feats[0]
                else:
                    entry_combo = []
                    for i in range(len(feats[0])):
                        add = '"' + '\n'.join([f[i] for f in feats]) + '"'
                        entry_combo.append(add)
                    entry += entry_combo
                fut_ans.append(entry)
        fut_ans.sort(key = lambda x: (x[1].replace('#', '').replace('**', ''), x[1]))
        for fut in fut_ans:
            fh.write('\t'.join(fut) + '\n')
        
                
                
                
                
            
            
            