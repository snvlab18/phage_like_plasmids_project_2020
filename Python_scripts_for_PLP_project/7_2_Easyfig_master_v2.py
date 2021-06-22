#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 20:42:18 2020

@author: lutra

colors
https://www.rapidtables.com/web/color/RGB_Color.html

main colors - first 16
"""
import sqlite3
import glob
import re
import PLP_main_functions as myf
import main_annotate


def generate_gb(new_gb, name, host, size, feats, dna, easyfig=False):
    global main
    global dropbox
    global lost_hash
    
    gb=open(new_gb,'w')
    first=f'LOCUS       {name}              {size} bp    DNA\nDEFINITION  {host} {name}, complete sequence\nFEATURES             Location/Qualifiers\n     source          1..{size}\n'
    gb.write(first)
    
    color_code = {}
    colors = f'{dropbox}Kira/iTol_templates/rgb_colors.txt'
    colors = open(colors).read().strip().split('\n')
    colors = [c.split('\t') for c in colors]
    
    for c in colors:
        if easyfig:
            color_code[c[0]] = c[2].replace('(', '').replace(')', '').replace(',',' ') #easyfig
        else:
            color_code[c[0]] = c[0] #genoplotr
        
    
    for f in feats:
        el, c1, c2, strand, color = f
        
        el = el.replace('**', '#')
        if el[:1] == '#':
            el = el[1:] + '#'
        
        sign = ''
        if '#' in el:
            sign = '#'
            
        if 'block' not in el:
            if color == 'dark green' and len(el.split('_')[0]) > 3:
                el = (el.split('_'))[0]
            elif color == 'blue':
                el = el.split(' ')[0].split('-Phage')[0].replace('(', '_').replace(')', '')
                if 'IncFIB_pKPHS1_1__CP003223' in el:
                    el = el.replace('IncFIB_pKPHS1_1__CP003223', 'IncFIB_pKPHS1')
                elif 'IncFIB_pHCM2_1__AL513384' in el:
                    el = el.replace('IncFIB_pHCM2_1__AL513384', 'IncFIB_pHCM2')
                elif 'p0111' in el or 'IncY' in el:
                    el = el.split('_')[0]
            elif color not in ('dark green', 'blue'):
                el = (el.split('_'))[0]
            if el[:3] == 'bla':
                el = el[3:]
        else:
            if not easyfig:
                el = ''
            else:
                pass
        
        if sign not in el:
            if el not in lost_hash:
                lost_hash.append(el)
            el += sign
            
        if color == 'blue':
            print(el)
        
        
        coords = '{}..{}'.format(c1, c2)
        
        if strand == '-':
            coords='complement('+coords+')'
            
        if color not in ():
            entry='     CDS             '+coords+'\n                     /locus_tag="'+el+'"\n                     /colour='+color_code[color]+'\n'
            gb.write(entry)
        else:
            entry='     CDS             '+coords+'\n                     /colour='+color_code[color]+'\n'
            gb.write(entry)
    
    dna = [dna[i*60:(i+1)*60] for i in range(int(len(dna)/60)+1)]
    dna = [d for d in dna if d]
    text='ORIGIN\n'
    count=1
    for f in dna:
        reg=f.strip()
        tab3=' '*(9-len(str(count)))+str(count)+' '
        for i in range(int(len(reg)/10)+1):
            piece=reg[10*i:10*(i+1)]
            tab3=tab3+piece+' '
        tab3+='\n'
        text+=tab3
        count+=60
    gb.write(text)
    gb.write('//\n\n')
    gb.close()
    return


def Easyfig(query, figure, fold, comp):
    easyfig_cmd = 'Easyfig.py -svg -width 8000 -ann_height 250 -blast_height 500 -f1 T -f2 10000 -min_length 1500 -aln left -blast_col 0 191 255 0 0 255 -blast_col_inv 255 215 0 255 140 0 -bo F -f CDS 0 0 0 rect -glt 5 -genet 20 -legend single -leg_name locus_tag -uncomp T -bo F '
    
    if comp:
        choose = 0
    else:
        choose = 1
    
    easyfig_cmd = 'python2.7 {}Easyfig-lutra/{}'.format(['/home/shiri/plasmid_project/tools/', '/data/Bioinformatics/'][choose], easyfig_cmd)
    
    cmd = f'{easyfig_cmd} -o {figure} {query}'
    myf.run(cmd)
    
    for q in query.split()[:0]:
        q_name = q.split('/')[-1].split('.')[0]
        cmd = f'{easyfig_cmd} -o {fold}{q_name}.svg {q}'
        myf.run(cmd)
    return


PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
    dropbox = '/home/shiri/Dropbox/'
else:   
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'
    dropbox = '/home/lutra/Dropbox/'

lost_hash = []
proj = 'D6'
print('PROJECT!', proj, main, '\n')

with_args = main + 'all_ARG_encoding.txt'
with_args = open(with_args).read().strip().split('\n')

not_only_ARG_encoding = True

plas_fold = f'{main}prophages_groups/{proj}/'
cluster_fold = plas_fold + 'HomBlocks_per_cluser/'
clusters = glob.glob(cluster_fold + '*')
clusters = [c +'/' for c in clusters if '.' not in c and c.split('/')[-1] not in ('Easyfig', 'Gepard_dotplot', 'iqtree', 'Viz_GenBank', 'Viz_extra_matches')]
clusters.sort()

# clusters = [c for c in clusters if 'ARG_extra' not in c and 'mcr' not in c and 'rotated_ARG' not in c and 'problem' not in c]
# clusters = [c for c in clusters if 'pHCM' in c]
# clusters = [c for c in clusters if 'rotated_ARG' in c]
clusters = [c for c in clusters if 'other' in c]

tree_fold = cluster_fold + 'iqtree/'

viz_fold = cluster_fold + 'Viz_GenBank/'
# viz_fold = cluster_fold + 'Viz_extra_matches/'
myf.mkdir(viz_fold)

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
    
expand_feats['D6_putative_replicon_orf42'] = {'D6_putative_replicon_orf42_c45699_46596': 'D6_orf42'}
    
for e in expand_feats:    
    print(e, len(expand_feats[e]))
print('--\n\n')


for clust in clusters[:]:
    print(clust)
    fasta = clust + 'fasta/'
    plasmids = glob.glob(fasta + '*')

    gb_fold = viz_fold + 'genbank/'
    blast_fold = viz_fold + 'blast_genoplotr/'
    
    for f in (gb_fold, blast_fold, ):
        myf.mkdir(f)
        
    # mauve_out = [g for g in glob.glob(clust + '*') if g[-10:] == '.mauve.out'][0]
    # group_name = mauve_out.split('/')[-1].split('.')[0]
    group_name = clust.split('/')[-2]
    print(group_name)
    
    tree = tree_fold + group_name + '_extract_backbones.fasta.treefile'
    if tree not in glob.glob(tree_fold + '*'):
        print('No tree!!')
        tree = ''
        
    annotate = clust + 'annotate/'
    backbone = clust + 'backbone/'
    prokka = clust + 'prokka/'
    
    rerun = False
    
    #annotate!!!
    print('annotate first')
    
    db = main + 'db/'
    db = glob.glob(db + '*')
    # db = [d for d in db if 'oriT' not in d]
    db = [f for f in db if '.fsa' in f] + ' '.join([' '.join(glob.glob(f+'/*')) for f in db if '.fsa' not in f]).split()
    db = [f for f in db if '.fsa' in f]
    
    
    print('show files in my local db', *db, '--\n\n', sep = '\n')
    
    count = 0
    # print(*plasmids)
    
    plasmids = [p for p in plasmids if 'CP042620' in p] + [p for p in plasmids if 'CP042620' not in p]
    for query in plasmids:
        count += 1
        if rerun and (query.split('/')[-1].split('_')[0] in with_args or not_only_ARG_encoding):
            print(query, count, str(round(100*count/len(plasmids), 1))+'%')
            main_annotate.annotate_query(query, clust, db, annotate, rerun=rerun)
            
            print('--\n\n')
        
       
    describe = []
    for fold in (annotate, backbone, prokka, ):
        describe += [g for g in glob.glob(fold + '*') if '.gff' not in g]
    print(len(describe), describe[0])
    

    #collect_species
    table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'
    
    database = f'{main}{table}.sqlite3'
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    
    plasmid_host = {}
    task = 'SELECT nucleotide, organism FROM ' + table
    for row in cur.execute(task):
        nt, org = [str(r) for r in row]
        org = ' '.join(org.split(' ')[:2])
        org = org.replace('virus', 'phage')
        
        plasmid_host[nt] = org
        
    matches = '/data/Current_work/Phage_like_plasmids/PLP_final/ARG_maps/cut_ARG_extra_Incs/matches1.txt'
    matches = open(matches).read().strip().split('\n')
    matches = [m.split('\t') for m in matches]
    
    for m in matches:
        plasmid_host[m[3]] = m[4]
    
    
    #start to create color db
    my_colors = {'PlasmidFinder':'blue', 'ISel_db' : 'orange', 'ARG' : 'red', 'oriTDB': 'light blue', 'VFDB_setB_nt': 'dark green', 'CDS': 'gray', 'hypothetical': 'light gray', 'tRNA': 'medium orchid', 'toxin': 'lawn green', 'tail': 'sienna', 'phage': 'dark golden rod'}
    pale_colors = 'pale green,dark khaki,medium aqua marine,khaki,powder blue,burly wood,olive drab,light pink,medium purple,light sky blue,bisque,rosy brown,lavender,dark golden rod'.split(',')
    print('pale colors', len(pale_colors))
    
    colors = f'{dropbox}Kira/iTol_templates/rgb_colors.txt'
    colors = open(colors).read().strip().split('\n')
    colors = [c.split('\t')[0] for c in colors[40:]]
    color_pointer = 0
    
    arg_classes = ['beta-lactam', 'aminoglycoside', 'trimethoprim', 'sulphonamide', 'phenicol', 'tetracycline', 'fosfomycin', 'macrolide', 'colistin', 'quinolone', 'rifampicin']
    oriTDB_classes = ['relaxase_all_blastx', 'auxiliary_all_blastx', 'oriT_all', 't4cp_all_blastx']

    
    for file in describe:
        file = open(file).read().strip().split('\n')
        file = [f.split('\t')[:2] for f in file if f]
        file = [f for f in file if 'ref_CDS' not in f[0]]
        for f in file:
            db, el = f
            if 'ref_CDS' not in db:
                if db in arg_classes:
                    db = 'ARG'
                elif db in oriTDB_classes:
                    db = 'oriTDB'
                elif 'backbone' in db:
                    db = el
                elif 'D6_putative_replicon' in db or 'enterobacteriaceae' in db:
                    db = 'PlasmidFinder'
                elif db in my_colors:
                    pass
                else:
                    print('LOST CLASS', db)
                
                if db not in my_colors:
                    color = pale_colors[color_pointer]
                    my_colors[db] = color
                    print(color_pointer, db, color)
                    color_pointer += 1
                    if len(pale_colors) <= color_pointer:
                        color_pointer -= len(pale_colors)
    
    lists_of_feats = {}
    if 1:        
        for plas in plasmids[:]:
            name = plas.split('/')[-1].split('.')[0]
            gb_name = name.split('_')
            if len(gb_name[0]) > 4:
                gb_name = gb_name[0]
            else:
                gb_name = '_'.join(gb_name[:2])
                                             
            host = plasmid_host[gb_name]
            
            dna = ''.join(open(plas).read().strip().split('\n')[1:])
            size = len(dna)
            
            feat_files = [a for a in describe if name in a]
            feats = []
            
            feats_back = [['', '1', '20', '+', 'black'], ['', f'{size - 20}', str(size), '+', 'black']]
            
            for file in feat_files:
                file = open(file).read().strip().split('\n')
                file = [f.split('\t')[:7] for f in file if f]
                for f in file:
                    if len(f) == 5:
                        db, el, c1, c2, nap = f
                        if 'backbone' in db:
                            db = el
                        color = my_colors[db]
                        feats_back.append(f[1:] + [color])
                        
                    elif len(f) >= 7:
                        db, el, c1, c2, nap, pident, qcovs = f
                        
                        if 'ref_CDS' not in db:
                            if db in expand_feats:
                                for element in expand_feats[db]:
                                    if element in el:
                                        if db == 'D6_putative_replicon_orf42':
                                            print(db, el, expand_feats[db][element])
                                        el = el.replace(element, expand_feats[db][element])
                            
                            if db in arg_classes:
                                db = 'ARG'
                            elif db in oriTDB_classes:
                                db = 'oriTDB'
                            elif 'backbone' in db:
                                db = el
                            elif 'D6_putative_replicon' in db or 'enterobacteriaceae' in db:
                                db = 'PlasmidFinder'
                        
                            if db not in my_colors:
                                print('LOST COLOR!!', db)
                                            
                            color = my_colors[db]
                            # if color not in ('gray', 'light gray',):
                            if 1:
                                feats.append([el + '+++' + color, int(c1), int(c2), nap, float(pident), float(qcovs)])
                    else:
                        print('LOST string', f)
            
            
            filter_feats, report = myf.find_overlaps(feats)
            # print(name)
            # print(report)
            
            feats = feats_back[::]
            feats = []
            for ff in filter_feats:
                ff = [str(f) for f in ff[:4]]
                el, color = ff[0].split('+++')
                feats.append([el] + ff[1:] + [color])
                
            feats.sort(key = lambda x: int(x[1]))
            
            lists_of_feats[name] = [host, size, feats, dna]
            print(name, len(feats))
            print('-----------------\n\n\n')
            
    #Easyfig      
    if 0 and lists_of_feats:
        for name in lists_of_feats:
    #        feats.sort(key = lambda x: int(x[1])) #do not sort so ARGs and IS elements will be in the top of the backbone blocks
            gb_file = gb_fold + name + '.gb'
            host, size, feats, dna = lists_of_feats[name]
            generate_gb(gb_file, name, host, size, feats, dna, easyfig=True)
        
        
        gb = glob.glob(gb_fold + '*')
        if tree:
            tree = open(tree).read().strip()
            order = re.findall('[A-Z][A-Z][0-9_][0-9]+', tree)
        else:
            order = [f.split('/')[-1].split('_')[0] for f in plasmids]
        print(*order, sep = '\n')
        
        query = []
        for o in order:
            query += [g for g in gb if o in g]
        print(query)
        print(len(query))
        query = ' '.join(query)
        
        figure = viz_fold + group_name + '.svg'
        
        Easyfig(query, figure, viz_fold, PC_lab)
        
        
    if 1 and lists_of_feats:
        #geneplotr
        for name in lists_of_feats:
    #        feats.sort(key = lambda x: int(x[1])) #do not sort so ARGs and IS elements will be on the top of the backbone blocks
            gb_file = gb_fold + name + '.gb'
            host, size, feats, dna = lists_of_feats[name]
            generate_gb(gb_file, name, host, size, feats, dna, easyfig=False)
            
print('-----\n\n', len(lost_hash), *lost_hash[:5], sep = '\n')