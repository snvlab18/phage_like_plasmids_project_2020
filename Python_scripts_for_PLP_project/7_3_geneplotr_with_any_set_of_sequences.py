#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 12:55:56 2020

@author: lutra
"""
import glob
import sqlite3
import subprocess


def run(cmd):
#    print(cmd)
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    return


def blast(query, subject, output):
    cmd='blastn -query %s -subject %s -out %s -outfmt 6 -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.000001' %(query, subject, output)
    run(cmd)
    return
    

PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
else:
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'
    


html_data = main + 'Charts/html_data.txt'
if html_data not in glob.glob(main + 'Charts/*'):
    with open(html_data, 'w') as fh:
        pass


#collect all files
PLP_groups = ['D6', 'P1', 'SSU5']

genbank_files = []
fasta_files = []
for PLP in PLP_groups:
    fold_rotate = f'{main}prophages_groups/{PLP}/rotate/'
    fasta_files += glob.glob(fold_rotate + '*')
    
    fold_gb = f'{main}prophages_groups/{PLP}/HomBlocks_per_cluser/Viz_GenBank/genbank/'
    genbank_files += glob.glob(fold_gb + '*')

print(len(fasta_files), len(genbank_files))

#geneplotr folders
viz_fold = main + 'Charts/genbank_show/'
blast_fold = main + 'Charts/genbank_show/genbank_blastn/'

blast_null = main + 'Charts/genbank_show/blast_null.txt'
blast_null = ''
        
#collect_species
table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'

database = f'{main}{table}.sqlite3'
conn = sqlite3.connect(database)
cur = conn.cursor()

plasmid_host = {}
improve_seq_names = {}
collect_data = []
task = 'SELECT project_ID, nucleotide, organism, completeness, slen, PLP_status, project_ID_number, Major_replicon FROM ' + table
for row in cur.execute(task):
    proj_ID, nt, org, compl, slen, plp_status, proj_ID, major_rep = [str(r) for r in row]
    org = ' '.join(org.split(' ')[:2])
    org = org.replace('virus', 'phage')
    plasmid_host[nt] = org
    if compl == 'complete' and int(slen) < 2000000 and '_' in plp_status:
        improve_seq_names[nt] = plp_status.split('_')[1] + '_' + proj_ID
        collect_data.append([nt, org, int(slen), plp_status, proj_ID, major_rep])

    
#prepare long lists of IDs
clusters = []
if 1:
    Inc_types = {'D6': ['D6_putative_replicon', 'other'], 'P1': ['IncY_1__K02380', 'p0111_1__AP010962'], 'SSU5': ['IncFIB_pHCM2_1__AL513384', 'IncFIB_pKPHS1_1__CP003223', 'IncFIB_H89-PhagePlasmid_1__HG530657', 'IncFIB_pLF82-PhagePlasmid_1__CU638872', 'other']}
    
    proj = 'D6'
    order = [c[0] for c in collect_data if c[0] in c[3] and proj in c[3]]
    upd_names = {order[0]: proj + '_1_' + order[0]}
    
    #SSU5: [3:], [2:3], [1:2] :50/ 50:, [:1] :50/ 50:;
    #P1: [:1] :50/ 50:, [1:]
    #D6: [:]
    for rep in Inc_types[proj][:]:
        clusters.append(rep)
        add = [c for c in collect_data if c[0] not in c[3] and proj in c[3] and rep == c[-1]]
        add.sort(key = lambda x: (int(x[-2])))
        order += [a[0] for a in add][:]
        
        for a in add:
            upd_names[a[0]] = f'{proj}_{a[-2]}_{a[0]}'
        # order += [a[0] for a in add if 'pestis' in a[1]][:]
        
    if 0:
        order += ['CP012254', 'CP014707', 'CP053410', 'KP763470', 'KY515226', 'CP030187', 'CP044180', 'CP057631', 'CP057633']
   
if 0:
    # proj = 'P1_mcr1_rep_quest'
    # order = ['AF234172', 'MK256965', 'KX880944', 'CP047090']
    
    # proj = 'SSU5_rearrengments'
    # order = ['JQ965645', 'CP052390', 'CP052311', 'CP052393', 'CP052352', 'CP052672', 'CP041426']
    
    # proj = 'SSU5_viruses'
    # order = ['JQ965645', 'CP054387', 'CP053388', 'MK422451']
    
    
    #D6
    proj = 'D6_replicons'
    order = ['MF356679', 'AP018804', 'CP042620']
    
    # proj = 'D6_CP055609'
    # order = ['MF356679', 'CP055609']
    
    
    #P1
    # proj = 'P1_ARG_phages'
    # order = ['AF234172', 'AF503408', 'KU760857', 'FO818745', 'MH445380']
    
    # proj = 'P1_mcr1_single'
    # order1 = ['AF234172', ] + ['KX880944', 'KX518745', 'MH781719', 'LC511659', 'MK256965', 'MG288678', 'CP047090', 'MF510496', 'MN232192', 'MH844525', 'CP035313']
    # order = [o for o in order if o in order1]
    # print(*order, sep = ' ')
    
    
    #SSU5
    # proj = 'SSU5_Ypestis'
    # order= ['JQ965645', 'AL117211', 'CP000670', 'CP045155', 'CP033694']
    # order = ['JQ965645', 'CP001587', 'AL117211', 'CP015119', 'CP000900', 'CP009714', 'CP006807']
        
    # proj = 'SSU5_various_replicons'
    # order = ['JQ965645', 'CP040567', 'CP056704', 'CP050499', 'CP029421']
    # others = ['CP057631', 'CP057633', 'CP012254', 'CP014707', 'KY515226', 'KP763470', 'CP053410', 'CP030187', 'CP044180']
    # order += others
    
    
    
    
# order = order[:4]
     
print(len(order), order[:5])
fut_name = f'{proj}_{"_".join(clusters)}_n{len(order)}'
print(fut_name)


#input for Flexidot
if 0:
    flexidot_fold = main + 'Charts/flexidot/'
    gepard_fold = main + 'Charts/gepard/'
    gepard_file = gepard_fold + fut_name + '.fasta'
    gepard = open(gepard_file, 'w')
    
    save = []
    for o in order[:]:
        file = [f for f in fasta_files if o in f][0]
        upd_o_name = upd_names[o]
        upd_file = f'{flexidot_fold}fasta/{upd_o_name}.fasta'
        with open(upd_file, 'w') as fh:
            prev_file = open(file).read().strip().split('\n')
            new_file = '>' + upd_o_name + '\n' + '\n'.join(prev_file[1:]) + '\n'
            fh.write(new_file)
            gepard.write(new_file)
        save.append(upd_file)
        
    with open(f'{flexidot_fold}flexidot_in_{fut_name}.txt', 'w') as fh:
        fh.write(','.join(save) + '\n')
        
    gepard.close()

#geneplotr
if 1:    
    super_fasta = ''
    if len(order) < 10:
        fut_name = '_'.join([proj] + order)
        print(order)
    else:
        fut_name = f'{proj}_{"_".join(clusters)}_n{len(order)}'
        print(fut_name)
        # super_fasta = viz_fold + fut_name + '.fasta'
        # with open(super_fasta, 'w') as fh:
        #     pass
    
    gb_load = []
    gb_fill = []
    gb_list = []
    offsets = []
    
    names_separately = []
    names = []
    blast_comp = []
    blast_list = []
    
    
    count = 1
    gbl = 'gb{} <- read_dna_seg_from_file("{}"))'
    gbf = 'gb{}$fill <- gb{}$col'
    
    collect_anchors = []
    for o in order:
        new_o = improve_seq_names[o] + ' ' +  o + ' ' + plasmid_host[o]
        collect_anchors += [improve_seq_names[o], o]
        
        gb_name = f'gb{count}'
        if len(order) < 50:
            names.append(f'"{new_o}"')
        else:
            names_separately.append(f'name{count} <- c("{new_o}")')
            names.append(f'name{count}')
        
        o_file = [g for g in genbank_files if o in g][0]
        # print(o_file)
        
        gbl = f'{gb_name} <- read_dna_seg_from_file("{o_file}")'
        gb_load.append(gbl)
        
        gbf = f'{gb_name}$fill <- {gb_name}$col'
        gb_fill.append(gbf)
        
        gb_list.append(gb_name)
        offsets.append('0')
        
        if super_fasta:
            add_fasta = [f for f in fasta_files if o in f][0]
            with open(super_fasta, 'a') as fh:
                fh.write(open(add_fasta).read())
            
        
        if o != order[-1]:
            query = [f for f in fasta_files if order[count-1] in f][0]
            subject = [f for f in fasta_files if order[count] in f][0]
            output = blast_fold + order[count-1] + '_vs_' + order[count] + '.blast'
            
            if False or output not in glob.glob(blast_fold + '*') or 'CP056516' in fut_name:
                blast(query, subject, output)
                output_proc = open(output).read().strip().split('\n')
                output_proc = [o.split('\t') for o in output_proc] 
                output_proc = ['\t'.join(o) for o in output_proc if int(o[3]) > 500]
                
                with open(output, 'w') as fh:
                    fh.write('\n'.join(output_proc) + '\n')
                print(output)
            
            if blast_null:
                output = blast_null
            
            blast_c = f'blast{count} <- read_comparison_from_blast("{output}")'
            blast_comp.append(blast_c)
            blast_list.append(f'blast{count}')
        count += 1
        
    
    gb_load = '\n'.join(gb_load) + '\n\n'
    gb_fill = '\n'.join(gb_fill) + '\n\n'
    
    gb_list = ", ".join(gb_list)
    gb_list = f'dna_seqs <- list({gb_list})\n\n'
    
    
    names = ", ".join(names)
    names = f'dna_seg_labels <- c({names})\n\n'
    if names_separately:
        names_separately = '\n'.join(names_separately)
        names = names_separately + '\n\n' + names
    
    offsets = ', '.join(offsets)
    offsets= f'offsets <- list({offsets})\n\n'
    
    rot_vars = {'P1': '90', 'D6': '90', 'SSU5': '90'} #D6 was 45
    
    annots = '''annots <- lapply(dna_seqs, function(x){
      mid <- middle(x)
      annot <- annotation(x1=mid, text=x$name, rot=90)\n})\n\n'''
    
    blast_comp = '\n'.join(blast_comp) + '\n\n'
    
    blast_list = ", ".join(blast_list)
    blast_list = f'blast <- list({blast_list})\n\n'
    
        
    
    start = f'''##
    ## Figure of {fut_name}
    ##
    
    #install.packages("ade4")
    #install.packages("grid")
    #install.packages("genoPlotR")
    library(ade4)
    library(grid)
    library(genoPlotR)
    '''
    
    width = 20
    if len(order) > 20:
        width = 50
    height = len(order)*1.15
    height = round(height, 2)

    annotation_cex=1
    
        
    main_name = f'{fut_name}'
    if blast_null:
        main_name += '_null'
        
    print(main_name)
    
    plot = f'''plot_gene_map(dna_segs=dna_seqs,
                  main="{main_name}", main_pos="right",
                  gene_type="side_blocks",
                  comparisons=blast, offsets=offsets,
                  global_color_scheme = c("auto", "auto", "blue_red", 0.2), override_color_scheme=TRUE,
                  annotations=annots, annotation_cex={annotation_cex}, annotation_height=0,
                  dna_seg_labels=dna_seg_labels,
                  dna_seg_label_cex=1.1,
                  dna_seg_scale=TRUE, scale=TRUE)'''
    
    save_image = f'''svg(file.path("{viz_fold}", "{main_name}_geneplotr.svg"),
          width={width}, height={height})
    
    {plot[:-1]},
                  plot_new=FALSE)
    
    dev.off()
    
    '''
    
    blocks = [start, gb_load + gb_fill, names, gb_list + offsets + annots, blast_comp + blast_list, save_image]
    blocks = [b for b in blocks if b]
    if len(order) < 150:
        blocks = ['\n\n'.join(blocks)]
    
    blocks = ['```{r}\n' + b + '\n\n```\n\n' for b in blocks]
    r = ''.join(blocks)
    
    figure = viz_fold + main_name + '.rmd'
    
    print(figure)
    with open(figure, 'w') as fh:
        fh.write(r + '\n')
    
    print('--\n')
    
    if main_name not in open(html_data).read():
        with open(html_data, 'a') as fh:
            collect_anchors = [main_name] + collect_anchors
            fh.write('\t'.join(collect_anchors) + '\n')
    
    