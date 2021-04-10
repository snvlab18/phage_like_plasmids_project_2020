#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:19:05 2020

@author: lutra
"""
import glob
import subprocess
import sqlite3


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
#    else:
#        print(fold, 'exists')
    return


def species_name(txt):
    txt = txt.split()
#    if txt[1] not in ('virus', 'phage', 'sp.'):
#        return f'{txt[0][:1]}.{txt[1]}'
    return f'{txt[0]}_{txt[1]}'


def rev_comp(dna):
    trantab = "".maketrans("ACTGactg","TGACtgac")
    return dna[::-1].translate(trantab)


def dna_to_fasta(dna):
    dna = [dna[i*60:(i+1)*60] for i in range(int(len(dna)/60)+1)]
    dna=[d for d in dna if d]
    fasta = '\n'.join(dna)
    return fasta
    

def rotate_around_coord(seq, c1, rev = False):
    '''Specified coord will become first in the new seq.
    If not rev dir c1 is included into transfered part (I move the match with its first coord), 
    else excluded (I move evth after the match)'''
    if not rev:
        new = seq[c1-1:] + seq[:c1-1]
    else:
        new = seq[c1:] + seq[:c1]
        new = rev_comp(new)
    
    check = new + new
    if flag:
        check = rev_comp(check)
    if seq not in check:
        return 'ERROR!!'
    
    return new
    

def iTol_color_leaves(filename, title, dataset_color, dataset_shape, data, color_db=False, PC='shiri'):
    if '/data/' in filename:
        PC = 'lutra'
    
    sp_colors = '''
    Citrobacter_freundii (255,140,0)
    Citrobacter_sp. (255,140,0)
    Cronobacter_sakazakii (128,0,128)
    Enterobacter_cloacae (184,134,11)
    Enterobacter_hormaechei (184,134,11)
    Enterobacteria_phage (245,222,179)
    Escherichia_albertii (255,0,0)
    Escherichia_coli (255,0,0)
    Escherichia_fergusonii (255,0,0)
    Escherichia_marmotae (255,0,0)
    Escherichia_phage (250,128,114)
    Escherichia_virus (250,128,114)
    Klebsiella_aerogenes (0,0,255)
    Klebsiella_grimontii (0,0,255)
    Klebsiella_oxytoca (0,0,255)
    Klebsiella_phage (65,105,225)
    Klebsiella_pneumoniae (0,0,255)
    Klebsiella_variicola (0,0,255)
    Kluyvera_cryocrescens (139,69,19)
    Salmonella_enterica (0,128,0)
    Salmonella_phage (50,205,50)
    Salmonella_sp. (0,128,0)
    Shigella_sonnei (255,0,0)
    Yersinia_pestis (255,255,0)
    uncultured_bacterium (220,20,60)
    '''
    
    print(filename, PC)
    if not color_db:
        sp_colors = sp_colors.strip().split('\n')
        sp_colors = [s.strip().replace(')', ',0.5)').replace('(', 'rgba(').split( ) for s in sp_colors]
        color_db = dict(sp_colors)
        print(color_db)
        
    template = f'/home/{PC}/Dropbox/Kira/iTol_templates/dataset_styles_template.txt'
    tmp = open(template).read().strip()
    
    tmp = tmp.replace('placeholder_title', title)
    tmp = tmp.replace('placeholder_dataset_color', dataset_color)
    
    diversity = list(set([d[1] for d in data]))
    diversity = [d.replace(' ', '_') for d in diversity]
    diversity.sort()
    colors = [color_db[c] for c in diversity]
    
    shapes = ' '.join([dataset_shape] * len(diversity))
    tmp = tmp.replace('placeholder_shapes', shapes)
    tmp = tmp.replace('placeholder_colors', ' '.join(colors))
    tmp = tmp.replace('placeholder_labels', ' '.join(diversity))
    tmp += '\n'
    
    for d in data:
        d[1] = d[1].replace(' ', '_')
        entry = [d[0], 'label', 'node', 'rgb(0,0,0)', '1', 'normal', color_db[d[1]]]
        tmp += ' '.join(entry) + '\n'
    
    with open(filename, 'w') as fh:
        fh.write(tmp)
    print(title, 'is ready!\n')
    return
    
    
def iTol_strip(filename, title, dataset_color, dataset_shape, data, reverse_colors = False, rotate_colors = 0, PC='shiri'):
    if '/data/' in filename:
        PC = 'lutra'
    
    colors0 = '#0000FF #00FF00 #FF0000 #FF00FF #2F4F4F #800000 #9b9bff #ffc080 #ff7373 #008080 #FF8C00 #808000 #8FBC8F #00BFFF #8A2BE2 #EE82EE #FFE4C4 #D2691E #778899'.split()
    if reverse_colors:
        colors0 = colors0[::-1]
    if rotate_colors:
        colors0 = colors0[rotate_colors:] + colors0[:rotate_colors]
        print('rotated for ', rotate_colors)
    
    template = f'/home/{PC}/Dropbox/Kira/iTol_templates/iTol_strip_template.txt'
    tmp = open(template).read().strip()
    
    tmp = tmp.replace('placeholder_title', title)
    tmp = tmp.replace('placeholder_dataset_color', dataset_color)
    
    diversity = list(set([d[1] for d in data]))
    diversity = [d.replace(' ', '_') for d in diversity]
    diversity.sort()
    
    shapes = ' '.join([dataset_shape] * len(diversity))
    colors = colors0[:len(diversity)]
    
    print(len(diversity), len(colors))
    
    div_colors = dict(zip(diversity, colors))
    
    tmp = tmp.replace('placeholder_shapes', shapes)
    tmp = tmp.replace('placeholder_colors', ' '.join(colors))
    tmp = tmp.replace('placeholder_labels', ' '.join(diversity))
    tmp += '\n'
    
    for d in data:
        d[1] = d[1].replace(' ', '_')
        entry = [d[0], div_colors[d[1]], d[1]]
        tmp += ' '.join(entry) + '\n'
    
    with open(filename, 'w') as fh:
        fh.write(tmp)
    print(title, 'is ready!\n')
    return len(diversity)


def iTol_binary(filename, title, dataset_color, dataset_shape, header, data, PC='shiri'):
    if '/data/' in filename:
        PC = 'lutra'
    
    colors0 = '#0000FF #00FF00 #FF0000 #FF00FF #2F4F4F #800000 #9b9bff #ffc080 #ff7373'.split()
    shapes0 = [str(i) for i in range(1, 7)]
    template = f'/home/{PC}/Dropbox/Kira/iTol_templates/iTol_binary_template.txt'
    tmp = open(template).read().strip()
    
    tmp = tmp.replace('placeholder_title', title)
    tmp = tmp.replace('placeholder_dataset_color', dataset_color)
    
    colors = []
    shapes = []
    while len(colors) < len(header):
        colors = colors[::] + colors0[::]
    while len(shapes) < len(header):
        shapes = shapes[::] + shapes0[::]
    colors = colors[:len(header)]
    shapes = shapes[:len(header)]
    
    tmp = tmp.replace('placeholder_shapes', ','.join(shapes))
    tmp = tmp.replace('placeholder_colors', ','.join(colors))
    tmp = tmp.replace('placeholder_lables', ','.join(header))
    
    data = [','.join(d) for d in data]
    data = '\n' + '\n'.join(data)
    tmp += data
    
    with open(filename, 'w') as fh:
        fh.write(tmp)
    print(title, 'is ready!\n')
    return


PC_lab = True

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
else:
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'


copy = False
rotate = False

fasta_files = main + 'fasta/'
coords_files = main + 'coverage_per_nt_results/coverage_per_nt_run/'

fasta_files = glob.glob(fasta_files + '*')
coords_files_check = glob.glob(coords_files + '*')

table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'
database = main + table + '.sqlite3'
conn = sqlite3.connect(database)
cur = conn.cursor()

Inc_types = {'D6': [], 'P1': ['p0111_1__AP010962', 'IncY_1__K02380'], 'SSU5': ['IncFIB_pKPHS1_1__CP003223', 'IncFIB_pLF82-PhagePlasmid_1__CU638872', 'IncFIB_H89-PhagePlasmid_1__HG530657', 'IncFIB_pHCM2_1__AL513384']}

project = ''

mkdir(main + 'prophages_groups')


#second switch
extra_Incs_count = 0
ARG_extra_Incs_count = 0
for phage in ['D6', 'P1', 'SSU5'][-1:]:
    print(phage)
    
    phage_fold = main + 'prophages_groups/' + phage + '/'
    mkdir(phage_fold[:-1])
    
    select = []
    
    simple_types = []
    
    ARG_plus_extra_Incs = []
    no_ARGs_extra_Incs = []
    
    record = []
    task = 'SELECT nucleotide, organism, completeness, genome, slen, shape, PLP_status, enterobacteriaceae, ResFinder FROM '+ table
    save_len = {}
    for row in cur.execute(task):
        gid, host, compl, genome, slen, shape, plp, plasmidfinder, resfinder = [str(r) for r in row]
        if compl == 'complete' and '_' + phage in plp:
            source = genome
            if source != 'plasmid':
                if source == 'chromosome' or int(slen) > 1000000:
                    source = 'chromosome'
                else:
                    source = 'phage'
            
            save_len[gid] = int(slen)
            species = species_name(host)
            
            inc_type = []
            for inc in Inc_types[phage]:
                if inc in plasmidfinder:
                    inc_type.append(inc)
            
            if not inc_type:
                inc_type = ['other']
#                print(gid, host)
                    
            inc_type = '_'.join(inc_type)
            
            cool_ARGs = ['NDM', 'KPC', 'IMP', 'mcr']
            check_cool_ARGs = [c for c in cool_ARGs if c in resfinder]
            
            if resfinder == '0':
                res = '-'
            else:
                resfinder_N = set([r.split(' ')[0] for r in resfinder.split('\n')])
                if len(resfinder_N) >= 3:
                    res = 'MDR'
#                    print(gid, *resfinder)
                else:
                    res = '+'
            
            extra_incs = []
            if plasmidfinder != 0:
                extra_incs = [p.split()[0].replace('#', '').replace('**', '') for p in plasmidfinder.split('\n')]
                extra_incs = [e for e in extra_incs if e not in Inc_types[phage] and e != '0']
                
            if source != 'chromosome' and extra_incs:
                extra_Incs_count += 1
            
            if source != 'chromosome' and extra_incs and resfinder != '0':
                ARG_extra_Incs_count += 1
                
            if source != 'chromosome' and extra_incs and resfinder == '0':
                no_ARGs_extra_Incs.append([gid] + extra_incs)
            
            
            # if (source != 'chromosome' and extra_incs and resfinder != '0') or (gid in plp) or (source != 'chromosome' and check_cool_ARGs):
            #     arg_cluster_name = 'ARG_extra' 
            if (source != 'chromosome' and extra_incs and resfinder != '0') or (gid in plp):
                arg_cluster_name = 'ARG_extra_Incs'
            # if (source != 'chromosome' and 'mcr' in check_cool_ARGs and not extra_incs) or (gid in plp):
            #     arg_cluster_name = 'mcr'
                
                # print(gid + ' ' + source, extra_incs, resfinder, '----\n', sep = '\n')
                add_order = 2
                if gid in plp:
                    add_order = 0
                elif extra_incs and resfinder != '0':
                    add_order = 1
                ARG_plus_extra_Incs.append([gid, add_order, inc_type, host, slen])
            
            record.append([gid, species, res, source, inc_type, shape])
            
    ARG_plus_extra_Incs.sort(key = lambda x: (x[2], x[1], x[3], int(x[4]), x[0]))
    print('here')
    print(*ARG_plus_extra_Incs, sep = '\n')
    print(len(ARG_plus_extra_Incs))
    ARG_plus_extra_Incs = [a[0] for a in ARG_plus_extra_Incs]
    
    print('\n* * *')
    print('no ARG extra Incs', len(no_ARGs_extra_Incs), *no_ARGs_extra_Incs, '* * * *\n', sep = '\n')
    
            
    print('N', len(record))
    
    select = [r[0] for r in record if r[3] != 'chromosome']
    fasta = [g for g in fasta_files if g.split('/')[-1].split('.')[0] in select]
    
    select_prophages = [r[0] for r in record if r[3] == 'phage']
    fasta_prophages = [g for g in fasta_files if g.split('/')[-1].split('.')[0] in select_prophages]
    
    phage_fasta = phage_fold + 'fasta/'
    mkdir(phage_fasta[:-1])
    
    rotate_fasta = phage_fold + 'rotate/'
    mkdir(rotate_fasta[:-1])
    
    
    group_by_Rep = phage_fold + 'IDs_grouped_by_Rep.txt'
    with open(group_by_Rep, 'w') as fh:
        unique_Incs = list(set([r[4] for r in record]))
        unique_Incs.sort()
        print(*unique_Incs)
        total_count = 0
        for uniq in unique_Incs:
            IDs = [r[0] for r in record if r[4] == uniq and r[3] != 'chromosome']
            total_count += len(IDs)
            fh.write(uniq + ' ' + ' '.join(IDs) + '\n')
        
        fh.write(arg_cluster_name + ' ' + ' '.join(ARG_plus_extra_Incs) + '\n')
        ypestis = [r[0] for r in record if 'pestis' in r[1] and r[3] != 'chromosome']
        print(len(ypestis))
        ypestis = ['Ypestis'] + ypestis
        fh.write(' '.join(ypestis) + '\n')
        print('total count', total_count)
        
    print('--\n')
    
    
    
    for f in fasta[:]:
        name = f.split('/')[-1].split('.')
        new_name = name[0] + '.' + name[-1]
        if copy:
            run(f'cp {f} {phage_fasta}{new_name}')
        
        if rotate: 
            coords = f'{coords_files}{name[0]}_{phage}.txt'
            if coords not in coords_files_check:
                print(name[0], 'lost coords!', coords)
            else:
                seq = ''.join(open(f).read().strip().split('\n')[1:])
                slen = save_len[name[0]]
                if len(seq) != slen:
                    print(name[0], 'seq length problem!', len(seq), slen)
                coords = open(coords).read().strip().split('\n')
                coords = [c.split('\t') for c in coords]
                coords = [c for c in coords if int(c[6]) > 100]
                coords.sort(key = lambda x: int(x[9]))
                first = coords[0]
                
                if int(first[9]) >= 100:
                    print(name, 'start problem!', *first)
                    print(*[c[9:11] for c in coords[0:5]], sep = '\n')
                    print('')
                    
    #            print(*first)
                new = ''
                flag = ''
                c1, c2 = [int(f) for f in first[11:13]]
                if c2 < c1:
                    flag = 'rev'
                    
                new = rotate_around_coord(seq, c1, rev=flag)
                if len(new) != slen:
                    print(name[0], 'new length problem!', len(new), slen)
                    
                with open(f'{rotate_fasta}{name[0]}_{phage}_c{c1}{flag}.fasta', 'w') as fh:
                    fh.write(f'>{name[0]}_{phage}_c{c1}{flag}\n')
                    fh.write(dna_to_fasta(new) + '\n')
    
    
    prophage_fold = phage_fold + 'prophages_only/'
    mkdir(prophage_fold[:-1])
    
    print('Here are prophages')
    ready_rotate = glob.glob(rotate_fasta + '*')
    for prophage in select_prophages:
        print(prophage)
        ready_fasta = [r for r in ready_rotate if prophage in r][0]
        cmd = f'cp {ready_fasta} {prophage_fold}'
        run(cmd)
    
    if False:
        #create iTol docs
        mkdir(phage_fold + 'iTol_docs')
        
        main_color = '#0000FF #00FF00 #FF0000 #FF00FF #2F4F4F #800000 #9b9bff #ffc080 #ff7373'.split()
        main_shape = [str(i) for i in range(1, 7)]
        file = f'{phage_fold}iTol_docs/host_iTol.txt'
        data = [[r[0], r[1]] for r in record]
        print(data[:5])
        iTol_color_leaves(file, 'host', main_color[0], main_shape[0], data)
        
        pointer1 = 1
        pointer2 = 0
        
        mkdir(phage_fold + 'iTol_docs')
        
        order = ['Res_status', 'DNA_source', 'Inc_class', 'Shape']
        for o in order[:]:
    #        iTol_strip(host_file, 'plasmid_host', '#ff7373', '1', hosts, True)
            file = f'{phage_fold}iTol_docs/{o}_iTol.txt'
            data = [[r[0], r[pointer1+1]] for r in record]
            print(data[:5])
            to_rotate = iTol_strip(file, o, main_color[pointer1], main_shape[pointer1], data, rotate_colors = pointer2)
            pointer1 += 1
            pointer2 += to_rotate
            
        
        print(len(record))
        print(len([r for r in record if r[4] != '-']))
        record.sort(key = lambda x: x[1])
        record = ['\t'.join(r) for r in record]
        with open(phage_fold + 'seq_description.csv', 'w') as fh:
            fh.write('\t'.join(['nucleotide', 'DNA_source', 'host', 'Inc_class', 'Res_status']) + '\n')
            fh.write('\n'.join(record) + '\n')
            
    
    print('--\n\n\n')

print('all_extra_Incs_count', extra_Incs_count)
print('ARG_extra_Incs_count', ARG_extra_Incs_count)