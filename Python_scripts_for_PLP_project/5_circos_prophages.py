#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 12:45:12 2020

@author: lutra
"""
import re
import sqlite3
import pandas as pd
from PIL import ImageColor


def species_name(txt):
    txt = txt.split()
    # if txt[1] != 'sp.':
    #     return f'{txt[0][:1]}. {txt[1]}'
    # else:
    #     return f'{txt[0]} {txt[1]}'
    
    txt = ' '.join(txt[:2])
    txt = txt.replace('virus', 'phage')
    return txt
    

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


def PlasmidFinder(txt, txt1):
    if txt == '0':
        txt = ''
    if txt1 == '0':
        txt1 = ''
    
    txt = (txt + '\n' + txt1)
    if txt == '\n':
        txt = '0'
    txt = txt.strip().split('\n')
    
    plas = [t.strip().split()[0].replace('D6_putative_replicon_orf42_c45699_46596', 'D6orf42').split('_')[0].split('(')[0].replace('#', '').replace('*', '').replace('/','') for t in txt]
    
    if plas == ['0']:
        plas = ['UD']
    plas = list(set(plas))
    
    order = [o for o in ('D6_orf42', 'IncY', 'p0111', 'IncFIB', 'IncFIA', 'IncFII', 'IncHIA', 'IncHI1A', 'IncHI1B') if o in plas]
    plas = [p for p in plas if p not in order]
    plas.sort()
    plas = order + plas
    
    
    # if len(txt) > 1:
    #     print(txt, *plas, '\n')
    return '_'.join(plas)


def PLP_group(SSU5, P1, D6):
#    SSU5, P1, D6 = float(SSU5). float(P1), float(D6)
    thresh = 40
    if SSU5 >= thresh:
        return 'SSU5'
    elif P1 >= thresh:
        return 'P1'
    elif D6 >= thresh:
        return 'D6'
    else:
        return 'lost'


PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
    dropbox = '/home/shiri/Dropbox/'
else:
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'
    dropbox = '/home/lutra/Dropbox/'
    
    
if 0:
    table = 'Phage_like_plasmids_SSU5_P1_D6_8Sep20'
    database = main + table + '.sqlite3'
    print(database)
    
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    table = pd.read_sql('SELECT * FROM ' + table, conn)
    print(list(table))
    print(table.shape)
    df = table[(table.completeness == 'complete') & (table.slen < 2000000)]
    print(df.shape)
    
    df = df.fillna('-')
    
    df['PLP_group'] = df.apply(lambda row: PLP_group(row.JQ965645_SSU5_ref_cov, row.AF234172_P1_ref_cov, row.MF356679_D6_ref_cov), axis = 1)
    df['sp_names'] = df.apply(lambda row: species_name(row.organism), axis = 1)
    df['plas_list'] = df.apply(lambda row: PlasmidFinder(str(row.enterobacteriaceae), str(row.D6_putative_replicon_orf42)), axis = 1)
    df['arg_list'] = df.apply(lambda row: ResFinder(str(row.ResFinder)), axis = 1)
    
    print(df)
    
    group = df.groupby(["PLP_group", "sp_names", "plas_list", 'arg_list'], as_index=False)["nucleotide"].count()
    print(group)
    
    plas_colors = list(group.plas_list.unique())
    plas_colors.sort(key = lambda x: x.upper())
    print(*plas_colors, len(plas_colors), sep = ',')
    
    group.to_csv(main + 'Circos/Circos_group.csv', sep = '\t', index = False)



if 0:
    plasmid_colors = ['IncA/C2', 'IncFIB', 'IncFII', 'IncQ1', 'IncX1', 'IncY', 'p0111' , 'UD']
    colors = ['0,255,255', '0,0,204', '255,0,255', '255,255,0', '153,0,204', '28,145,236' , '238,43,49', '255,102,0']
    
    plasmid_colors = 'UD,p0111_IncFII,p0111_IncFII_IncN,IncY_IncX1,IncY_IncHI2_IncHI2A_IncN,D6orf42_IncFII,IncY_IncFIA_IncFII,D6orf42_IncFIB_IncFIA_IncFII,IncFIB_IncFII,p0111,IncFIB_IncN,p0111_IncQ1,IncFIB,IncFIB_IncC,D6orf42,p0111_IncFIB,IncY,IncY_IncFIB_IncFIA_IncFII,IncFIB_IncHI1B'
    plasmid_colors = plasmid_colors.split(',')
    
    
    colors = dropbox + 'Kira/iTol_templates/glasbey_colors_n19.txt'
    colors = open(colors).read().strip().split('\n')
    colors = [c.split() for c in colors if c[:1] == '#']
    plasmid_colors = [c[1] for c in colors]
    colors = [c[0] for c in colors]
    colors = [str(ImageColor.getrgb(c)).replace('(', '').replace(')', '').replace(' ', '') for c in colors]
    
    print(*colors, sep = '\n')
    
    
    if 1:
        plasmid_colors = dict(zip(plasmid_colors, colors))
        print(plasmid_colors['UD'])
        
        work_file = main + 'Circos/Circos_group.csv'
        for prophage in ('SSU5', 'P1', 'D6'):
            work = open(work_file).read().strip().split('\n')
            work = [w.split('\t') for w in work]
            work = [w[1:] for w in work if w[0] == prophage]
    #        work = [[w[0]] + [w[1].split('_')[0].split('(')[0]] + w[2:] for w in work]
            
            bacteria = list(set([w[0] for w in work]))
            plasmids = list(set([w[1] for w in work]))
            resistance = list(set([w[2] for w in work]))
            
            bacteria.sort()
            plasmids.sort()
            resistance.sort()
            
            if '-' in resistance:
                resistance = ['-'] + [r for r in resistance if r != '-']
            
            print(prophage, *plasmids, '', sep = '\n')
        #    print(*bacteria, '\n')
        #    print(*plasmids, '\n')
            # print(*resistance, '\n')
            
            table = []
            row_number = len(resistance) + 1 + 2
            first = ['data'] * 4 + [str(i) for i in range(1, row_number)]
            second = ['data'] * 4 + ['bln1', 'No_ARGs'] + [r.replace('-', '_').replace(',', '_') for r in resistance[1:]] + ['bln2']
            third = [row_number, 5, '0,0,0', 'bln3'] + ['-'] * (len(first)-5) + [1]
            table += [first, second, third]
            
            collect = {}
            
            for b in bacteria:
                collect[b] = {}
                for p in plasmids:
                    for r in resistance:
                        coll = [w for w in work if [b, p, r] == w[:-1]]
                        if coll:
                            if p not in collect[b]:
                                collect[b][p] = {}
                            collect[b][p][r] = sum([int(c[3]) for c in coll])
            
            row_number +=1
            for b in bacteria:
                for p in plasmids:
                    if p in collect[b]:
                        row = [row_number, 0, plasmid_colors[p], b.replace(' ','_').replace('._','.') + '_' + p]
                        row.append('-')
                        row_width = 0
                        for r in resistance:
                            if r in collect[b][p]:
                                row.append(collect[b][p][r])
                                row_width += int(collect[b][p][r])
                            else:
                                row.append('-')
                        row.append('-')
                        row[1] = row_width
                        table.append(row)
                        row_number += 1
            
            last = [row_number, 5, '0,0,0', 'bln4', 1] + ['-'] * (len(first)-5)
            table.append(last)
            
            table = ['\t'.join([str(tt) for tt in t]) for t in table]
            table = '\n'.join(table)
            table = table.replace('/', '').replace('(', '').replace(')', '').replace("'", '').replace(' ', '').replace('.','_').replace('__', '_')
            table = table.replace('_', 'ZzZ').replace('#', 'UuU').replace('**', 'YyY')
            
            with open(main + 'Circos/table_' + prophage + '.txt', 'w') as fh:
                fh.write(table)
                
if 0:
    for prophage in ('SSU5', 'P1', 'D6'):
        svg = main + 'Circos/svg/' + prophage + '.svg'
        svg_upd = main + 'Circos/svg/' + prophage + '_upd.svg'
        plain = open(svg).read()
        plain = plain.replace('ZzZ', '_').replace('UuU','#').replace('YyY', '**')
        font_size = list(set(re.findall('font-size="([0-9.]+)px"', plain)))
        font_size = max([float(f) for f in font_size])
        
        curr_font = f'font-size="{font_size}px"'
        new_font = f'font-size="{font_size*1.2}px"'
        
        
        print(prophage, font_size, new_font)
        plain = plain.replace(curr_font, new_font)
        
        with open(svg_upd, 'w') as fh:
            fh.write(plain)

#circos legend            
if 0:
    import matplotlib.pylab as plt
    import matplotlib.patches as patches
    
    colors = dropbox + 'Kira/iTol_templates/glasbey_colors_n19.txt'
    colors = open(colors).read().strip().split('\n')
    colors = [c.split() for c in colors if c[:1] == '#']
    plasmid_colors = [c[1] for c in colors]
    
    colors = [c[0] for c in colors]
    # colors = [np.array(ImageColor.getrgb(c)).reshape(-1,3) for c in colors]
    
    # Create a color palette
    palette = dict(zip(plasmid_colors, colors))
    print(plasmid_colors)
    
    handles = [patches.Patch(color=palette[x], label=x) for x in palette.keys()]
    
    # Create legend
    plt.legend(handles=handles, frameon = False)
    # Get current axes object and turn off axis
    plt.gca().set_axis_off()
    
    legend_file = main + 'Circos/svg/legend.svg'
    # plt.savefig(legend_file, format = 'svg', dpi = 600)
    plt.show()
    
    
#geneplotr legend
if 1:
    import matplotlib.pylab as plt
    import matplotlib.patches as patches
    
    colors = dropbox + 'Kira/iTol_templates/rgb_colors.txt'
    colors = open(colors).read().strip().split('\n')
    colors = [c.split('\t') for c in colors]
    
    my_colors = {'PlasmidFinder':'blue', 'ISel_db' : 'orange', 'ARG' : 'red', 'oriTDB': 'light blue', 'VFDB_setB_nt': 'dark green', 'CDS': 'gray', 'hypothetical': 'light gray', 'tRNA': 'medium orchid', 'toxin': 'lawn green', 'tail': 'sienna', 'phage': 'dark golden rod'}
    palette = {}
    my_sort = ['ARG', 'PlasmidFinder', 'VFDB_setB_nt', 'ISel_db','tRNA', 'toxin', 'oriTDB',  'phage', 'tail', 'CDS', 'hypothetical']
    mc_sort_name = ['ARG', 'Replication region', 'Virulence factor',  'IS element', 'tRNA', 'Toxin/ antitoxin', 'Conjugation gene', 'Phage protein', 'Tail protein', 'CDS', 'Hypothetical protein']
    for mc in my_sort:
        print(mc, my_colors[mc])
        rgba = [c[1] for c in colors if c[0] == my_colors[mc]][0]
        # rgba = tuple([int(r) for r in re.findall('[0-9]+', rgba)])
        palette[mc_sort_name[my_sort.index(mc)]] = rgba
    #     my_sort.append(mc)
    # print(my_sort)
    
    handles = [patches.Patch(color=palette[x], label=x) for x in palette.keys()]
    
    plt.legend(handles=handles, frameon = False)
    # Get current axes object and turn off axis
    plt.gca().set_axis_off()
    
    legend_file = main + 'Charts/geneplotr_legend.svg'
    plt.savefig(legend_file, format = 'svg', dpi = 600)
    plt.show()
    
