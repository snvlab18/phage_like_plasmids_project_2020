#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 12:55:56 2020

@author: lutra
"""
import glob
import subprocess
import sqlite3

def run(cmd):
#    print(cmd)
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    return


def blastn_cov(ref, plasmid, outp):
    cmd = f'blastn -task blastn -query {ref} -subject {plasmid} -out {outp} -outfmt "6 qaccver saccver qlen slen qcovs pident length mismatch gapopen qstart qend sstart send evalue bitscore"'
    run(cmd)
    return


def rev_comp(dna):
    trantab = "".maketrans("ACTGactg","TGACtgac")
    return dna[::-1].translate(trantab)


def dna_to_fasta(dna):
    dna = [dna[i*60:(i+1)*60] for i in range(int(len(dna)/60)+1)]
    dna=[d for d in dna if d]
    fasta = '\n'.join(dna)
    return fasta


def blast(query, subject, output):
    cmd='blastn -query %s -subject %s -out %s -outfmt 6 -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.000001' %(query, subject, output)
    run(cmd)
    return
    

def rotate_around_coord(seq, c1, flag, rev = False):
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



main0 = '/data/Current_work/Phage_like_plasmids/PLP_final/'


main = main0 + 'ARG_maps/cut_ARG_extra_Incs/'

inserts = main + 'ready/'
inserts = glob.glob(inserts + '*')

to_process = main + 'ARG_extra_Incs_matched_plasmids_turn_to_fasta/'
to_process = glob.glob(to_process + '*')

ref_seqs = main + 'ref_phages/'
ref_seqs = glob.glob(ref_seqs + '*')

PLP_groups = ['D6', 'P1', 'SSU5']

genbank_files = []
PLP_rotate = []
for PLP in PLP_groups:
    fold_rotate = f'/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/{PLP}/rotate/'
    PLP_rotate += glob.glob(fold_rotate + '*')
    
    fold_gb = f'/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/{PLP}/HomBlocks_per_cluser/Viz_GenBank/genbank/'
    genbank_files += glob.glob(fold_gb + '*')
genbank_files += glob.glob('/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/D6/HomBlocks_per_cluser/Viz_extra_matches/genbank/*')

blastn = main + 'local_blastn_ARG_extra_Incs/'
outfold = main + 'rotated_ARG_extra_Incs_matched_plasmids/'

#geneplotr
viz_fold = main0 + 'prophages_groups/D6/HomBlocks_per_cluser/Viz_extra_matches/'
# viz_fold = main0 + 'Charts/genbank_show/'
blast_fold = main0 + 'prophages_groups/D6/HomBlocks_per_cluser/Viz_extra_matches/blast_genoplotr/'
        


#improve_names
table = 'Phage_like_plasmids_SSU5_P1_D6_8Sep20'
table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'

database = f'{main0}{table}.sqlite3'
conn = sqlite3.connect(database)
cur = conn.cursor()

plasmid_host = {}
improve_seq_names = {}
collect_data = []
task = 'SELECT project_ID, nucleotide, organism, completeness, slen, PLP_status, Enterobacteriaceae FROM ' + table
for row in cur.execute(task):
    proj_ID, nt, org, compl, slen, plp_status, plasmidfinder = [str(r) for r in row]
    org = ' '.join(org.split(' ')[:2])
    org = org.replace('virus', 'phage')
    # plasmid_host[nt] = org
    improve_seq_names[nt] = proj_ID
    plasmid_host[nt] = proj_ID + ' ' + nt + ' ' + org
    if compl == 'complete' and int(slen) < 2000000:
        collect_data.append([nt, org, int(slen), plp_status, plasmidfinder])
print('MF356679', plasmid_host['MF356679'], improve_seq_names['MF356679'])

rotated_matches = glob.glob(outfold + '*')

matches = main + 'matches.txt'
matches = open(matches).read().strip().split('\n')
matches = [m.split('\t') for m in matches]

#rotate
if 0:
    for m in matches[:]:
        ref_name,  insert_name, insert_host, to_rotate_name, rotate_host = m
        print(*m)
        insert = [i for i in inserts if insert_name in i][0]
        to_rotate = [p for p in to_process if to_rotate_name in p][0]
        
        # insert_upd = insert.replace('(', '').replace(')', '').replace("'", '')
        # if insert_upd != insert:
        #     cmd = f'mv "{insert}" {insert_upd}'
        #     run(cmd)
        #     insert = insert_upd
        
        # print(insert, to_rotate, sep = '\n')
        
        coords = blastn + to_rotate_name + '_vs_' + insert_name + '.txt'
        blastn_cov(insert, to_rotate, coords)
        
        
        if coords not in glob.glob(blastn + '*'):
            print(insert_name, 'lost coords!', coords)
        else:
            seq = ''.join(open(to_rotate).read().strip().split('\n')[1:])
            slen = len(seq)
            
            coords = open(coords).read().strip().split('\n')
            coords = [c.split('\t') for c in coords]
            coords = [c for c in coords if int(c[6]) > 1500]
            coords.sort(key = lambda x: int(x[9]))
            first = coords[0]
            
            if int(first[9]) >= 100:
                print(insert_name, to_rotate_name, 'start problem!', *first, sep = '\n')
                print(*[c[9:11] for c in coords[0:5]], sep = '\n')
                print('')
                
    #       print(*first)
            new = ''
            flag = ''
            c1, c2 = [int(f) for f in first[11:13]]
            if c2 < c1:
                flag = 'rev'
                
            new = rotate_around_coord(seq, c1, flag, rev=flag)
            if len(new) != slen:
                print(insert_name, to_rotate_name, 'new length problem!', len(new), slen)
            
            name = f'{to_rotate_name}_vs_{insert_name}_c{c1}{flag}'
            print(name)
            with open(f'{outfold}{name}.fasta', 'w') as fh:
                fh.write(f'>{name}\n')
                fh.write(dna_to_fasta(new) + '\n')
        
        print('--\n')
        
#geneplotr
if 1:
    for m in matches:
        ref_name, insert_name, insert_host, to_rotate_name, rotate_host = m
        
        fut_name = '_'.join([ref_name, insert_name, to_rotate_name])
        
        # plasmid_host = {}
        # plasmid_host[insert_name] = 'PLP_' + insert_name + '_' + insert_host.replace(' ', '_')
        # plasmid_host[to_rotate_name] = to_rotate_name + '_' + rotate_host.replace(' ', '_')
        
        
        plasmid_host[to_rotate_name] = to_rotate_name + ' ' + rotate_host
        
        
        ref_seq = [r for r in ref_seqs if '_' + ref_name + '_' in r]
        if ref_seq:
            ref_seq = ref_seq[0]
        else:
            ref_seq = '/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/P1/rotate/CP035313_P1_c78487.fasta'
        PLP_plas = [i for i in PLP_rotate if insert_name in i][0]
        print(to_rotate_name)
        add_match = [p for p in rotated_matches if to_rotate_name in p][0]
        
        fasta_files = [ref_seq, PLP_plas, add_match]
        
        print(ref_seq, PLP_plas, add_match, sep = '\n')
        ref_name = ref_seq.split('/')[-1].split('_ref')[0]
        ref_name = ref_name.split('_')[0]
        
        order = [ref_name, insert_name, to_rotate_name]
        print(order)
        
        
        gb_load = []
        gb_fill = []
        gb_list = []
        offsets = []
        names = []
        blast_comp = []
        blast_list = []
        
        
        count = 1
        gbl = 'gb{} <- read_dna_seg_from_file("{}"))'
        gbf = 'gb{}$fill <- gb{}$col'
        for o in order:
            if o in plasmid_host:
                new_o = plasmid_host[o]
            else:
                new_o = o
            
            gb_name = f'gb{count}'
            names.append(f'"{new_o}"')
            # o_file = [g for g in genbank_files if o in g][0]
            if o == order[-1]:
                o_file = [g for g in genbank_files if o in g and viz_fold in g][0]
            else:
                print(o)
                o_file = [g for g in genbank_files if o in g and viz_fold not in g][0]
            print(o_file)
            
            gbl = f'{gb_name} <- read_dna_seg_from_file("{o_file}")'
            gb_load.append(gbl)
            
            gbf = f'{gb_name}$fill <- {gb_name}$col'
            gb_fill.append(gbf)
            
            gb_list.append(gb_name)
            offsets.append('0')
            
            if o != order[-1]:
                query = [f for f in fasta_files if order[count-1] in f][0]
                subject = [f for f in fasta_files if order[count] in f][0]
                output = blast_fold + order[count-1] + '_vs_' + order[count] + '.blast'
                
                if True or output not in glob.glob(blast_fold + '*') or 'CP056516' in fut_name:
                    blast(query, subject, output)
                    
                    output_proc = open(output).read().strip().split('\n')
                    output_proc = [o.split('\t') for o in output_proc] 
                    output_proc = ['\t'.join(o) for o in output_proc if int(o[3]) > 1000]
                    
                    with open(output, 'w') as fh:
                        fh.write('\n'.join(output_proc) + '\n')
                    print(output)
                
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
        
        tree_width=3
        width = 20
        height = 3*1.1
        annotation_cex=1
        
            
        main_name = f'{fut_name}'
        
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
        
        blocks = [start, gb_load + gb_fill + names + gb_list + offsets + annots, blast_comp + blast_list, save_image]
        blocks = [b for b in blocks if b]
        blocks = ['\n\n'.join(blocks)]
        
        blocks = ['```{r}\n' + b + '\n\n```\n\n' for b in blocks]
        r = ''.join(blocks)
        
        
        figure = viz_fold + main_name + '.svg'
        
        with open(f'{viz_fold}{main_name}.rmd', 'w') as fh:
            fh.write(r + '\n')
        
        print('--\n')