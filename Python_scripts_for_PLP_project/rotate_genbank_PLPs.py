#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 03:06:18 2019

@author: lutra
"""
from Bio import SeqIO
import subprocess
import glob
import re
import os

def run(cmd):
    process = subprocess.Popen(cmd, shell=True)
    process.wait()


def rev_comp(dna):
    trantab = "".maketrans("ACTGactg","TGACtgac")
    return dna[::-1].translate(trantab)


def dna_to_fasta(fasta_file, name, dna):
    dna = [dna[i*60:(i+1)*60] for i in range(int(len(dna)/60)+1)]
    dna=[d for d in dna if d]
    dna = '\n'.join(dna)
    with (open(fasta_file, 'w')) as fh:
        fh.write('>{}\n{}\n'.format(name, dna))
    return


def blast(params, query, subject, output):
    cmd='blastn -query %s -subject %s -out %s -outfmt "6 %s" -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.000001' %(query, subject, output, params)
    run(cmd)
    result = ((open(output).read()).strip()).split('\n')
    return result


def regular_blast(query, subject, output):
    cmd='blastn -query %s -subject %s -out %s -outfmt 6 -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 0.000001' %(query, subject, output)
    run(cmd)
    return


def tblastn(params, query, subject, output):
    cmd='tblastn -query %s -subject %s -out %s -outfmt "6 %s" -evalue 0.000001' %(query, subject, output, params)
    run(cmd)
    result = ((open(output).read()).strip()).split('\n')
    return result

def blastx(params, query, subject, output):
    cmd='blastx -query %s -subject %s -out %s -outfmt "6 %s" -evalue 0.000001' %(query, subject, output, params)
    run(cmd)
    result = ((open(output).read()).strip()).split('\n')
#    print('BLASTX', subject)
#    print(*result, '\n', sep = '\n')
    return result


def blast_the_start(query, subject, output):
    params = 'qstart qend sstart send qcovs'
    result = blast(params, query, subject, output)
    result = [r.split() for r in result]
    result = [[int(rr) for rr in r] for r in result]
    result.sort()
    print (result[0])
    c1, c2, qcovs = result[0][2:]
    if c1 < c2:
        return (c1, 0 ) #cut from the left side, add
    else:
        return (c1, 1 ) #cut from the right side, add and revert
   
    
def blast_db(query, subject, output, program='blastn'):
    params = 'pident slen length qstart qend sstart send saccver qseq'
    if program == 'blastn':
        result = blast(params, query, subject, output)
    elif program == 'blastx':
        result = blastx(params, query, subject, output)
    best = []
    if result != ['']:
        result = ((open(output).read()).strip()).split('\n')
        result0 = [r.split() for r in result]
        result = [r for r in result0 if float(r[0])>=90 and int(r[2])/int(r[1])>=0.6]
        if 'isfinder' in subject:
            result = [r for r in result0 if float(r[0])>=80 and int(r[2])/int(r[1])>=0.2]
        
#        print (*result, sep='\n')
        matches = []
        for r in result:
            pident = float(r[0])
            slen, length, qstart, qend, sstart, send = [int(rr) for rr in r[1:-2]]
            saccver, qseq = r[-2:]
            cov = length/slen
            nap = '+'
            if sstart > send:
                nap = '-'
            matches.append([pident, cov, saccver, qstart, qend, nap, pident, round(100*cov, 1), qseq])
        
        best = [m for m in matches if m[0] == 100 and m[1] == 1]
        # print ('check No of matches', len(matches), end = ' ')
        trunk = [m[:2] + ['**' + m[2]] + m[3:] for m in matches if m not in best and m[1] <= 0.95]
        muts = [m[:2] + ['#' + m[2]] + m[3:] for m in matches if m not in best and m[1] > 0.95]
        matches = muts + trunk
        # print (len(trunk), len(muts), len(matches) + len(best))
        
        if 'isfinder' in subject:
            matches = best + matches
            best = []
            
        matches.sort(key=lambda x: [x[1], x[0]], reverse=True)
        
        for m in matches:
            m1, m2 = m[3:5]
            add = True
            for b in best:
                b1, b2 = b[3:5]
                if m[5] == b[5]:
                    if (b1<=m1 and m2<=b2) or (b1>=m1 and m2>=b2):
                        add = False
            if add:
                best.append(m)
#        print ('!!!!')
        # print (*best, sep='\n')
    best = [b[2:] for b in best]
    return best


def generate_gb(new_gb, name, host, size, feats, dna, easyfig=False):
    gb=open(new_gb,'w')
    first=f'LOCUS       {name}              {size} bp    DNA\nDEFINITION  {host} {name}, complete sequence\nFEATURES             Location/Qualifiers\n     source          1..{size}\n'
    gb.write(first)
    
    color_code = {}
    colors = '/data/Current_work/Phage_like_plasmids/mauve/fasta/trees/templates/rgb_colors.txt'
    colors = open(colors).read().strip().split('\n')
    colors = [c.split('\t') for c in colors]
    
    for c in colors:
        if easyfig:
            color_code[c[0]] = c[2].replace('(', '').replace(')', '').replace(',',' ') #easyfig
        else:
            color_code[c[0]] = c[0] #genoplotr
    
    for f in feats:
        el, c1, c2, strand, color = f
        if 'block' not in el:
            el = (el.split('_'))[0]
            if el[:3] == 'bla':
                el = el[3:]
        else:
            if not easyfig:
                el = ''
            else:
                pass
        
        if el[:1] == '#':
            el = el[1:] + '#'
        coords = '{}..{}'.format(c1, c2)
        
        if strand == '-':
            coords='complement('+coords+')'
            
        if color not in ():
            print(color, color_code[color])
            
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



def generate_gb_old(new_gb, name, size, feats, dna):
    gb=open(new_gb,'w')
    first='LOCUS       '+name+'              '+str(size)+' bp    DNA\nFEATURES             Location/Qualifiers\n     source          1..'+str(size)+'\n'
    gb.write(first)
    color_code = {}
    color_code['gb']='128 128 128' #gray
    color_code['arg']='255 0 0' #red
    color_code['Inc']='128 0 128' #purple
    color_code['phage'] = '255 204 102' #sand
    color_code['hypro']='230 230 230'#light gray
    color_code['IS']='255 140 0' #orange
    
    
    color_code['antitoxin']='0 255 0' #light green
    color_code['tra']='0 191 255' #light blue
    color_code['dna']='128 0 128' #purple
    color_code['UV_res']='0 255 255' #aqua
    color_code['edges']='180 180 180'
    for f in feats:
        el, c1, c2, strand, color = f
        el = (el.split('_'))[0]
        if el[:1] == '#':
            el = el[1:] + '#'
        coords = '{}..{}'.format(c1, c2)
        if strand == '-':
            coords='complement('+coords+')'
        if color not in ():
            entry='     CDS             '+coords+'\n                     /locus_tag="'+el+'"\n                     /colour='+color_code[color]+'\n'
            gb.write(entry)
        else:
            pass
#            entry='     CDS             '+coords+'\n                     /colour='+color_code[color]+'\n'
#            gb.write(entry)
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
    gb.write('//')
    gb.close()
    return


def Easyfig(query, figure, fold):
    easyfig_cmd = 'Easyfig.py -svg -width 8821 -ann_height 250 -blast_height 1000 -f1 T -f2 10000 -min_length 1500 -aln left -blast_col 0 191 255 0 0 255 -blast_col_inv 255 215 0 255 140 0 -bo F -f CDS 0 0 0 rect -glt 5 -genet 20 -legend both -leg_name locus_tag -uncomp T '
    easyfig_cmd = 'python2.7 {}Easyfig-lutra/{}'.format(['/home/shiri/plasmid_project/tools/', '/data/Bioinformatics'][0], easyfig_cmd)
    
    cmd = f'{easyfig_cmd} -o {figure} {query}'
    run(cmd)
    
    for q in query.split():
        q_name = q.split('/')[-1].split('.')[0]
        cmd = f'{easyfig_cmd} -o {fold}{q_name}.svg {q}'
        run(cmd)
    return



if __name__ == '__main__':
    if 0:
        print('Run rotate fasta!!!')
        main = '/data/Current_work/Phage_like_plasmids/mauve/fasta/'
        
        proj = 'SSU5'
        ref_name = {'P1': 'AF234172', 'D6': 'MF356679', 'SSU5': 'JQ965645'}
        ref = main + ref_name[proj] + '.fasta'
        
        all_fasta_fold = f'{main}complete_plasmids_{proj}/'
        all_fasta = glob.glob(all_fasta_fold + '*')
        
        new_fold = f'{main}/HomBlocks_{proj}/'
        cmd = 'mkdir ' + new_fold
        run(cmd[:-1])
        
        new_fold = new_fold + f'complete_plasmids_{proj}_rotate/'
        cmd = 'mkdir ' + new_fold
        run(cmd[:-1])
        
        output = main +'output.txt'
        
        for fasta in all_fasta:
            name = fasta.split('/')[-1].split('.')[0] + '_rotate'
            coord, reverse = blast_the_start(ref, fasta, output)
            
            new_fasta = new_fold + name + '.fasta'
            
            if reverse:
                dna = ''.join(open(fasta).read().strip().split('\n')[1:])
                print('rotate!', len(dna), dna[:10])
                dna = rev_comp(dna)
                dna_to_fasta(new_fasta, name, dna)
                coord, reverse = blast_the_start(ref, new_fasta, output)
            print ('')
            
            dna_final = ''
            feats_final = []
            
            dna_final = dna[coord-1:] + dna[:coord-1]
            dna_to_fasta(new_fasta, name, dna_final)
        
        os.remove(output)
        cmd = f'cp {ref} {new_fold}{ref_name[proj]}_{proj}_ref.fasta'
        run(cmd)

        
    
    
    if 0:
        print ('Run rotate_genbank!!!')
        
        main = '/data/Current_work/Phage_like_plasmids/Nucleotide_db/rotate_genbank/'
        main2 = '/home/shiri/plasmid_project/Phage_like_plasmids/Nucleotide_db/rotate_genbank/'
        main = '/home/shiri/plasmid_project/Phage_like_plasmids/SSU5/'
        main = '/home/shiri/plasmid_project/Phage_like_plasmids/prophages/rotate_genbank/'
    
        all_gb = main + 'all_sequence_20Nov19.gb'
        all_gb = main + 'all_gb.gb'
        all_gb = main + 'ARG_prophages.gb'
        
        task_file = main + 'Group3_dist_heatmap_clusters_17.txt'
        task_file = main + 'SSU5_task.txt'
        task_file = main + 'ARG_prophages_task.txt'
        
        task = open(task_file).read().strip().split('\n')
        task = [t.split() for t in task]
        task = [[tt.split('.')[0] for tt in t] for t in task]
        for t in task:
            print(*t)
        
        task_name = task_file.split('/')[-1].split('.')[0]
        gbk_fold = main + task_name + '/'
        if gbk_fold[:-1] not in glob.glob(main + '*'):
            run(f'mkdir "{gbk_fold[:-1]}"')
        
        
        output = main +'output.txt'
        
        for t in task[2:]:
            ref_name = t[0]
            
            t_fold = gbk_fold + ref_name + '/'
            run(f'mkdir "{t_fold[:-1]}"')
            
            figure = gbk_fold + ref_name + '_as_' + task_name + '.svg'
            query = []
            
            figure_labels = gbk_fold + ref_name + '_as_' + task_name + '.txt'
            labels = []
            
            for gid in t[:]:
                for gb_record in SeqIO.parse(open(all_gb, "r"), "genbank") :
                    if gid == gb_record.id.split('.')[0]:
                        name = gb_record.id + ' ' + gb_record.description
                        print ("Name %s, %i features" % (name, len(gb_record.features)))
                        
                        dna = str(gb_record.seq)
                        size = len(dna)
                        
                        gb_fasta = t_fold + gid + '.fasta'
                        new_gb = t_fold + gid + '.gb'
                        query.append(new_gb)
                        
                        
                        dna_to_fasta(gb_fasta, name, dna)
                        ref = t_fold + ref_name + '.fasta'
                        
                        coord, reverse = blast_the_start(ref, gb_fasta, output)
                        
                        feats = []
                        for feat in gb_record.features[1:]:
                            if feat.type == 'CDS':
                                gene = 'orf'
                                try:
                                    gene = feat.qualifiers['gene'][0]
                                except:
                                    try:
                                        gene = feat.qualifiers['product'][0]
                                    except:
                                        pass
                                
                                loc = str(feat.location)
                                c1, c2 = (re.findall('[0-9]+', loc))[0], (re.findall('[0-9]+', loc))[-1]
                                c1, c2 = int(c1)+1, int(c2)                
                                strand = (re.findall('[+-]', loc))[0]
                                if reverse:
                                    rev_strand = {'+':'-', '-':'+'}
                                    strand = rev_strand[strand]
                                    c1 = size - c1 + 1
                                    c2 = size - c2 + 1
                                    c1, c2 = c2, c1    
                                feats.append([gene, c1, c2, strand])
                                
                        if reverse:
                            print ('rotate!')
                            dna = rev_comp(dna)
                            dna_to_fasta(gb_fasta, name, dna)
                            coord, reverse = blast_the_start(ref, gb_fasta, output)
                        print ('')
                        
                        dna_final = ''
                        feats_final = []
                        
                        dna_final = dna[coord-1:] + dna[:coord-1]
                        dna_to_fasta(gb_fasta, name, dna_final)
                        
                        shift_size = len(dna[coord-1:])
                #        print('Use BLAST!', dna_final[:240])
                            
                        feats_final = [] 
                        for f in feats:
                            gene, c1, c2, strand = f
                            c1 += shift_size
                            c2 += shift_size
                            if c1 > size:
                                c1 -= size
                            if c2 > size:
                                c2 -= size
                            if c2 > c1:
                                feats_final.append([gene, c1, c2, strand, 'gb'])
                            else:
                                print ('cut', f, c1, c2)
                                feats_final.append([gene, c1, size, '>'+strand, 'gb'])
                                feats_final.append([gene, 1, c2, '<'+strand, 'gb'])
                        
                        
                        label = [gb_record.id, gb_record.description, str(size) + 'bp']
                        all_spec = []        
                        args = main2 + 'resfinder_db/*'
                        args = [a for a in glob.glob(args) if '.fsa' in a]
                        add_to_label = []
                        for arg in args[:]:
                            color = 'arg'
                            add = blast_db(gb_fasta, arg, output)
                            add = [[a[0].replace('bla', '')] + a[1:] for a in add]
                            add_to_label += [a[0] for a in add]
                            add = [a+[color] for a in add]
                            all_spec += add
                        label.append(' '.join(add_to_label))
                        print ('Total ARGs:', len(all_spec))
                        print ('-\n\n')
                        
                        incs = main2 + 'replicon_db/*'
                        incs = [a for a in glob.glob(incs) if '.fsa' in a]
                        add_to_label = []
                        for inc in incs[:]:
                            color = 'Inc'
                            add = blast_db(gb_fasta, inc, output)
                            add_to_label += [a[0] for a in add]
                            add = [a + [color] for a in add]
                            all_spec += add
                        label.append(' '.join(add_to_label))
                            
                        isels = main + 'isfinder_db/*'
                        isels = [a for a in glob.glob(isels) if '.fsa' in a]
                        add_to_label = []
                        for isel in isels:
                            color = 'IS'
                            add = blast_db(gb_fasta, isel, output)
                            add_to_label += [a[0] for a in add]
                            add = [a+[color] for a in add]
                            all_spec += add
                        label.append(' '.join(add_to_label))
    #                        
    #                    tra = main + 'groups_db/*'
    #                    tra = [a for a in glob.glob(tra) if '.fsa' in a]
    #                    for tr in tra:
    #                        color = (((tr.split('/'))[-1]).split('.'))[0]
    #                        add = blast_db(gb_fasta, tr, output)
    #                        add = [a+[color] for a in add]
    #                        all_spec += add
    
                        label = '\n'.join(label)
                        labels.append(label)
    
                        for to_add in (feats_final,): 
                            for feat in to_add:
                                m1, m2 = feat[1:3]
                                add = True
                                for b in all_spec:
                                    b1, b2 = b[1:3]
                                    if feat[3] == b[3]:
                                        if (m1>=b1 and m2<=b2):
                                            add = False
                                            break
                                        elif (m1<=b1 and m2>=b2):
                                            pass
                                        elif (m1<=b1 and m2>=b1):
                                            feat[2] = b1
                                        elif (m1<=b2 and m2>=b2):
                                            feat[1] = b2
                                if add:
                                    if 'phage' in feat[0].lower():
                                        feat[-1] = 'phage'
                                    if 'hypothetical protein' in feat[0].lower():
                                        feat[-1] = 'hypro'
                                        
                                    print ('add', feat)
    #                                remove all gene names!!!
                                    feat[0] = ' '
                                    all_spec.append(feat)
                
                        all_spec.sort(key = lambda x: x[1])
                        feats_final.sort(key = lambda x: x[1])
                        
                        if False:
                            generate_gb(new_gb, name, size, feats_final, dna_final) #all features
                        else:
                            generate_gb(new_gb, name, size, all_spec, dna_final) #only sprecial functional genes from databases
                        print ('--\n\n')
            
            query = ' '.join(query)
            Easyfig(query, figure, t_fold)
            
            with open(figure_labels, 'w') as fh:
                fh.write('\n\n'.join(labels))
                