#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 00:32:48 2020

@author: lutra
"""
import glob
import os
from rotate_genbank_PLPs import blast_db, run
from datetime import datetime

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


def annotate_query(query, main, db, outf, suff = '', rerun=False, this_file = False):
    global collect_diffs    
    mkdir(outf)
    
    check = glob.glob(outf + '*')
    name = query.split('/')[-1].split('.')[0]
    final_output = outf + name + suff +'.txt'
    
    
    if final_output not in check or rerun:
    
        output = main + 'output.txt'
        prolonged_query = main + 'long_query.fasta'
        
#        print(*db, '--', sep = '\n')
#        print(query)
        record = []
        
        query_seq = open(query).read().strip().split('\n')
        dna = ''.join(query_seq[1:])
        orig_size = len(dna)
#        print(orig_size)
        
        dna += dna[:10000]
        with open(prolonged_query, 'w') as fh:
            fh.write(f'>{query_seq[0]}\n{dna}')
    
        for subj in db[:]:
            # print(subj)
            group = subj.split('/')[-1].split('.')[0]
            
            program = 'blastn'
            if 'blastx' in subj:
                program = 'blastx'
            
            res = blast_db(prolonged_query, subj, output, program)
            if res:
                more = [r for r in res if r[2] > orig_size]
                # print('\n\nafter end', *more, '---\n\n', sep = '\n')
                to_add = [r for r in more if r[1] <= orig_size]
                
                res = [r for r in res if r[2] <= orig_size]
                for add in to_add:
                    a1, a2 = add[1], add[2]
                    after_end = f'decircle {a1}..{orig_size} 1..{int(a2) - orig_size}'
                    res.append(add + [after_end])                
                
                res = [[group] + r for r in res]
                record += res
          
        record = ['\t'.join(map(str, r)) for r in record]  
        if this_file and rerun and final_output in check:
            curr_records = open(final_output).read().strip().split('\n')
            lost = ['LOST\t' + c for c in curr_records if c and c not in record]
            found = ['FOUND\t' + r for r in record if r and r not in curr_records]
            entry = ''
            if lost:
                entry += ('\n'.join(lost)+'\n')
            if found:
                entry += ('\n'.join(found)+'\n')
            if entry:
                with open(collect_diffs, 'a') as fh:
                    fh.write(name + '\n' + entry + '-----------------\n\n')
        
        # print(final_output)
        # print(*record, sep = '\n')
        with open(final_output, 'w') as fh:
            fh.write('\n'.join(record) + '\n')
    
        os.remove(output)
        os.remove(prolonged_query)
        # print('annotate_query has finished')
    return
    
    
def coverage_query(subj, main, outf, suff = '', rerun=False):
    db = main + 'coverage_db/'
    mkdir(outf)
    
    check = glob.glob(outf + '*')
    name = subj.split('/')[-1].split('.')[0]
    final_output = outf + name + suff +'.txt'
    
    if final_output not in check or rerun:
        db = glob.glob(db + '*')
        db = [d for d in db if '.fasta' in d]
#        print(*db, '--', sep = '\n')
        
        suff = ''
        
        output = main + 'output.txt'
        
        collect = []
        for ref_pl in db:
            ref_pl_name = ref_pl.split('/')[-1].split('.')[0]
            
            cmd = f'blastn -query {ref_pl} -subject {subj} -out {output} -outfmt "6 qcovs" -word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2'
#            print(cmd)
            run(cmd)
            
            res = open(output).read().strip().split('\n')[0]
            collect.append('\t'.join([ref_pl_name, res]))
        
        with open(outf + name + suff +'.txt', 'w') as fh:
            fh.write('\n'.join(collect) + '\n')
            
        os.remove(output)
    return


if __name__ == '__main__':
    main = '/home/shiri/kpn_project/All_Kpn_plasmids/'
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/IncRep_relationships/SSU5_associated_replicons/'
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
    
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'
    
    fasta = main + 'fasta/'
    annotate = main + 'annotate/'
    coverage = main + 'coverage/'
    
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    collect_diffs = main + f'Lost_n_found_{dt_string}.txt'.replace('/','_').replace(':','_').replace(' ','_')
    with open(collect_diffs, 'w') as fh:
        fh.write(dt_string + '\n')
    
    
    fasta = glob.glob(fasta + '*')
    fasta = [f for f in fasta if 'test' not in f]
    #fasta = [f for f in fasta if 'CP025870' in f]
    fasta = [f for f in fasta if f.split('/')[-1].split('.')[0] in ['CP012254', 'CP014707', 'CP053410', 'KP763470', 'CP030187', 'CP044180', 'CP057631', 'CP057633', 'KY515226']]
    
    
    db = main + 'db/'
    db = glob.glob(db + '*')
    db = [f for f in db if '.fsa' in f] + ' '.join([' '.join(glob.glob(f+'/*')) for f in db if '.fsa' not in f]).split()
    db = [f for f in db if '.fsa' in f]
    
    print(len(db), *db, sep = '\n')
    
    print(len(fasta))
    count = 0
    for f in fasta[:]:
        count += 1
        print(f, count, str(round(100*count/len(fasta), 1))+'%')
        annotate_query(f, main, db, annotate, rerun=False, this_file=True)
        coverage_query(f, main, coverage, rerun=False)
        print('--\n\n')
        
    
    #for plas in ['P1', 'SSU5']:
    #    fold = f'/data/Current_work/Phage_like_plasmids/IncRep_relationships/replicon_blastn/{plas}/'
    #    fold = f'/home/shiri/plasmid_project/Phage_like_plasmids/IncRep_relationships/replicon_blastn/{plas}/'
    #    annotate_all(fold)
    
