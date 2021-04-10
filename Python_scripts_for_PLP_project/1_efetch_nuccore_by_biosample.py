#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 18:11:45 2019

@author: lutra
"""
import subprocess
import json
import glob


def run(cmd):
    print(cmd)
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
    print('-\n')
    

def mkdir(fold):
    cmd = 'mkdir ' + fold[:-1]
    run(cmd)

main = '/home/shiri/plasmid_project/Other_projects/Salmonella_p373/PCR_targets/Nucleotide_db/'
main = '/home/shiri/plasmid_project/Phage_like_plasmids/Nucleotide_db/'
main = '/home/shiri/plasmid_project/Phage_like_plasmids/all_EC_phage_like_plasmids/'
main = '/home/shiri/plasmid_project/Phage_like_plasmids/prophages/'
main = '/data/Current_work/Phage_like_plasmids/IncRep_relationships/replicon_blastn/P1/'
main = '/home/shiri/plasmid_project/Phage_like_plasmids/IncRep_relationships/replicon_blastn/P1/'
main = '/home/shiri/kpn_project/All_Kpn_plasmids/'
main = '/home/shiri/plasmid_project/Phage_like_plasmids/IncRep_relationships/replicon_blastn/P1/'
main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'

main = '/home/shiri/plasmid_project/Phage_like_plasmids/IncRep_relationships/SSU5_associated_replicons/'
main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'

main = '/home/shiri/plasmid_project/ST131_genomes/'
main = '/data/Current_work/ST131_genomes/'

all_nuccore_seqs = False

# if all_nuccore_seqs:
#     file = main + 'group2_biosamples.txt'
# else:
#     file = main + 'entrez.txt'
#     file = main + 'both_db_entrez.txt'
#     file = main + 'JQ965645l_entrez.txt'
#     file = [f for f in glob.glob(main + '*') if 'Entrez' in f]
#     print(*file)

file = [main + 'new_Entrez_26Aug20.txt']
file = [main + 'Entrez_add_missed_replicons_for_all_groups.txt']
file = [main + 'Entrez_all_3Sep20.txt']
file = [main + 'All_EC_1Mb_100Mb_13Dec20.seq']
file = [main + 'Entrez_Ecoli_complete_genomes.txt']
file = [main + 'Entrez_all_found_phages.txt']
# file = [main + 'add.txt']


if False:    
    outf = main + 'fasta/'
    mkdir(outf)
    
    ready = glob.glob(outf + '*')
    
    accs = []
    for f in file:
        accs += ((open(f).read()).strip()).split('\n')
    print('going to retrive this number of seqs: ', len(accs))
    
    for acc in accs[:]:
        if '{}{}.fasta'.format(outf, acc) not in ready:
            cmd = 'esearch -db nuccore -query "{}" | efetch -format fasta > {}{}.fasta'.format(acc, outf, acc)
            print(acc)
            run(cmd)
    
        
if True:
    metadata = main + 'matched_phage_metadata/'
    mkdir(metadata)
    
    ready = glob.glob(metadata + '*')
    
    
    qaccvars = []
    for f in file:
        qaccvars += ((open(f).read()).strip()).split('\n')
    qaccvars.sort()
    for qacc in qaccvars:
        if '{}{}.json'.format(metadata, qacc) not in ready:
            cmd = 'esearch -db nuccore -query "{}" | esummary -mode json > {}{}.json'.format(qacc, metadata, qacc)
            run(cmd)
    
    res = {}
#    feats = ['nucleotide', 'biosample', 'organism', 'BioSample_title', 'completeness', 'genome', 'slen']
    feats = ['nucleotide', 'biosample', 'organism', 'completeness', 'genome', 'slen']
    friz = len(feats)
    jsons = [j for j in glob.glob(metadata + '*') if j[-5:] == '.json']
    for js in jsons[:]:
        qacc = (((js.split('/'))[-1]).split('.'))[0]
        res[qacc] = {}
        for f in feats[1:friz]:
            res[qacc][f] = []
        
        with open(js) as js_file:
            try:
                load = json.load(js_file)
            except:
                print(js)
                exit
            
            uids = load['result']['uids']
            for uid in uids:
                # print(qacc, uid)
                for f in feats[1:friz]:
                    if f in load['result'][uid]:
                        res[qacc][f].append(str(load['result'][uid][f]))
                
                # print(load['result'][uid])
                if 'subtype' in load['result'][uid] and 'subname' in load['result'][uid]:
                    cols = load['result'][uid]['subtype']
                    data = load['result'][uid]['subname']
                    pairs = zip(cols.split('|'), data.split('|'))
                    
                    for p in pairs:
                        # print(p)
                        col, dat = p
                        col = col.replace(' ', '_')
                        col = col.lower()
                        if col not in feats:
                            feats.append(col)
                        if col not in res[qacc]:
                            res[qacc][col] = []
                        res[qacc][col].append(dat)
                    # print('--\n\n')
        
        for bios in res[qacc]['biosample']:
            if bios:
                runinfo = '{}{}_{}.txt'.format(metadata, qacc, bios)
                if runinfo not in ready:
                    cmd = 'esearch -db biosample -query "{}" | efetch -format runinfo > {}'.format(bios, runinfo)
                    run(cmd)
                
                addinfo = [a for a in ((open(runinfo).read()).strip()).split('\n') if '   /' in a or '1: ' in a]
                for add in addinfo:
                    if ': ' in add:
                        pass
    #                    res[qacc]['BioSample_title'].append((add.split(': '))[1])
                    else:
                        col = (((add.split('/'))[1]).split('='))[0]
                        col = col.replace(' ', '_')
                        col = col.lower()
                        dat = (add.split('"'))[1]
                        if col not in feats:
                            feats.append(col)
                        if col not in res[qacc]:
                            res[qacc][col] = []
                        res[qacc][col].append(dat)
#        print('')
    
    with open(main + 'results.csv', 'w') as fh:
        fh.write('\t'.join(feats) + '\n')
        for qacc in qaccvars:
            qacc = (qacc.split('.'))[0]
            entry = [qacc]
            info = res[qacc]
            for feat in feats[1:]:
                if feat in info:
                    add = ';'.join(set(info[feat]))
                else:
                    add = ''
                entry.append(add)
            fh.write('\t'.join(entry) + '\n')
                
