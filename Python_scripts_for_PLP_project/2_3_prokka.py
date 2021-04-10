#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 19:44:58 2020

@author: lutra
"""
import glob
import PLP_main_functions as myf


PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
else:   
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'


proj = 'SSU5'
print('PROJECT!', proj, main, '\n')


plas_fold = f'{main}prophages_groups/{proj}/'
cluster_fold = plas_fold + 'HomBlocks_per_cluser/'
clusters = glob.glob(cluster_fold + '*')
clusters = [c +'/' for c in clusters if '.' not in c and c.split('/')[-1] not in ('Easyfig', 'Gepard_dotplot', 'iqtree', 'Viz_GenBank')]
clusters.sort()

TA_systems = ['higA', 'higB', 'hic', 'ccd', 'phd', 'pem', 'vap', 'yaf', 'hip', 'rel', 'doc', 'higB-1'] 

phage_ref = ['AF234172', 'MF356679', 'JQ965645']

# clusters = [c for c in clusters if 'rotated_ARG' not in c and 'mcr' not in c]
# clusters = [c for c in clusters if 'rotated_ARG' not in c and 'ARG_extra' not in c and 'mcr' not in c]
# clusters = [c for c in clusters if 'rotated_ARG' in c]


for clust in clusters[:]:
    print(clust)
    plas_group = clust.split('/')[-2]
    fasta = glob.glob(clust + 'fasta/*')
    
    print(plas_group, len(fasta))
    
    prokka = f'{clust}prokka/'
    myf.mkdir(prokka)
    
    if glob.glob(prokka + '*'):
        for f in fasta[:]:        
            name = f.split('/')[-1].split('.')[0]
            phage_ref_check = [p for p in phage_ref if p in name]
            
            cmd = f'prokka --outdir {prokka}{name} --genus PLP_genbanks_upd --prefix {name} --kingdom Bacteria --locustag {name} {f}'
            
            cmd = f'prokka --outdir {prokka}{name} --prefix {name} --kingdom Bacteria --locustag {name} {f}'
            myf.run(cmd)
            
            cmd = f'mv {prokka}{name}/{name}.gff {prokka}{name}.gff'
            myf.run(cmd)
            
            cmd = f'rm -r {prokka}{name}'
            myf.run(cmd)
            
            file = f'{prokka}{name}.gff'
            feats = [f.split('\t') for f in open(file).read().split('##FASTA')[0].strip().split('\n') if f[:1] != '#']
            new_feats = []
            for f in feats:
                c1, c2 = f[3:5]
                nap = f[6]
                prod = [p for p in f[8].split(';') if 'product=' in p][0].split('=')[-1]
                gene = [p for p in f[8].split(';') if 'gene=' in p]
                if gene:
                    gene = gene[0].split('=')[-1]
                else:
                    gene = ''
                    
                prod_type = ''
                
                if prod == 'hypothetical protein':
                    prod_type = 'hypothetical'
                    
                elif f[2] == 'tRNA':
                    prod_type = 'tRNA'
                    gene = prod.split('-')[-1]
                    print(f[2], prod, gene)
                    
                elif 'toxin' in prod and 'enterotoxin' not in prod and 'Dermonecrotic' not in prod or prod in TA_systems or gene in TA_systems: #or prod[:3].lower() in TA_systems or 
                    prod_type = 'toxin'
                    if 'toxin ' in prod and not gene:
                        gene = prod.split('toxin ')[-1]
                    print('toxin', prod, gene)
                
                # elif '___replicon' in prod:
                #     print(prod, '->', gene)
                #     prod_type = 'PlasmidFinder'
                
                elif '___tail' in prod:
                    prod_type = 'tail'
                    if not phage_ref_check:
                        gene = ''
                        
                elif '___transfer' in prod or '___conjug' in prod:
                    prod_type = 'oriTDB'
                    
                    
                elif '___phage' in prod:
                    prod_type = 'phage'
                    if not phage_ref_check:
                        gene = ''
                else:
                    prod_type = 'CDS'
                    
                    if gene != 'glnS' and gene.lower() != 'pap2':
                        if not phage_ref_check:
                            gene = ''
                    else:
                        gene = gene[:3].lower() + gene[3:]
                        print(gene)
                new_feats.append([prod_type, gene, c1, c2, nap, '50', '50'])
                
            new_feats = ['\t'.join(n) for n in new_feats]
            # print(*new_feats[:10], sep = '\n')
            with open(file[:-4] + '.txt', 'w') as fh:
                fh.write('\n'.join(new_feats) + '\n')
            