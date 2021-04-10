#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 17:14:42 2020

@author: lutra

CREATE TABLE "Phage_like_plasmids_SSU5_P1_D6" (
	"id"	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
	"nucleotide"	TEXT,
	"biosample"	TEXT,
	"organism"	TEXT,
	"completeness"	TEXT,
	"genome"	TEXT,
	"slen"	INTEGER,
	"ResFinder"	INTEGER,
	"JQ965645_SSU5"	INTEGER,
	"AF234172_P1"	INTEGER,
	"p721"	INTEGER,
	"plasmidfinder"	TEXT,
	"MF356679_D6"	INTEGER,
	"strain"	INTEGER,
	"host"	INTEGER,
	"plasmid"	TEXT,
	"country"	INTEGER,
	"isolation_source"	INTEGER,
	"lat_lon"	INTEGER,
	"collection_date"	INTEGER,
	"collected_by"	INTEGER,
	"host_disease"	INTEGER,
	"latitude_and_longitude"	INTEGER,
	"geographic_location"	INTEGER,
	"Other_features"	TEXT
"""
import os
import glob
from csv_in_tables_to_sqlite3 import csv_to_sqlite3
from PLP_main_functions import find_overlaps


fold = '/home/shiri/plasmid_project/Phage_like_plasmids/IncRep_relationships/SSU5_associated_replicons/'

fold = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
fold = '/data/Current_work/Phage_like_plasmids/PLP_final/'

name = 'Phage_like_plasmids_SSU5_P1_D6_30Jun20'
name = 'SSU5_associated_replicons_12Jul20'
name = 'Phage_like_plasmids_SSU5_P1_D6_8Sep20'
name = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'

ref_list = []

res_file = fold + 'results.csv'

annotate = fold + 'annotate/'
annotate_files = glob.glob(annotate + '*')

coverage = fold + 'coverage/'
coverage_files = glob.glob(coverage + '*')

print(len(annotate_files), len(coverage_files))

res = open(res_file).read().strip().split('\n')
res = [r.split('\t') for r in res]

header = []
create = '''"project_ID"	TEXT,
    "project_ID_number"	INTEGER,
    "nucleotide"	TEXT,
	"biosample"	TEXT,
	"organism"	TEXT,
	"completeness"	TEXT,
	"genome"	TEXT,
	"slen"	INTEGER,
    "shape" TEXT,
    "PLP_status" TEXT,
	"MF356679_D6_ref_cov"	INTEGER,
	"AF234172_P1_ref_cov"	INTEGER,
	"JQ965645_SSU5_ref_cov"	INTEGER,
    "Major_replicon" TEXT,
    "Major_replicon_variant" TEXT,
    "Major_replicon_sequence" TEXT,
	"add_feats"	TEXT,
	"strain"	INTEGER,
	"host"	INTEGER,
	"plasmid"	TEXT,
	"country"	INTEGER,
	"isolation_source"	INTEGER,
	"lat_lon"	INTEGER,
	"collection_date"	INTEGER,
	"collected_by"	INTEGER,
	"host_disease"	INTEGER,
	"latitude_and_longitude"	INTEGER,
	"geographic_location"	INTEGER,
	"Other_features"	TEXT'''
create = create.replace('\n', ' ').split()
create = [c.replace('"', '') for c in create if '"' in c]
add_feats = []
print(create)

resfinder_db = [g.split('/')[-1].split('.')[0] for g in glob.glob(fold + 'db/resfinder_db/*') if 'fsa' in g]
#print(resfinder_db)



expand_feats = {}
for oritdb in glob.glob(fold + 'db/oriTDB/*') + [fold + 'db/VFDB_setB_nt.fsa']:
    oritdb_name = oritdb.split('/')[-1].split('.')[0]
    oritdb_db = [v[1:].replace(' ', '--') for v in open(oritdb).read().strip().split('\n') if v[:1] == '>']
    oritdb_dict = dict(zip([v.split('--')[0] for v in oritdb_db], oritdb_db))
    expand_feats[oritdb_name] = oritdb_dict
    
for e in expand_feats:    
    print(e, len(expand_feats[e]))


prot_db = {}
prot_descr = [g for g in glob.glob(fold + 'db/*') if '_ref' in g]
for prot0 in prot_descr:
    prot0 = open(prot0).read().strip().split('\n')
    prot0 = [p[1:].replace(']', '').split(' [') for p in prot0 if p[:1] == '>']
    
    for prot in prot0:
        lcl = prot[0]
        func = [p for p in prot if 'protein=' in p]
        if func:
            func = func[0].split('=')[1]
            prot_db[lcl] = func
        
print(len(prot_descr), len(prot_db))

db = {}
header_init = res[0]
for r in res[1:]:
    gid = r[0]
    
    while len(r) != len(header_init):
        r.append('')
    db[gid] = dict(zip(header_init, r))
    
    annot = annotate + gid + '.txt'
    if annot not in annotate_files:
        print(gid, 'Annotation file is missing!!!')
    else:
        annot = [a.split('\t') for a in open(annot).read().strip().split('\n')]
        if annot != [['']]:
            annot = [a[:-1] if 'decircle' not in a[-1] else a[:-2] + [a[-1]] for a in annot]
            annot.sort(key = lambda x: int(x[2]))
            count_annot = len(annot)
            feats = list(set([a[0] for a in annot]))
            for f in feats:
                feat = f
                if feat in resfinder_db:
                    feat = 'ResFinder'
                
                if feat not in db[gid]:
                    db[gid][feat] = []
                    
                if feat not in create and feat not in add_feats:
                    add_feats.append(feat)
                
                vals = [a for a in annot if a[0] == f]
                decircle = [v[-1].split()[-1] for v in vals if 'decircle' in v[-1]]
                remove_dups = [v for v in vals if f'{v[2]}..{v[3]}' in decircle]                
#                if remove_dups:
#                    print(gid, *remove_dups, sep = '; ')
                vals = [v for v in vals if v not in remove_dups]
                
                if feat in expand_feats:
                    for v in vals:
                        for short_name in expand_feats[feat]:
                            if short_name in v[1]:
                                v[1] = v[1].replace(short_name, expand_feats[feat][short_name])
                
                count_annot -= (len(vals) + len(remove_dups))
                if 'ref' in f:
                    vals_upd = []
                    for v in vals:
                        vid = v[1].replace('*', '').replace('#', '')
                        if vid in prot_db:
                            v = [v[0], v[1] + ' ' + prot_db[vid]] + v[2:]
                        else:
                            print(vid)
                        vals_upd.append(v)
                    vals = vals_upd[::]
                
                if feat == f:
                    vals = [a[1:] for a in vals]
                db[gid][feat] += [' '.join(a) for a in vals]
            
            if count_annot:
                print('Feat problem!', gid, count_annot)
                
    cov = coverage + gid + '.txt'
    if cov not in coverage_files:
        print(gid, 'Coverage file is missing!!!')
    else:
        cov = [a.split('\t') for a in open(cov).read().strip().split('\n')]
        status = []
        for ref in ['MF356679_D6_ref', 'AF234172_P1_ref', 'JQ965645_SSU5_ref']:
            cov_ref = [c for c in cov if ref in c and len(c)>1 and c[1]]
            if not cov_ref:
                cov_ref = '0'
            else:
                cov_ref = cov_ref[0][1]
                if float(cov_ref) >= 40:
                    status.append(ref.split('_ref')[0])
            db[gid][ref + '_cov'] = cov_ref
        db[gid]['PLP_status'] = '<<AND>>'.join(status)
        db[gid]['shape'] = '0'
        db[gid]['project_ID'] = '0'
        db[gid]['project_ID_number'] = '0'
        db[gid]['Major_replicon'] = '0'
        db[gid]['Major_replicon_variant'] = '0'
        db[gid]['Major_replicon_sequence'] = '0'
        
        

shapes = fold + 'GenBank_filtered_seqs_3Aug20.gb'
shapes = [s for s in open(shapes).read().strip().split('\n') if s[:5] == 'LOCUS']
shapes = [[ss for ss in s.split(' ') if ss] for s in shapes]
shapes = [[s[1], s[-3]] for s in shapes]

for s in shapes:
    if s[0] in db:
        db[s[0]]['shape'] = s[1]
    else:
        db[s[0]] = {}
        db[s[0]]['shape'] = 'no data'


add_feats = ['enterobacteriaceae', 'ResFinder'] + [a for a in add_feats if a not in ('enterobacteriaceae', 'ResFinder', )]
print(add_feats)      
collect = []

add_feats = ' '.join([a + '_N ' + a for a in add_feats]).split(' ')


header = create[:create.index('add_feats')] + add_feats + create[create.index('add_feats')+1:]
print(header)
print(len(create), len(add_feats), len(header), len(create)-1 + len(add_feats))

extra_feats = [h for h in header_init if h not in create]
#print(extra_feats)

reports = []
collect.append('\t'.join(header))
for r in res[1:]:
    gid = r[0]
    entry = []
    
    header_feats = [h for h in header if h[-2:] != '_N']
    for feat in header_feats:
        if feat in create and feat != 'Other_features':
            entry.append(db[gid][feat])
        elif feat in add_feats:
            if feat in db[gid]:
                vals = db[gid][feat]
                
                vals = [v.split() for v in vals]
                new_vals = []
                for v in vals:
                    if 'decircle' not in v:
                        v = ['_SPACE_'.join(v[:-5])] + v[-5:]
                    else:
                        dec_index = v.index('decircle')
                        last = '_'.join(v[dec_index:])
                        v = v[:dec_index]
                        v = ['_SPACE_'.join(v[:-5]) + '_SPACE_' + last] + v[-5:]
                    new_vals.append(v)
                
                vals = new_vals[::]
                try:
                    vals = [[aa[0], int(aa[1]), int(aa[2]), aa[3], float(aa[4]), float(aa[5])] for aa in vals]                       
                except:
                    print(feat)
                    print(vals)
                    vals = []
                vals_upd, report = find_overlaps(vals)
                vals_upd = [' '.join(map(str, v)).replace('_SPACE_', ' ') for v in vals_upd]
                num = len(vals_upd)
                vals_upd = '<<AND>>'.join(vals_upd)
                entry += [str(num), vals_upd]
                
                if report:
                    report = report.replace('_SPACE_', ' ')
                    report = f'gid\n{report}\n\n'
                    reports.append(report)
            else:
                entry += ['0', '0']
        elif feat == 'Other_features':
            others = []
            for extra in extra_feats:
                if db[gid][extra]:
                    others.append(f'{extra}: {db[gid][extra]}')
            others = '<<AND>>'.join(others)
            entry.append(others)
        else:
            print('MISSED_FEATURE', feat)
#    print(len(header), len(entry))
#    print(*entry, sep = '\n')
    
    collect.append('\t'.join(entry))
    

with open(f'{fold}Report_overlaps_{name}.txt', 'w') as fh:
    fh.write(''.join(reports))


upd_table = fold + name + '.csv'
with open(upd_table, 'w') as fh:
    fh.write('\n'.join(collect) + '\n')
    
csv_to_sqlite3(upd_table, fold + name + '.sqlite3')
os.remove(upd_table)