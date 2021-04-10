#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 16:33:21 2021

@author: lutra
"""
import sqlite3
import re

PC_lab = False

if PC_lab:
    main = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/'
else:
    main = '/data/Current_work/Phage_like_plasmids/PLP_final/'

table = 'Phage_like_plasmids_SSU5_P1_D6_12Nov20'
database = main + table + '.sqlite3'
print(database)

conn = sqlite3.connect(database)
cur = conn.cursor()

data = []
host_variety = []

missed_host = ['', 'missing', 'N/A', 'NA', 'unknown', 'Unknown', 'not available: to be reported later', 'not collected', 'not available', 'ucc isolate', 'UCC isolate',  'UCC strain', 'Missing']
bacterial_host = ['high concentration of fluoride', 'Electroporation of Y pestis Kim6+ in pCD1Ap+', 'lab strain', 'Klebsiella pneumoniae ST13-OXA48', 'Salmonella enterica subsp. enterica serovar Typhimurium strain LT2']

#modify[''] = ''
modify = {'Homo sapiens': 'human'}
modify['clinical isolate'] = 'human'
modify['biological fluid;biological fluid-human'] = 'human'
modify['Diarrheal patient'] = 'human'
modify['blood'] = 'human'
modify['human fecal sample'] = 'human'

modify['Marmota baibacina'] = 'marmot'
modify['Marmota sibirica'] = 'marmot'
modify['Urocitellus undulatus'] = 'ground squirrel'
modify['Spermophilus sp.'] = 'ground squirrel'
modify['Otospermophilus beecheyi _California ground squirrel'] = 'ground squirrel'

modify['Neopsylla setosa'] = 'flea'

modify['Sus scrofa'] = 'wild pig'
modify['Swine']='pig'
modify['swine']='pig'
modify['swine cecum']='pig'
modify['Pooled pig faecal samples collected from floor of farm']='pig'
modify['pig manure'] = 'pig'
modify['pig feces'] = 'pig'

modify['Corvus brachyrhynchos'] = 'crow'
modify['Gallus gallus domesticus'] = 'chicken'
modify['Canis lupus familiaris'] = 'dog'

modify['Bos taurus'] = 'cattle'
modify['horse gut colonizing bacterium'] = 'horse'
modify['Pooled sheep faecal samples collected from floor of farm'] = 'sheep'
modify['Pooled cattle faecal samples collected from floor of farm'] = 'cattle'

modify['beef burger'] = 'food'
modify['peanut butter'] = 'food'
modify['Pork'] = 'food'
modify['milk'] = 'food'

modify['wastewater treatment plant effluent'] = 'sewage'
modify['Humber waste water treatment plant'] = 'sewage'
modify['Hospital sewage'] = 'hospital'
modify['Shower 3'] = 'hospital'
modify['Wastewater influent sample'] = 'sewage'
modify['Wastewater effluent sample'] = 'sewage'

modify['Freshwater sample from downstream of wastewater treatment plant'] = 'environment'
modify['Freshwater sample from upstream of wastewater treatment plant'] = 'environment'

good = ['pig', 'chicken', 'dog', 'ground squirrel', 'cattle', 'lettuce', 'turkey', 'giant panda', 'environment', 'forest soil']

collect = {}
check = {}
task = 'SELECT project_ID, nucleotide, Major_replicon, organism, host, isolation_source, country, collection_date FROM ' + table
for row in cur.execute(task):
    projID, gid, major, organism, host, source, country, date = [str(r) for r in row]
    if 'SSU5' in projID:
        data.append(projID)
        if major not in collect:
            collect[major] = [[], [], [], []]
            check[major] = 0
        
        check[major] += 1
        
        organism = ' '.join(organism.split()[:2])
        collect[major][3].append(organism)
        
        year = re.findall('[0-9][0-9][0-9][0-9]', date)
        if year:
            year = year[0]
        else:
            year = 'UD'
        collect[major][2].append(year)
        
        if country == '':
            country = 'UD'
        else:
            country = country.split(':')[0]
        collect[major][1].append(country)
        
        spot = ''
        host = host.split(';')[0]
        good_host = False
        if host not in missed_host + bacterial_host:
            if host in modify:
                host = modify[host]
                good_host = True
            else:
                host = host.lower()
                if host not in good:
                    host_variety.append(host)
                else:
                    good_host = True
                    
        if not good_host:
            host = source
            if host not in missed_host + bacterial_host:
                if host in modify:
                    host = modify[host]
                    good_host = True
                else:
                    # host = host.lower()
                    if host.lower() not in good:
                        host_variety.append(host)
                    else:
                        good_host = True
        
        if not good_host:
            # host_variety.append(host)
            collect[major][0].append('UD')
        else:
            collect[major][0].append(host)
            
   
print(*set(host_variety), sep = '\n')
print(len(data))


for coll in collect:
    coll_hosts, coll_countries, coll_years, coll_orgs = collect[coll]
    gr_n = check[coll]
    print(coll, gr_n)
    for gr in [coll_hosts, coll_countries, coll_years, coll_orgs][-1:]:
        if len(gr) != gr_n:
            print('LEN PROBLEM!!!', len(gr), gr_n)
        
        gr_vars = list(set(gr))
        gr_vars = [[g, gr.count(g)] for g in gr_vars]
        if gr != coll_years:
            gr_vars = [g + [0] if g[0] != 'UD' else g + [-1] for g in gr_vars]
        else:
            gr_vars = [g + [int(g[0])] if g[0] != 'UD' else g + [0] for g in gr_vars]
        gr_vars.sort(key = lambda x: (x[2], x[1]), reverse = True)
        gr_vars = [f'{g[0]} ({g[1]})' for g in gr_vars]
        gr_vars = '\n'.join(gr_vars)
        print(gr_vars, '\n--\n\n')
        