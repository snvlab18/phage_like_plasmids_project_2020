#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 17:41:18 2019

@author: lutra
"""
import sqlite3

def csv_to_sqlite3 (file_in, database, new_line = ''):
    conn=sqlite3.connect(database)
    cur=conn.cursor()
    
    table = (((file_in.split('/'))[-1]).split('.'))[0]
    table = table.replace('-', '_')
    print (table)
    ex1='DROP TABLE IF EXISTS ' + table
    cur.execute(ex1)
    conn.commit()
    
    work = ((open(file_in)).read()).split('\n')
    work = [w.replace('|','_') for w in work]
    work = [w.replace('"','_') for w in work]
    work = [w.replace("'",'_') for w in work]
    work = [w.replace('(','_') for w in work]
    work = [w.replace(")",'') for w in work]
    work = [w.replace(",",'') for w in work]
    
    work = [w.split('\t') for w in work if w]
    
    check_size = max([max([len(ww) for ww in w]) for w in work])
    print (check_size)
    field_limit = check_size+1
    
    header = work[0]
    header = [h.replace(' ','_') for h in header]
    header = [h.replace('-','_') for h in header]
    header = [h.replace('+','') for h in header]
    header = [h.replace('[','') for h in header]
    header = [h.replace(']','') for h in header]
    header = [h.replace('.','') for h in header]
    header = [h.replace('/','_') for h in header]
    header = [h if h != 'group' else 'group_whatever' for h in header]
    print (*header,'\n')
    
    column = 'id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE'
    insert = []
    for i in range(len(header)):
        check_type = work[1]
        insert.append(header[i])
        if '.' in check_type[i]:
            try:
                check = float(check_type[i])
                column += (', ' + header[i] + ' REAL')
            except:
                column += (', ' + header[i] + ' TEXT')
                if 'sequence' in header[i]:
                    column += ('({})'.format(field_limit))
        else:
            try:
                check = int(check_type[i])
                column += (', ' + header[i] + ' INTEGER')
            except:
                column += (', ' + header[i] + ' TEXT')
                if 'sequence' in header[i]:
                    column += ('({})'.format(field_limit))
    
    ex2='CREATE TABLE {} ({})'.format(table, column)
#    print(ex2)
#    print(*ex2.split('group'), sep ='\n\n')
    cur.execute(ex2)
    conn.commit()
    
    print (len(work[1:]), end = ', ')
    count = 0
    
    insert = ', '.join(insert)
    for w in work[1:]:
        w = ['"'+ww.replace('<<AND>>', '\n')+'"' for ww in w]  
        
        if new_line:
            w = [ww.replace(new_line, '\n') for ww in w]
        
        w = ', '.join(w)
        task = 'INSERT INTO {} ({}) VALUES ({})'.format(table, insert, w)
#        print(task)
        cur.execute(task)
        count +=1
    print(count)
    conn.commit()
    conn.close()
    print('\n\n')
    return

