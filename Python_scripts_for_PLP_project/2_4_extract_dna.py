#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:52:49 2019

@author: lutra
"""
def rev_comp(dna):
    trantab = "".maketrans("ACTGactg","TGACtgac")
    return dna[::-1].translate(trantab)


def find_reg(file_in, c1, c2):
    dna=''.join((((open(file_in).read()).strip()).split('\n'))[1:])
    piece=dna[c1-1:c2]
    return piece

def fold_dna(dna, size = 60):
    dna = [dna[i*size:(i+1)*size] for i in range(int(len(dna)/size)+1)]
    dna=[d for d in dna if d]
    dna = '\n'.join(dna)
    return dna

def dna_to_fasta(fasta_file, name, dna):
    dna = [dna[i*60:(i+1)*60] for i in range(int(len(dna)/60)+1)]
    dna=[d for d in dna if d]
    dna = '\n'.join(dna)
    with (open(fasta_file, 'w')) as fh:
        fh.write('>{}\n{}\n'.format(name, dna))
    return


file_in = '/home/shiri/kpn_project/U95/Houston_run/pU95.fasta'
c1, c2 = 1, 16684

file_in = '/home/shiri/plasmid_project/Phage_like_plasmids/PLP_final/coverage_db/MF356679_D6_ref.fasta'
c1, c2 = 45699, 46596

file_in = '/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/D6/HomBlocks_per_cluser/other/fasta/CP056699_D6_c49378rev.fasta'
c1, c2 = 49373, 50514 #D6_CP056699_repE_plus_upstream
c1, c2 = 49759, 50514 #D6_CP056699_repE_only
c1, c2 = 49373, 49758 #D6_CP056699_upstream_repE_only

file_in = '/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/D6/HomBlocks_per_cluser/other/fasta/CP042620_D6_c42481rev.fasta'
c1, c2 = 46352, 47613 #D6_CP042620_repB_plus_upstream
c1, c2 = 46729, 47613 #D6_CP042620_repB_only
c1, c2 = 46352, 46728 #D6_CP042620_upstream_repB_only

file_in = '/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/D6/HomBlocks_per_cluser/other/fasta/AP018804_D6_c7937rev.fasta'
c1, c2, name = 44240, 45454, 'D6_AP018804_repFIB_plus_upstream'
# c1, c2, name = 44594, 45454, 'D6_AP018804_repFIB_only'
# c1, c2 = 44240, 44593 #D6_AP018804_upstream_repFIB_only



#collect SSU5 other replicons
file_in = ''
c1, c2, name, rev = 1, 2, '', False

# file_in = '/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/SSU5/HomBlocks_per_cluser/other/fasta/CP057631_SSU5_c104093.fasta'
# c1, c2, name, rev =  65186, 66259, 'repA_Citro', True


# file_in = '/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/SSU5/HomBlocks_per_cluser/other/fasta/CP030187_SSU5_c66011rev.fasta'
# c1, c2, name, rev = 64956, 66011, 'repA_FIB_Sent', True


# file_in = '/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/SSU5/HomBlocks_per_cluser/other/fasta/CP012254_SSU5_c98290.fasta'
# c1, c2, name, rev = 63719, 64819, 'repA_Crono', True


# file_in = '/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/SSU5/HomBlocks_per_cluser/other/fasta/CP014707_SSU5_c63071rev.fasta'
# c1, c2, name, rev = 61854, 63071, 'repA_Sent', True


# file_in = '/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/SSU5/HomBlocks_per_cluser/other/fasta/KY515226_SSU5_c80639.fasta'
# c1, c2, name, rev = 62006, 63223, 'repA_Sent_1', True


# file_in = ''
# c1, c2, name, rev = 1, 2, 'repA_Sent_2', False


# file_in = ''
# c1, c2, name, rev = 1, 2, 'repA_Sent_3', False


fold_out = '/data/Current_work/Phage_like_plasmids/PLP_final/prophages_groups/SSU5/HomBlocks_per_cluser/other/extract_repA_replication_proteins/'

name =  name + ' ' + file_in.split('/')[-1][:-6]
print(name)
fasta_file = fold_out +  name.replace(' ', '_') + '.fasta'

piece = find_reg(file_in, c1, c2)
if rev:
    piece = rev_comp(piece)

fin = fold_dna(piece)
dna_to_fasta(fasta_file, name, piece)
#print (fin)
