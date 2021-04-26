1) Search with BLASTn
- run blastn with three extrachromosomal lysogenic phages
megablast, 20000 aligned seqs
download all results as text

- filter for seqs covering the phage for >=20%, download these results as text -> **_cov20plus_**
for file in *cov20plus_*; do grep -v '#' $file | cut -f2 | sort -u; done | sort | uniq -c | awk '{if ($2) print $2}' > Entrez_cov20plus.txt

1_efetch_nuccore_by_biosample.py
    - collect fasta files and metadata



2) Annotate fasta files
- prepare folders db/ and coverage_db/

Download databases from ResFinder, PlasmidFinder (enterobacteriaceae.fsa, save other replicons into all_other_replicons.fsa), oriTDB, VFDB. 
All IS elements found in the reference phages (D6, P1, SSU5) by the ISfinder online.
Get reference CDSs from the reference GenBank pages.

!!! Replicons will be expanded during the second round of the annotation. See below!

DATABASES:
db/
├── AF234172_P1_ref_CDS.fsa
├── all_other_replicons.fsa
├── D6_putative_replicon_orf42.fsa
├── enterobacteriaceae.fsa
├── ISel_db.fsa
├── JQ965645_SSU5_ref_CDS.fsa
├── MF356679_D6_ref_CDS.fsa
├── oriTDB
│   ├── auxiliary_all_blastx.fsa
│   ├── oriT_all.fsa
│   ├── relaxase_all_blastx.fsa
│   └── t4cp_all_blastx.fsa
├── resfinder_db
│   ├── aminoglycoside.fsa
│   ├── antibiotic_classes.txt
│   ├── beta-lactam.fsa
│   ├── CHECK-entries.sh
│   ├── colistin.fsa
│   ├── config
│   ├── fosfomycin.fsa
│   ├── fusidicacid.fsa
│   ├── glycopeptide.fsa
│   ├── INSTALL.py
│   ├── macrolide.fsa
│   ├── nitroimidazole.fsa
│   ├── notes.txt
│   ├── oxazolidinone.fsa
│   ├── phenicol.fsa
│   ├── phenotype_panels.txt
│   ├── phenotypes.txt
│   ├── quinolone.fsa
│   ├── README.md
│   ├── rifampicin.fsa
│   ├── sulphonamide.fsa
│   ├── tetracycline.fsa
│   └── trimethoprim.fsa
└── VFDB_setB_nt.fsa

COVERAGE REFS
coverage_db/
├── AF234172_P1_ref.fasta
├── JQ965645_SSU5_ref.fasta
└── MF356679_D6_ref.fasta

- first round of the annotation
2_1_main_annotate.py

- collect metadata (3_1_process_metadata.py), divide seqs into groups (7_1_divide_groups_into_clusters.py)  and return
- using lists of grouped sequences IDs, collect all relevand GenBank files via NCBI Batch Entrez

- collect all annotated proteins using prokka-toolkit to use for prokka annotation

prokka-genbank_to_fasta_db D6_genbanks.gb P1_genbanks.gb SSU5_genbanks.gb ARG_extra_Incs_matched_plasmids.gb > PLP_genbanks.faa
cd-hit -i PLP_genbanks.faa -o PLP_genbanks -T 0 -M 0 -g 1 -s 0.8 -c 0.9

- update PLP_genbanks with improved proteins names
2_2_cd_hit_wrapper.py
    - modify names (shorten, unify)
    - add labels to color proteins with GenePlotR later on
    + generate table with all info about proteins in clusters
    
- fake Backteria kingdom in Prokka with a new annotation file
cp PLP_genbanks_upd /path/to/prokka/db/kingdom/Bacteria_PLP/sprot
makeblastdb -hash_index -dbtype prot -in /path/to/prokka/db/kingdom/Bacteria_PLP/sprot -logfile /dev/null
cp /path/to/prokka/db/kingdom/Bacteria_PLP/* /path/to/prokka/db/kingdom/Bacteria/

- run Prokka
2_3_prokka.py

- based on the literature data extract D6 replicon in a new file D6_putative_replicon_orf42.fsa
- based on the Porkka annotation extract SSU5 replicons and add them to the file enterobacteriaceae.fsa
2_4_extract_dna.py

run the second round of the annotation
2_1_main_annotate.py



3) Combine metadata into an SQL database
3_1_process_metadata.py
3_2_process_metadata_add_project_IDs.py

Collect lists of grouped PLPs and return to step 2
Then rerun step 3 and continue

3_3_process_metadata_add_major_replicon_variant.py



4) Prepare tables about groups
Supplementary Table S1
4_1_PLP_prepare_TableS1.py

Summary for Table 1
4_2_prepare_table_about_major_replicons.py

Supplementary Table S3
4_3_TableS3_lists_of_ARGs.py
    - remove overlapping ARGs from the report
 
Сolor the Table S3 manually



5) Figure 2. Circos diagrams
5_circos_prophages.py
    - intermediate table for all groups
    - special table with Circos format
    !! get svg images http://mkweb.bcgsc.ca/tableviewer/visualize/
    - remove repeated ARGs
    - upd svg
    - add legend
    - make legend for GenePlotR figures too
    

    
6) Chi square test with Yates correction
6_PLP_chi_square_test.py
    - ARGs are associated with minor replicons
    - Y.pestis is associated with virulence genes

    
    
7) Divide sequences into groups, shift seqs to the first nt of the reference, prepare GenBank files for the vizualisation
7_1_divide_groups_into_clusters.py
    - init group folders
    - divide seqs into groups
    - shift seqs to the first nt of the reference
    - prepare iTol documents (not used in the publication)
    
7_2_Easyfig_master_v2.py
    - combine main and Prokka annotation into GenBank file for each grouped sequence
    - run Easyfig (not used in the publication)

7_3_geneplotr_with_any_set_of_sequences.py
    - specify a subgroup of any PLP group
    - prepare an input files for gepard
    - make GenePlotR R scripts


    
8) ARG-encoding regions from plasmids with additional replicons
mkdir ARG_maps
mkdir ARG_maps/cut_ARG_extra_Incs

8_1_ARG_surroundings.py
    - make list of ARG-encoding plasmids with additional replicons (AEPWAR)
    - check the file how each AEPWAR aligns to the reference
    - find coordinates of inserted blocks which encode both ARGs and additional replicons between regions aligned to the reference
    - cut these inserted blocks and save them separately
    
BLASTn the separated inserted blocks (megablast, 100 aligned seqs)
choose the best matches
collect fasta and genbank (ARG_extra_Incs_matched_plasmids.gb) files
annotate plasmids returned the best matches with the separated inserted blocks (Step 2. Annotate. Note that Prokka database is created including ARG_extra_Incs_matched_plasmids.gb)
    
8_2_rotate_ARG_extra_Inc_around_its_insert.py
    - shift plasmids returned the best matches with the separated inserted blocks so they would have the best visual alignment with ARG-encoding PLPs with additional replicons
    !! (run 7_2_Easyfig_master_v2.py to prepare GenBank files)
    - make GenePlotR R scripts
    
    
    
7+8) Vizualise plasmid annotations and comparisons
Run GenePlotR R scripts, get svg images
Process svg images, add legend (5_circos_prophages.py), save jpg

Run Gepard
cd /path/to/gepart_input_files
ls | cut -d. -f1 > gepard_task.txt
suff_size=1500

while read suff; do /data/Bioinformatics/gepard/gepard-1.30/gepardcmd.sh -seq1 $suff.fasta -seq2 $suff.fasta -matrix /data/Bioinformatics/gepard/gepard-1.30/matrices/edna.mat -word 15 -maxwidth $suff_size -maxheight $suff_size -outfile $suff.png; done < gepard_task.txt



9) Distribution of PLP-related replicons
Collect PLP major replicons as separate files
BLASTn fasta files (blastn, 20000 aligned seqs)
filter results for >60% of query coverage, >90% of identity (these sequences carry PLP-related replicons)

get list of GenBank IDs of plasmids and chromosomes carrying PLP-related replicons
efetch fasta files and metadata for these sequences if they are absent in our PLP database
save fasta files with all original fasta files in these project, annotate newly added sequences (as in the Step 2)

expand the working database (as in the Step 3)

9_1_distribution_of_PLP_replicons.py
    - subcluster all sequences with the PLP-related relicons per group (grouped/ungrouped, complete/not, chromosome/not)
    - draw plots
    - save svg
    
Process svg, save jpg

Prepare Supplementary Table S2
9_2_PLP_prepare_TableS2.py



10) Genetic comparison and clustering within each phage group (D6, P1 and SSU5) based on the distribution of the reference CDSs
mkdir Heatmap_CDSs in each PLP group dir
save the reference CDSs into a special folder (i.e., /path/to/D6/Heatmap_CDSs/MF356679_CDSs/MF356679_CDSs.fasta)

10_CDS_matrix_blastn_muscle_any_seq.py
    - run blastn locally for all grouped sequences against the reference CDSs
    - collect blastn results into a heatmap 
    - draw heatmap, add dendrograms, ARGs lists, legens
    - save svg
    
Process svg, save jpg
    


