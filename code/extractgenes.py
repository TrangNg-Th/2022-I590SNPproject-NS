# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
from Bio import SeqIO


prj_path = 'C:\\Users\\ttran\\OneDrive - Indiana University\\2022-I590-SNP&PopulationGenetics-project\\'
code_path = prj_path + 'Git\\code\\'
data_path = prj_path + 'data\\Oryza-sativa\\'
output_path = data_path + 'Indica\\full-genome\\'


os.chdir(code_path)

#genome = 'Indica\\full-genome\\mtDNA-fullgenome_AP011077.1-Oryza-sativa-Indica-Group-mitochondrial-DNA2C-cultivar_-Lead-rice_Gen.gb'
#genome = 'Indica\\full-genome\\mtDNA-fullgenome_CP018169.1-Oryza-sativa-Indica-Group-cultivar-Shuhui498-mitochondrion%2C-complete-genome_GenBank.gb'
genome = 'Indica\\full-genome\\OS-zenshan97-cpDNA_CP056064.1-Oryza-sativa-Indica-Group-cultivar-Zhenshan-97-chloroplast%2C-complete-genome_GenBank.gb'

# =============================================================================
# # import gene list
# with open(gnamecpDNA, 'r') as input_file:
#     gene_names = [line.strip('\n') for line in input_file]
#     
# =============================================================================



# open file for genome
gb_object = SeqIO.read(data_path + genome, 'gb')
all_CDS = []


for feature in gb_object.features :
    if feature.type == 'CDS':
        start = feature.location.start
        end = feature.location.end
        strand = feature.location.strand
        all_CDS.append([feature, start, end, strand])
            
            
gene_sequences = []
gene_sequences1 = []
for seq in all_CDS:
    gene = seq[0]  

    if 'gene' in gene.qualifiers.keys():
        gene_name = gene.qualifiers['gene'][0]

        if '-fragment' in gene_name:
            pass
        else:
            if seq[3] == -1:
                #print(gene.extract(gb_object))
                extract = gene.extract(gb_object).reverse_complement()
                extract1 = gene.extract(gb_object).reverse_complement()
                #print(extract)
                #print('='*50)
            else:
                extract = gene.extract(gb_object)
                extract1 = gene.extract(gb_object)
                
            #extract.seq = gene.extract(gb_object).seq
            extract.id = gene_name
            extract1.id = gene_name
            extract.description = ''
            extract1.description = '(' + str(seq[1]) +':' + str(seq[2]) + ')' + ', strand=' + str(seq[3])
            gene_sequences.append(extract)
            gene_sequences1.append(extract1)
            #print('gene %s has been found'%gene_name)

             
# =============================================================================
# outputfile = output_path + 'extracted-cds\\OS-cpDNA-Zhenshan-97-CDS.fasta'
# SeqIO.write(gene_sequences, outputfile, format ='fasta')   
# 
# outputfile = output_path + 'extracted-cds\\OS-cpDNA-Zhenshan-97-CDS-full.fasta'
# SeqIO.write(gene_sequences1, outputfile, format ='fasta')    
# =============================================================================



# =============================================================================


