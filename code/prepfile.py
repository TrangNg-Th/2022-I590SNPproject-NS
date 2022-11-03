import numpy as np
import pandas as pd


def prepfile(filename, path=''):
    # Add all the content from alignment to a list of file
    
    rawfile = []
    idxbgblock, idxendblock = 0, 0 # index of beginblock and endblock
    idxbgset, idxendset = 0, 0 # index of beginset and endset

    
    print(f'The input file is: {filename}')
    print()
    
    with open(path + filename, "r") as f:
        i = 0
        for l in f:
            rawfile.append(l)
            if "beginblock" in l:
                idxbgblock = i+1
            elif "endblock" in l:
                idxendblock = i
            elif "beginset" in l:
                idxbgset = i+1
            elif "endset" in l:
                idxendset = i
            i += 1
            
            
    # extract the sequence sample
    stitles, scontents = [], []  # sample title, sample content
    
    for j in range(idxbgblock, idxendblock):
        s_j = rawfile[j].split()
        if len(s_j) > 1:  # only take the line that has a title and a seq (for now)
            stitles.append(s_j[0])
            scontents.append(s_j[1])

    # Extract the tiltes of the genes and the position of the genes
    gtitles, gcontents = [], []  # gene title, gene content
    for j in range(idxbgset, idxendset):
        g_j = rawfile[j].split()  
        if len(g_j) > 1:  # extract gene position ex:'psbD:1-1062'
            gname = g_j[1].split(':')[0]  # extract gene name
            gpos = g_j[1].split(':')[1].split('-')
            
            gtitles.append(gname)
            gcontents.append([int(gpos[0]), int(gpos[1])])

  
            
    # Prepare dictionary 
    # d will be of the following shape:
    # {gene1 : [sample1, sample2, sample3,...],
    #  gene2 : [sample1, sample2, sample3,...],...}
    # where sample1 : a string of nucleotides
    d = {gtitles[k]: [] for k in range(len(gtitles))}
    for j in range(len(stitles)):
        for k in range(len(gtitles)):
            s = scontents[j][gcontents[k][0] - 1: gcontents[k][1]]
            
            # Convert the sequence to an array --> Not really necessary I think
            #s_array = np.array(list(s))
            
            s_array = s
            d[gtitles[k]].append(s_array)

    print('Number of genes (before removing missing bases)=', len(d))
    print('Number of sample=', len(list(d.values())[0]))
    
    
    # Cleaning the dictionary d
    # d will be of the following shape:
    # {gene1 : [sample1, sample2, sample3,...],
    #  gene2 : [sample1, sample2, sample3,...],...}
    # Remove the positions at which all the samples are of bad quality 'N' or '-'
    
    
    # d_rm: contains all positions to be removed of each gene
    d_rm = {gene: [] for gene in d.keys()}  
    
    
    
    for gene in d.keys():  # check for each gene
        
        genelen = len(d[gene][0])  # Ask Alex if all sample is of same length ?
        
        for i in range(genelen):  # for each position
            # Since all gene starts with ATG
            # The first position is supposed to be an 'A' avery where
            
            rmpos = 0  # a count to check if at the current position, 
                      # every sample contains "N" or "-"
            
            for s in d[gene]:  # for each sequence sample1, sample2, sample3,..
                if s[i] not in ['A', 'T', 'G', 'C']:
                    rmpos += 1  # add 1 if this position is bad                   

            if rmpos == len(d[gene]) : # if all position is bad, remove this base
                d_rm[gene].append(i)
        d_rm[gene].reverse()
 
    
    # Go back to dictionary d, clean out all po
    d_clean = {k: d[k].copy() for k in d}
    for gene in d.keys():
        
        if len(d_rm[gene]) > 0: # If there's any position to be removed
        
            for i in d_rm[gene]:  # go through index to be removed
                d_clean[gene] = [sample[:i] + sample[i+1:] for sample in d[gene]]
            
    # Create data frame
    df = pd.DataFrame(d_clean, columns=d.keys(), index=stitles)

    return(df)

    
    


import os
path = os.getcwd()  # in local computer
filename =  'Alignment_format_example.TXT'
df = prepfile(filename=filename)