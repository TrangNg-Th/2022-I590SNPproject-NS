import numpy as np


# Next, we search for SNP positions

def findSNPpos(df_clean, gene: str):
    """
    Arg : 
        df: dataframe of aligned genes 
        gene : string
        
    WARNING: 
        df has to be cleaned of missing bases 
        df can contain '-' at some position but not at the same site for all samples
    
    Output:
        dictSNP : dictionary of SNP positions for a given gene
    
    """
    
    df = df_clean[gene]

    #print()
    #print('-'*50)
    
    #print(f'Finding SNP positions for the given aligned gene {gene} sequences')
    #print()
    
    
    
    # for i, j in the samples set 
    # i takes values in [Sample1, Sample2,...]
    # j takes values in [Sample_i+1, Sample...,..]
    
    # Define reference sequence:
    for i in range(len(df.values)):
        if ('-' not in df.values[i]) & ('N' not in df.values[i]):
            refseq = df.values[i]
            refseqid = i 
            break
        
    L = [i for i in range(len(df.values)) if i != refseqid ]
    
    # List of all SNP positions
    SNPlist = []
    
    # Compare refseq to all others, at each position   
    for idpos in range(len(refseq)):  # for each position in reference sequence
        SNP = False    
        for idseq in L:
            if idpos == 0:
                print(f'Comparing sequence {refseqid}-th vs {idseq}-th')
        
         
            # Compare at the idpos if the nucleotides are different
            SNP = (refseq[idpos] != df.iloc[idseq][idpos])
            
            if SNP == True:
                if df.iloc[idseq][idpos] in ['A', 'C', 'T', 'G']:
                    if (df.iloc[idseq][idpos] not in SNPlist):
                        print(f'Detected a SNP at position : {idpos} for gene {gene}')
                        SNPlist.append(idpos)
    SNPlist = list(set(SNPlist))
    if len(SNPlist) == 0:
        print(f'No SNP detected for gene {gene}') 
    return(SNPlist)


# Given a gene 

def calculatepi(df_clean, dictSNP, gene:str):
    """
    Argument : 
        df_clean : cleaned data frame
        gene : the gene name of interest
        v: if True then display all allelic frequencies for a given position of SNP
  
    Return:
        H : dictionary of heteorozigosities given a position for the current gene
        pi
    """
    
    
    # Define variables of return
    
    H = {}
    pi = 0
    
    # Define variables to use
    N = df_clean.shape[0]  # number of samples in the data
    df = df_clean[gene]
    nbSNP = len(dictSNP[gene])  # number of SNP
    
    # Print out the number
    #print()
    #print('Examnining current gene:', gene) 
    #print()

    #print('Number of SNP found:', len(dictSNP[gene]))
    #print()
    
    
    # If there is no SNP, return nothing
    if (nbSNP  == 0):
        print('No SNP to compute pi, nothing returned')
        #print('.'*30)
        return(None)
    
    
    # If there's SNP do the following
    for pos in dictSNP[gene]:  #  iterate through all SNP positions
        nbsample = N
        d = {'A': 0, 'C': 0, 'T': 0, 'G': 0} # define allelic frequency
        #print('.'*30)
        #print(f'Position {pos} of gene {gene}')
        #print()
        
        # Pick a reference allele:
        for i in range(N):
            if ('-' not in df.iloc[i][pos]) | ('N' not in df.iloc[i][pos]):
                a_ref = df.iloc[i]
                i_ref = i
                break
                    
        L = [i for i in range(N) if i != i_ref]

                 
        for i in L: 
            #print(df.iloc[i][pos])
            a_snp = df.iloc[i][pos]
            
            # count allelic frequency in different sequences
            if a_snp in d.keys():
                d[a_snp] += 1
            else:
                d[a_snp] = 1
                nbsample -= 1
                
        # Not forget the reference allele
        d[a_ref[pos]] += 1
        
        # Formula for heterozigosity
        #d_final = d.copy()   
        #print(f'At position {pos}, the allelic frequencies are: {d_final} (NOT NORMALIZED)')
        
        d_final = {nuc: np.round(d[nuc]/nbsample, 3) for nuc in d.keys()}
        #print(f'At position {pos}, the allelic frequencies are: {d_final}')
        #print('Number of comparable samples:', nbsample)
        #print()
        H[pos] = (nbsample/(nbsample-1))*(1 - sum([d_final[nuc]**2 for nuc in d_final.keys()]))
        
    pi = sum(H.values())/len(df.iloc[0])
    #print('Pi per site for gene', gene, '=', pi)
      
    return(H, pi)



# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 17:25:26 2022

@author: ttran
"""




def codonchart(codon:str):
  
    
  """
  Given 3 nucleotides, return the protein
  """
  
  assert len(codon) == 3, 'Please give exactly 3 nucleotides'
  assert 'T' not in codon, 'Please translate T to U'
  #https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables
  
  if codon[0] == 'U':
      if codon[1] == 'U':
          if (codon[2] == 'U') | (codon[2] == 'C'):
              return('Phe')
          elif (codon[2] == 'A') | (codon[2] == 'G'):
              return('Leu') 
      
      elif codon[1] == 'C':
          return('Ser')
      
      elif codon[1] == 'A':
          if (codon[2] == 'U') | (codon[2] == 'C'):
              return('Tyr')
          elif (codon[2] == 'A') | (codon[2] == 'G'):
              return('STOP')
      
      elif codon[1] == 'G': 
          if (codon[2] == 'U') | (codon[2] == 'C'): 
              return('Cys')
          elif (codon[2] == 'A') : 
              return('STOP')
      
          elif codon[2] == 'G': 
              return('Trp')

  elif codon[0] == 'C': 
      if codon[1] == 'U': 
          return('Leu') 
      
      elif codon[1] == 'C': 
          return('Pro') 
      
      elif codon[1] == 'A': 
          if (codon[2] == 'U') | (codon[2] == 'C'): 
              return('His') 
          
          elif (codon[2] == 'A') | (codon[2] == 'G'): 
              return('Gln') 
          
      elif codon[1] == 'G':
          return('Arg')

  elif codon[0] == 'A':
      if codon[1] == 'U': 
          if (codon[2] == 'U') | (codon[2] == 'C') | (codon[2] == 'A'): 
              return('Ile') 
          elif codon[2] == 'G': 
              return('Met-START')
    
      elif codon[1] == 'C':
          return('Thr')

      elif codon[1] == 'A':
          if (codon[2] == 'U') | (codon[2] == 'C'):
              return('Asn')
          elif (codon[2] == 'A') | (codon[2] == 'G'):
              return('Lys')

      elif codon[1] == 'G':
          if (codon[2] == 'U') | (codon[2] == 'C'):
              return('Ser')
          elif (codon[2] == 'A') | (codon[2] == 'G'):
              return('Arg')

  elif codon[0] == 'G':
      if codon[1] == 'U':
          return('Val')
      elif codon[1] == 'C':
          return('Ala')
      
      elif codon[1] == 'A':
          if (codon[2] == 'U') | (codon[2] == 'C'):
              return('Asp')
      
          elif (codon[2] == 'A') | (codon[2] == 'G'):
              return('Glu')
      
      elif codon[1] == 'G':
          return('Gly')
      
        
# Function that return if SNP is a synonymousSNP
def checksynonymousSNP(dictSNP, df):
    """
    Given the positions of the SNP, 
    find out if it's a synonymous 

    Parameters
    ----------
    dictSNP : {'gene1': [SNPpos1, SNPpos2], 
               'gene2': [SNPpos1, SNPpos2],...}
        DESCRIPTION.
    
    
    df : dataframe
        
        [             gene1      |       gene2      |    gene3  
         sample1       x                   x             x 
                                         
         samplen       x                   x             x  ]

    
    
    Returns
    -------
    output : dictionary
        { 'gene1' : {SNPpos1: 0, SNPpos2: 1} ,  # 0 if synonymous, 1 if not
          'gene2' }

    """
    
    output = {}
    for gene in dictSNP:
        
        n = len(dictSNP[gene])  # nb of SNP for a given gene

        # if there's a SNP for this gene
        if n > 0:
            res = {}
            
            for pos in dictSNP[gene]: # Go through the positions of SNP for this gene
                
                
                if (pos) % 3 == 0:      
                    
                    # c is the triplet alphabet from this Sample1
                    c = (''.join(df.loc['Sample1', gene][pos: pos+3])).replace('T','U')
                    amino_ref =  codonchart(c)  # Translate c to amino acid
                    
                    test = False  # a test param
                    
          
                    for s in list(df[gene].index)[1:]:
                        
                        # Extract from another sample2, sample3,...
                        c2 = (''.join(df.loc[s, gene][pos: pos + 3])).replace('T','U')
                        amino_cmp = codonchart(c2)
                        
                        if amino_ref != amino_cmp :
                            test = True
                            break
                
                elif (pos) % 3 == 1:
                    c = (''.join(df.loc['Sample1', gene][pos - 1 : pos + 2])).replace('T','U')
                    amino_ref = codonchart(c)
                    test = False  # a test param
                    
          
                    for s in list(df[gene].index)[1:]:
                        c2 = (''.join(df.loc[s, gene][pos - 1 : pos + 2])).replace('T','U')
                        amino_cmp = codonchart(c2)
                        if amino_ref != amino_cmp:
                            test = True
                            break
                
                elif (pos) % 3 == 2:
                    c = (''.join(df.loc['Sample1', gene][pos - 2 : pos + 1])).replace('T','U')
                    amino_ref = codonchart(c)
                    test = False  # a test param
                    
                    for s in list(df[gene].index)[1:]:
                        c2 = (''.join(df.loc[s, gene][pos - 2 : pos + 1])).replace('T','U')
                        amino_cmp = codonchart(c2)
                        if amino_ref != amino_cmp :
                            test = True
                            break
                        
                res[pos] = test
            output[gene] = res
    return(output)
          
        
        

          
    
def harmonic(n):
    """
    return harmonic number until from 1 to rank n
    """
    an = 0
    for i in range(1, n+1):
        an += 1/i
        
    return(an)

    

# Watterson theta
def WattersonTheta(df_clean, dictSNP, gene):
    
    """
    Arg:
        df_clean : clean dataframe
        gene : gene name
        dfSNP : dictionary of all segregating sites for all genes in df_clean
    Return
        Wt : Watterson theta gor a given gene
    """
    # n:number of samples
    n = df_clean[gene].shape[0]
    
    # K: Number of segregating sites
    K = len(dictSNP[gene])
    
    # Wt: watterson theta
    Wt = K/harmonic(n-1)
    
    return(Wt)
    