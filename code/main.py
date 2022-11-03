# Add the directory to find the codes

# Import some function
import os

# Set path
path = 'C:\\Users\\ttran\\OneDrive - Indiana University\\2022-I590-SNP&PopulationGenetics-project\\Codes\\'
os.chdir(path)


import pi
import prepfile as prep  # import the file prepfile.py 

#import watterson

# In GoogleDrive
#path = '/content/drive/MyDrive/SNP project background/Alignment example/'
#filename =  'Alignment_format_example.TXT'

# -----------------------------------------------
# Prepare the input file 
path = os.getcwd()  # in local computer
filename =  'Alignment_format_example.TXT'
df = prep.prepfile(filename=filename)

# Check the dataframe
print(df.head())

# Check the SNP positions
# dictSNP contains SNP positions for each gene
dictSNP = {gene: pi.findSNPpos(df, gene) for gene in df.columns}

# Calculate pi#
print('-'*50)
print()
print('Calculating heterozigosities and pi ...')
print()
print('-'*50)


for g in dictSNP:
    res = pi.calculatepi(df, dictSNP, gene=g)
    if res is not None:
        H, p = res[0], res[1]
        print('Pi per site for gene', g, '=', p)

print('-'*50)     
# Check synonymous SNP 
d = pi.checksynonymousSNP(dictSNP, df)

#print('-'*50)  
# Check Watterson theta
for g in dictSNP:
    Wt = pi.WattersonTheta(df, dictSNP, gene=g)
    print('Watterson theta for gene', g, '=', Wt)
