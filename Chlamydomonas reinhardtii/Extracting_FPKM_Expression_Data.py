### Necessary Imports
import csv
import numpy as np
import pandas as pd
from itertools import izip

### Import the raw data file as attained from the Join Genome Institute
data_from_file = pd.read_csv(r'gse25124_genes.fpkm_tracking', sep='\t')


### Dictionary to hold gene
d_genes = {}

### List of all fpkm values indices
fpkm_index_list = []
### Find all fpkm values indices
for i in range(0, len(data_from_file.columns)):
    if 'FPKM' in data_from_file.columns[i]:
        fpkm_index_list.append(i)

### Iterate through all rows and add fpkm data to dictionary
for i in data_from_file.iterrows():
    gene_identifier = i[1]['gene_id']

    for x in fpkm_index_list:
        fpkm_temp = float(i[1][x])

        if i[1]['tracking_id'] in d_genes.keys():
            d_genes[str(i[1]['tracking_id'])].append(fpkm_temp)
        else:
            d_genes[str(i[1]['tracking_id'])] = [fpkm_temp]


### Variable containing all keys within d_genes dictionary for reference
list_of_genes = d_genes.keys()

gene_names = []    # Gene IDs
gene_fpkm = []     # Gene IDs
### Iterate through dictionary and extract mean expression values for each gene
for i in list_of_genes:
    fpkm_average_temp = np.mean(d_genes[str(i)])
    gene_names.append(str(i))
    gene_fpkm.append(fpkm_average_temp)

### Save data
with open('Chlamydomonas_reinhardtii_Extracted_Expression_Data.csv', 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(izip(gene_names, gene_fpkm))
















