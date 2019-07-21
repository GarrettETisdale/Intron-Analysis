### Necessary Imports
import csv
import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde



gene_number_temp_1 = []               ### Gene names provided by the annotation data
total_intron_length_temp_1 = []       ### Intron lengths as calculated with the annotation data
exon_lengths_temp_1 = []              ### Exon lengths as caclulated with the annotation data (not used)
### Iterate through extracted intron data output from Extracting_Intron_Length_Data.py and add to list
with open('Dictyostelium_discoideum_Extracted_Intron_Data.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    for row in readCSV:
        gene_number_temp_1.append(row[0])
        total_intron_length_temp_1.append(float(row[1]))
        exon_lengths_temp_1.append(float(row[2]))


gene_number_temp_2 = []   ### Gene names as provided by the expression data
expression_temp_2 = []    ### Expression data
### Iterate through extracted expression data output from Extracting_FPKM_Expression_Data.py and add to list
with open('Dictyostelium_discoideum_Extracted_Expression_Data.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    for row in readCSV:
        gene_number_temp_2.append(row[0])
        expression_temp_2.append(float(row[17]))  ### Take the mean value


gene_number = []       ### Gene IDs after matching intron and expression data
expression = []        ### Expression data after matching intron and expression data
intron_length = []     ### Intron lengths after matching intron and expression data
exon_length = []       ### Exon lengths after matching intron and expression data
### Matching IDs and data from intron data file and expression data file
for i in range(0, len(gene_number_temp_1)):
    if gene_number_temp_1[i] in gene_number_temp_2:
        ### If intron length is greater than 0
        if total_intron_length_temp_1[i] > 0:
            gene_number_index = gene_number_temp_2.index(gene_number_temp_1[i])  ### Find index of desired gene in other data set
            ### If expression is greater than 0
            if expression_temp_2[gene_number_index] > 0:
                gene_number.append(gene_number_temp_1[i])
                intron_length.append(total_intron_length_temp_1[i])
                exon_length.append(exon_lengths_temp_1[i])
                expression.append(expression_temp_2[gene_number_index])

### Setting font size for all graphs in this program
size_font = 16

### Combine control data and intron length data for density operation
### and compute the density values as variable z
comb = [expression, intron_length]
z = gaussian_kde(comb)(comb)

### Plot Data
plt.figure(6, figsize=(9, 6.75))
fig = pylab.gcf()
plt.scatter(expression, intron_length, c=z, s=50, edgecolors='')
fig.canvas.set_window_title('')
plt.title('Dictyostelium discoideum - Total Intron Length vs Gene Expression', fontsize=size_font)
plt.ylabel('Total Intron Length (nucleotides)', fontsize=size_font)
plt.xlabel('Expression (Averaged FPKM)', fontsize=size_font)
plt.grid()
plt.show(block=True)

