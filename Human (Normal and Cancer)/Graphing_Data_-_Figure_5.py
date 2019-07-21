### Required imports
import csv
import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

gene_number = []                       # Entrez ID
expression_in_exp_form_control = []    # Control Expression
expression_in_exp_form_cancer = []     # Cancer Expression
total_intron_length = []               # Total intron length
transcript_lengths = []                # Total transcript/exon length

### Read in data aquired from Intron-and Exon-Length-Data-Parsing.Rmd
with open('Data_for_Graph_-_Human_(Control-and-Cancer)_-_unparsed.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    for row in readCSV:
        ## Do not use data if Control Expression is 0
        if row[2] != '0':
            ## Do not use data if Cancer Expression is 0
            if row[3] != '0':
                ## Do not use data if Intron Length is 0
                if row[4] != '0':
                    ## Do not use data if Gene ID is not present
                    if row[1] != 'NA':
                        gene_number.append(row[0])
                        expression_in_exp_form_control.append(float(row[2]))
                        expression_in_exp_form_cancer.append(float(row[3]))
                        total_intron_length.append(float(row[4]))
                        transcript_lengths.append(float(row[5]))

### Setting font size for all graphs in this program
size_font = 16


###############################################
#####Graphing Intron vs Control Expression#####
##############Graphing Figure 5 A##############
###############################################

### Combine control data and intron length data for density operation
### and compute the density values as variable z
comb = [expression_in_exp_form_control, total_intron_length]
z = gaussian_kde(comb)(comb)

### Plot Data
plt.figure(5, figsize=(9, 6.75))
fig = pylab.gcf()
plt.scatter(expression_in_exp_form_control, total_intron_length, c=z, s=50, edgecolors='')
fig.canvas.set_window_title('')
plt.title('Human (Normal) - Total Intron Length vs Gene Expression', fontsize=size_font)
plt.ylabel('Total Intron Length (nucleotides)', fontsize=size_font)
plt.xlabel('Expression (Intensity)', fontsize=size_font)
plt.grid()
plt.show(block=False)


###############################################
#####Graphing Intron vs Cancer Expression######
##############Graphing Figure 5 B##############
###############################################

### Combine cancer data and intron length data for density operation
### and compute the density values as variable z
comb = [expression_in_exp_form_cancer, total_intron_length]
z = gaussian_kde(comb)(comb)

### Plot Data
plt.figure(6, figsize=(9, 6.75))
fig = pylab.gcf()
plt.scatter(expression_in_exp_form_cancer, total_intron_length, c=z, s=50, edgecolors='')
fig.canvas.set_window_title('')
plt.title('Human (Cancer) - Total Intron Length vs Gene Expression', fontsize=size_font)
plt.ylabel('Total Intron Length (nucleotides)', fontsize=size_font)
plt.xlabel('Expression (Intensity)', fontsize=size_font)
plt.grid()
plt.show(block=True)
