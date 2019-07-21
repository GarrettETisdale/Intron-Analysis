### Required imports
import csv
import pylab
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

gene_number = []             # GI number
expression_in_exp_form = []  # Gene expression
total_intron_length = []     # Total intron length
transcript_lengths = []      # Transcript/exon length

### Read in data aquired from Cristian Castillo-Davis
with open('Human_-_Data_-_Castillo-Davis.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    for row in readCSV:
        ## Do not use data if Expression is 0
        if row[1] != '0':
            ## Do not use data if Intron Length is 0
            if row[3] != '0':
                gene_number.append(row[0])
                expression_in_exp_form.append(float(row[1]))
                total_intron_length.append(float(row[3]))
                transcript_lengths.append(float(row[4]))

### Setting font size for all graphs in this program
size_font = 16

### Combine control data and intron length data for density operation
### and compute the density values as variable z
comb = [expression_in_exp_form, total_intron_length]
z = gaussian_kde(comb)(comb)

### Plot Data
plt.figure(5, figsize=(9, 6.75))
fig = pylab.gcf()
plt.scatter(expression_in_exp_form, total_intron_length, c=z, s=50, edgecolors='')
fig.canvas.set_window_title('')
plt.title('Human (Normal) - Total Intron Length vs Gene Expression', fontsize=size_font)
plt.ylabel('Total Intron Length (nucleotides)', fontsize=size_font)
plt.xlabel('Expression (number of expression sequence tag hits)', fontsize=size_font)
plt.grid()
plt.show(block=True)
