### Required imports
import csv
import pylab
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit

### Function to be used for the inverse fit of Figures 2 A and B
def inverse_func(x, a, b):
    return a*(1/(b*(x)))

### function used to graph the best fit line for Figures 2 C, D, E, and F
def func_line(x, m, b):
    return m*x + b

### Setting font size for all graphs in this program
size_font = 16

gene_number = []         # Entrez ID
transcription_rate = []  # Transcription Rate Data
total_intron_length = [] # Intron Length Data
transcript_lengths = []  # Transcript/Exon length data

### Read in data aquired from Intron-and Exon-Length-Data-Parsing.Rmd
with open('Data_for_Graph_-_Mouse_-_unparsed.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    for row in readCSV:
        ## Do not use data if Transcription Rate is not present
        if row[4] != 'NA':
            ## Do not use data if Intron Length is not present
            if row[5] != 'NA':
                ## Do not use data if Intron Length is 0
                if row[5] != '0':
                    ## Do not use data if Gene ID is not present
                    if row[1] != 'NA':
                        gene_number.append(float(row[1]))
                        transcription_rate.append(float(row[4]))
                        total_intron_length.append(float(row[5]))
                        transcript_lengths.append(float(row[6]))


###############################################
#####Graphing Intron vs Transcription Rate#####
###########Graphing Figure 2 A and B###########
###############################################

### Produces a sequence of number to plotted with best fit
xdata = np.linspace(0, max(transcription_rate)+100, 100000)

### Produce betst fit for the inverse relationship as depicted in Figure 2 A and B
popt, pcov = curve_fit(inverse_func, transcription_rate, total_intron_length)

### Combine transcription rates data and intron length data for density operation
### and compute the density values as variable z
comb = [transcription_rate, total_intron_length]
z = gaussian_kde(comb)(comb)

### Plot Data and fits
plt.figure(1, figsize=(9, 6.75))
fig = pylab.gcf()
plt.scatter(transcription_rate, total_intron_length, c=z, s=50, edgecolors='')
plt.plot([], [], ' ', label="Inverse Fit: y=a/(b*x)")
plt.plot(xdata, inverse_func(xdata, *popt), 'r-', label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
fig.canvas.set_window_title('')
plt.title('Total Intron Length vs Transcription Rate', fontsize=size_font)
plt.ylabel('Total Intron Length (nucleotides)', fontsize=size_font)
plt.xlabel('Transcription Rate ([molecules]/(cell*hour))', fontsize=size_font)
plt.legend(loc=1, fontsize=10)
plt.grid()
plt.show(block=False)

### Suplementary
count_1 = 0
count_2 = 0
for i in range(0, len(transcription_rate)):
    if transcription_rate[i] < 7.83:
        count_1 += 1
    count_2 += 1

### Extra reference operations as can be found in the paper
print 'Number below 7.83:', count_1, 'Total number of data points:', count_2, 'Percent below 7.83', float(count_1)/float(count_2)


###############################################
#####Graphing Intron vs Transcript Length######
###########Graphing Figure 2 C and D###########
###############################################

### If desired use this code to view the data without the one outlier data point as mentioned in the text
### and replace all occurrences of variables transcript_lengths and total_intron_length in this section of code
### with intron_data_2 and script_data_2
# intron_data_2 = []
# script_data_2 = []
#
# for i in range(0, len(transcript_lengths)):
#     if transcript_lengths[i] < 40000:
#         intron_data_2.append(total_intron_length[i])
#         script_data_2.append(transcript_lengths[i])

### Produces a sequence of number to be plotted with the best fit
xdata = np.linspace(0, max(transcript_lengths), 50000)

### Produces all values for the best fit slope and intercept along with the fit's r value, p value, and standard err
slope, intercept, r_value, p_value, std_err = stats.linregress(transcript_lengths, total_intron_length)

### Combine transcript_lengths data and intron length data for density operation
### and compute the density values as variable z
comb = [transcript_lengths, total_intron_length]
z = gaussian_kde(comb)(comb)

print '#############################################'
print 'Intron vs Transcript Length Graph ### Slope:', slope, 'Intercept:', intercept, 'R-Value', r_value, 'P-Value', p_value, 'STD err', std_err
print 'Intron vs Transcript Length Graph ### R^2:', r_value**2

### Plot Data and fits
plt.figure(2, figsize=(9, 6.75))
fig = pylab.gcf()
plt.scatter(transcript_lengths, total_intron_length, c=z, s=50, edgecolors='')
plt.plot([], [], ' ', label="Line Fit: y=m*x+b")
plt.plot(xdata, func_line(xdata, slope, intercept), 'r-', label='fit: m='+str(np.round(slope, 3))+' b='+str(np.round(intercept, 3)))
plt.plot([], [], ' ', label='R^2:'+str(np.round(r_value**2, 5)))
fig.canvas.set_window_title('')
plt.title('Total Intron Length vs Total Exon Length', fontsize=size_font)
plt.legend(loc=1, fontsize=10)
plt.grid()
plt.ylabel('Total Intron Length (nucleotides)', fontsize=size_font)
plt.xlabel('Total Exon Length (nucleotides)', fontsize=size_font)
plt.show(block=True)


###############################################
#########Graphing Intron vs Predicted##########
##########Graphing Figure 2 E and F############
###############################################

### Function for predicting equation 12 in the text
def intron_func(expression, transcript, introns):
    intron_prediction = []
    actual_values = []

    for i in range(0, len(expression)):
        intron_prediction.append(transcript[i]*(1/expression[i]))
        actual_values.append(introns[i])

    return intron_prediction, actual_values

### Use prediciton function on transcription_rate,transcript_lengths, and total_intron_length]
predict, actual_values = intron_func(transcription_rate, transcript_lengths, total_intron_length)


### Produces a sequence of number to be plotted with the best fit
xdata = np.linspace(0, max(predict), 50000)

### Produces all values for the best fit slope and intercept along with the fit's r value, p value, and standard err
slope, intercept, r_value, p_value, std_err = stats.linregress(predict, actual_values)

### Combine predict data and actual_values data for density operation
### and compute the density values as variable z
comb = [predict, actual_values]
z = gaussian_kde(comb)(comb)

print '#############################################'
print 'Prediction Distribution Graph ### Slope:', slope, 'Intercept:', intercept, 'R-Value', r_value, 'P-Value', p_value, 'STD err', std_err
print 'Prediction Distribution ### R^2:', r_value**2

### Plot Data and fits
plt.figure(3, figsize=(9, 6.75))
fig = pylab.gcf()
plt.scatter(predict, actual_values, c=z, s=50, edgecolors='')
plt.plot([], [], ' ', label="Line Fit: y=m*x+b")
plt.plot(xdata, func_line(xdata, slope, intercept), 'r-', label='fit: m='+str(np.round(slope, 3))+' b='+str(np.round(intercept, 3)))
plt.plot([], [], ' ', label='R^2:'+str(np.round(r_value**2, 5)))
fig.canvas.set_window_title('')
plt.title('Total Intron Length vs Predicted', fontsize=size_font)
plt.legend(loc=1, fontsize=10)
plt.grid()
plt.ylabel('Total Intron Length (nucleotides)', fontsize=size_font)
plt.xlabel('Predicted', fontsize=size_font)
plt.show(block=True)

