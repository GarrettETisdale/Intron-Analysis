### Required imports
import csv
import pylab
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit

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

### Reproduction of Figure 2 E and F for transparency when data is carried over into Figure 3
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
plt.show(block=False)

###############################################
######Graphing Intron vs Binned Predicted######
##########Graphing Figure 3 A and B############
###############################################

### determining threshold limit values as mentioned in text
count = 0
for i in range(0, len(predict)):
    if predict[i] < 5000:
        count += 1

print 'Total number of data point:', len(predict), 'Number of values below 5,000 limit threshold:', count, 'Percent below 5,000 limit threshold:', float(count)/float(len(predict))

### The following code segregates the data into bins and saves all of the data within the bin into a
### dictonary to be called and averaged latter.
######Bin Code - Start######
d_bin = {}
bin_count = 0
bin_limit_1 = 0
bin_limit_2 = 500
delta = 100 # 100 & 500   ### Change this variable to produce either the 100 bin range or 500 bin range

flat_limit = 5000
while bin_limit_2 < max(predict):
    bin_temp = []
    bin_range_temp = []

    for i in range(0, len(predict)):
        if predict[i] >= bin_limit_1:
            if predict[i] < bin_limit_2:
                if predict[i] < flat_limit:
                    bin_temp.append(actual_values[i])
                    bin_range_temp.append(predict[i])

    d_bin["Bin {0}".format(bin_count)] = bin_temp
    d_bin["Bin {0} - Range".format(bin_count)] = bin_range_temp

    bin_limit_1 += delta
    bin_limit_2 += delta
    bin_count += 1
######Bin Code - End######

### The following code takes the mean values win the binned regions
######Bin Average - Start######
mean_bin_values = []    ### Contains mean values on the actual value axis
mean_bin_ranges = []    ### Contains mean values on the predicted axis
for i in range(0, len(d_bin.keys())/2):
    bin_temp = d_bin["Bin {0}".format(i)]
    bin_range = d_bin["Bin {0} - Range".format(i)]

    if len(bin_temp) > 0:
        mean_bin_values.append(float(np.mean(bin_temp)))
        mean_bin_ranges.append(float(np.mean(bin_range)))
    else:
        pass
######Bin Average - End######

### Produces a sequence of number to be plotted with the best fit
xdata = np.linspace(0, max(mean_bin_ranges), 50000)

### Produces all values for the best fit slope and intercept along with the fit's r value, p value, and standard err
slope, intercept, r_value, p_value, std_err = stats.linregress(mean_bin_ranges, mean_bin_values)

print '#############################################'
print 'Binned Graph ### Slope:', slope, 'Intercept:', intercept, 'R-Value', r_value, 'P-Value', p_value, 'STD err', std_err
print 'Binned Graph ### R^2:', r_value**2

### Plot Data and fits
plt.figure(8, figsize=(9, 6.75))
fig = pylab.gcf()
plt.scatter(mean_bin_ranges, mean_bin_values, s=50, edgecolors='', color='black')
plt.plot([], [], ' ', label="Line Fit: y=m*x+b")
plt.plot(xdata, func_line(xdata, slope, intercept), 'r-', label='fit: m='+str(np.round(slope, 3))+' b='+str(np.round(intercept, 3)))
plt.plot([], [], ' ', label='R^2:'+str(np.round(r_value**2, 5)))
fig.canvas.set_window_title('')
plt.title('Total Intron Length vs Predicted: Binned Every '+str(delta), fontsize=size_font)
plt.legend(loc=1, fontsize=10)
plt.grid()
plt.ylabel('Mean Total Intron Length (nucleotides)', fontsize=size_font)
plt.xlabel('Mean Prediction', fontsize=size_font)
plt.show(block=True)
