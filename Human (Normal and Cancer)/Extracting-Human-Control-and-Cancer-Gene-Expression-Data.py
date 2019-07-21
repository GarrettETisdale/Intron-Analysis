### Necessary Imports
import csv
import gzip
import numpy as np
from itertools import izip

f = gzip.open(r'GSE133624_reads-count-all-sample.txt.gz','rb')
file_content = f.read()

split_file_content = file_content.split('\n')

### Print the first 5 lines for view of format (N#### = Control and T#### = Cancer)
for i in split_file_content[:5]:
    print i


gene_ids = []
control_expression_average = []
cancer_expression_average = []
### Itterating though file, skipping the header and last blank line, while saving data in the above variables
for i in split_file_content[1:-1]:#:
    split_1 = i.split('\t')
    gene_ids.append(split_1[0])

    control_data = split_1[1:30]
    cancer_data = split_1[30:]

    control_list = []
    cancer_list = []

    for x in range(0, len(control_data)):
        control_list.append(float(control_data[x]))
    for x in range(0, len(cancer_data)):
        cancer_list.append(float(cancer_data[x]))

    control_expression_average.append(np.mean(control_list))   ### Mean of Control data using replicates
    cancer_expression_average.append(np.mean(cancer_list))     ### Mean of Cancer data using replicates


### Write Gene IDs, mean control expression data, and mean cancer expression data to csv file format
with open('Human-Normal-and-Cancer-Mean-Gene-Expression-Data.csv', 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(izip(gene_ids, control_expression_average, cancer_expression_average))

