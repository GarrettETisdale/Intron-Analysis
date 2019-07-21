### Required imports
import gzip
import csv
from itertools import izip

### Import exon annonataion information attained from Gene Expression Omnibus
f = gzip.open(r'Creinhardtii_281_v5.5.gene_exons.gff3.gz','rb')
file_content = f.read()

### Dictionary to hold and organize data for future calling
d_genes_exons = {}

### Split raw data file for iteration
split_file_content = file_content.split('\n')

### Iterate through all annotations
for i in split_file_content[2:-1]:
    split_1 = i.split('\t')           # split line by tab space ("\t") data sections
    split_2 = split_1[-1].split(';')  # split last tab section by ";" to isolate gene IDs
    split_3 = split_2[0].split('=')   # split the first ID by "=" to isolate specific identifier
    split_4 = split_3[1].split('.')   # split name for proper Gene ID

    gene_id = split_4[0]+'.'+split_4[1]  # recompile gene ID

    ### Append the start and end of the base pair locators if gene is already in dictionary
    if gene_id in d_genes_exons.keys():
        if split_1[2] == 'exon':
            start = min(float(split_1[3]), float(split_1[4]))
            end = max(float(split_1[3]), float(split_1[4]))

            d_genes_exons[gene_id].append([start, end])
        pass
    ### First ID encounter is not an exon, as such add gene ID before first encountered exon
    else:
        temp_max = max(float(split_1[3]), float(split_1[4]))
        temp_min = min(float(split_1[3]), float(split_1[4]))
        length = temp_max-temp_min                              # Total length of exon region

        d_genes_exons[gene_id] = []
        print length

        print split_1  # [0], gene_id
        # print split_1  # Print update to indicate where the program is at

gene_name_list = []      # List to contain the condensed dictionary gene names
exon_length_list = []    # List to contain the condensed dictionary exon lengths
intron_length_list = []  # List to contain the condensed dictionary intron lengths

### Iterate trough all dictionary items to condense and reformat data for saving
for i in d_genes_exons.keys():
    num_exons = len(d_genes_exons[i])  ### Number of exon regions for iteration reference

    intron_len_count_temp = 0  ### counter to keep track of intron length
    exon_len_count_temp = 0    ### counter to keep track of exon length

    ### Continue this iteration for all except the last exon region
    for x in range(0, num_exons-1):
        intron_max = max(d_genes_exons[i][x+1][0], d_genes_exons[i][x][1])  ### reference for end of intron (beggining of next exon)
        intron_min = min(d_genes_exons[i][x+1][0], d_genes_exons[i][x][1])  ### reference for start of intron (beggining of next exon)

        intron_len_count_temp += intron_max - intron_min                        ### Intron segment length
        exon_len_count_temp += d_genes_exons[i][x][1] - d_genes_exons[i][x][0]  ### Total intron length

        # For the second to last region perform operation detailing the last range in series
        if x == num_exons-2:
            exon_len_count_temp += d_genes_exons[i][x+1][1] - d_genes_exons[i][x+1][0]  ### Add last remaining exon length region

    ### Append gene data to save lists
    gene_name_list.append(i)
    intron_length_list.append(intron_len_count_temp)
    exon_length_list.append(exon_len_count_temp)

### Save data
with open('Chlamydomonas_reinhardtii_Extracted_Intron_Data.csv', 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(izip(gene_name_list, intron_length_list, exon_length_list))
