# Intron-Analysis
Walk through for the analysis of the kinetics of intron length.
## Reproduction of Results – Walkthrough
All python programs were programmed with the default libraries acquired from a base Anaconda installation.
All R scripts were programmed with the installation of R and acquired libraries as outlines in the Mus musculus walkthrough.
NOTE: Both the installation programs required for R and R Studio were too large for uploading to github, as such the exact installation programs will be avaliable upon request (tisdalegarrett@gmail.com)
### Mus musculus:
1.	Acquired supplementary file 3 from Global quantification of mammalian gene expression control by Bjorn Schwanhausser et al. (https://www.nature.com/articles/nature10098). File name downloaded as nature10098-s5.xls and saved into the Mus musculus folder.
2.	Copied data column wise selecting columns “MGI ID”, “mRNA copy number average [molecules/cell]”, “mRNA half-life average [h]”, and “transcription rate (vsr) average [molecules/(cell*h)]” into a new csv file for further parsing named “Schwanhausser et al – unparsed.csv” saved into the Mus musculus folder.
3.	To aid in the recreation of this analysis for third parties R and R studio were uninstalled.
4.	RStudio 1.2.1335 – Windows 7+ (64-bit) installer (release data 2019-04-08) was download and installed with general recommended setup.
    - Ran R studio and verified no current, unknown, installation of R was found.
    - Installer saved to ReAnalysis folder as RStudio-1.2.1335.exe
5. Downloaded R installer from https://cran.r-project.org/, selecting download R for Windows, install R for the first time, and finally Download R 3.6.1 for Windows.
    - Installer saved to ReAnalysis folder as R-3.6.1-win.exe
6.	Installed R on second hard drive, noting the old installation data was still present (but not registered by default by windows as previously established in step 4.a). Deleted old folder and installed the new version in its place.
    - Used all default settings including save version number in registry and associate R with .RData files
7.	 Printed “Hello World” to confirm that default R was installed and working in the RStudio environment.
8.	Created new R notebook and installed the necessary files as demonstrated in analysis-package-instillations.Rmd located in the Mus musculus folder (always downloaded or updated all (“a”) requested fields)
    - Installed BiocManager for utilizing Bioconductor packages
    - Installed TxDb.Mmusculus.UCSC.mm10.knownGene for accessing the mouse genome as described https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Mmusculus.UCSC.mm10.html
    -   Installed BiocParallel due to the default instillations within the other packages causing unnecessary struggle and program failure during the original analysis (installing directly seemed to fix these issues)
    -   Installed org.Mm.eg.db for accessing Entrez ID numbers for the conversion between MGI IDs to Entrez IDs as detailed in https://bioconductor.org/packages/release/data/annotaorg.Mm.eg.dbtion/html/org.Mm.eg.db.html
9.	All detailed steps required to retrieve intron and exon data are found within the code in Intron-and-Exon-Length-Data-Parsing.Rmd. Final output of which is Data_for_Graph_-_Mouse_-_unparsed.csv.
    -   Manually added headers to Data_for_Graph_-_Mouse_-_unparsed.csv
10.	All detailed steps required to produce the graphs in figure 2 of the paper are found within the code Graphing_Data_-_Figure_2.py
     -   Note: The one outlier data point as mentioned in the text for Figure 2 C and D was not excluded for transparency
11.	All detailed steps required to produce the graphs in figure 3 of the paper are found within the code Graphing_Data_-_Figure_3.py
### Human (Normal and Cancer):
Starting under the assumption that steps 1-8 are already carried out in the Mus musculus walk though.
1.	Created new R notebook the Human genome annotation R package as demonstrated in analysis-package-installations.Rmd located in the Human (Normal and Cancer) folder (updating all (“a”) requested fields)
    - Installed TxDb.Hsapiens.UCSC.hg19.knownGene as described in https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html
2.	Obtained the expression data for the Human control and Human cancer gene expression data from the “Supplementary file” section found on the Gene Expression Omnibus (GEO) at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133624
    - File was downloaded and saved as listed (GSE133624_reads-count-all-sample.txt.gz) into the Human (Normal and Cancer) folder
3.	Created new Python program named Extracting-Human-Control-and-Cancer-Gene-Expression-Data.py to extract data out of the compressed *.gz file. All detailed steps are found within the program. Program outputs Human-Normal-and-Cancer-Mean-Gene-Expression-Data.csv
    - Manually added headers to Human-Normal-and-Cancer-Mean-Gene-Expression-Data.csv for clarity
4.	Created new R notebook named ID_Transforms.Rmd to transform Ensembl IDs to Entrez IDs and match all shared IDs between the Human annotation package used and the IDs provided for the genes provided. All detailed steps are found within the program.
    - Data file was created named Human-Normal-and-Cancer-Mean-Gene-Expression-Data_With-Common-Ids.csv to then be used to obtain the intron (and unused exon) data
    - Manually added headers to Human-Normal-and-Cancer-Mean-Gene-Expression-Data_With-Common-Ids.csv
5.	Created new R notebook named Intron-and Exon-Length-Data-Parsing.Rmd to retrieve intron (and unnecessarily exon) data. All detailed steps are found within the program.
    - Manually added headers to Data_for_Graph_-_Human_(Control-and-Cancer)_-_unparsed.csv
6.	Created new Python program named Graphing_Data_-_Figure_5 to further parse and graph the results. All detailed steps are found within the program.
### Human (Castillo-Davis et al.)
1.	Obtained raw data files for Castillo-Davis et al. figure 1 A and B  (article found at: https://www.nature.com/articles/ng940z). Obtained from Cristian Castillo-Davis and saved as Raw-Data-Transcription-Selection-for-Garrett.xlsx in Human (Castillo-Davis et al) folder.
2.	Copied data from all Human data as provided in “Human – all data” section on the second tab of the provided *.xlsx spreadsheet and saved the entire section into a new file named “Human_-Data_-_Castillo-Davis.csv”
3.	Created new python program for graphing the results data in Human (Castillo-Davis et al) folder named Graphing_Data_-_Figure_4_B.py. All detailed steps are found within the program.
### C Elegan/Nematode (Castillo-Davis et al.)
Same analysis style as in Human (Castillo-Davis et al.)
1.	Obtained raw data files for Castillo-Davis et al. figure 1 A and B  (article found at: https://www.nature.com/articles/ng940z). Obtained from Cristian Castillo-Davis and saved as Raw-Data-Transcription-Selection-for-Garrett.xlsx in Human (Castillo-Davis et al) folder.
2.	Copied data from all Human data from the provided “Nematode – ALL DATA” section on the first tab of the provided *.xlsx spreadsheet and saved the entire section into a new file named “Nematode_-Data_-_Castillo-Davis.csv”
3.	Created new python program for graphing the results data in C Elegan (Nematode) (Castillo-Davise et al) folder named Graphing_Data_-_Figure_4_A.py. All detailed steps are found within the program.
### Chlamydomonas reinhardtii
1.	Obtained raw genome annotation data file and expression data file from Joint Genome Institute (Chlamydomonas reinhardtii v5.5).
    - Required a login to access 
    - At https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Creinhardtii selected “Bulk data”
    - Selected “OK, proceed to data” at the Data Usage Policy page
    - Under Creinhardtii, v5.5, annotation folder directory downloaded Creinhardtii_281_v5.5.gene_exons.gff3.gz file and saved to Chlamydomonas reinhardtii folder
    - Under Creinhardtii, v5.5, expression folder directory downloaded GSE25124_genes.fpkm_tracking file and saved to Chlamydomonas reinhardtii folder
2.	Created new Python program to extract the FPKM data from the expression data set named Extracting_FPKM_Expression_Data.py saved in Chlamydomonas reinhardtii folder. All detailed steps are found within the program. Saved extracted data in the Chlamydomonas reinhardtii folder as Chlamydomonas_reinhardtii_Extracted_Expression_Data.csv.
3.	Created new Python program to interpolate and save intron length data named Extracting_Intron_Length_Data.py saved in the Chlamydomonas reinhardtii folder. All detaied steps are found within the program. Saved extracted data in Chlamydomonas rienhardtii folder as Chlamydomonas_reinhardtii_Extracted_Intron_Data.csv.
    - Manually added headers to csv for clarity
4.	Created new Python program named Graphing_Data_-_Figure_4_C.py in the Chlamydomonas reinhardtii folder to merge the expression and intron data sets to graph the collected data. All detailed steps are found within the program.
### Dictyostelium discoideum
1.	Obtained genome annotations of Disctyostelium discoideum from NCBI (https://www.ncbi.nlm.nih.gov/assembly/GCF_000004695.1/#/def)
    - Using RefSeq assembly accession: GCF_000004695.1
    - In the Access the data list on the right side of the page clicked Download the RefSeq assembly
    - Downloaded GCF_000004695.1_dicty_2.7_genomic.gff.gz and saved file within Dictyostelium discoideum folder
2.	Obtained gene expression data from Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71036)
    - From the samples drop down list manually selected and transferred data from the to be listed sample links by clicking the “View full table…” button at the bottom of the linked pages and selecting all of the data on the page and transferring it to a csv file named Dictyostelium_discoideum_Extracted_Expression_Data.csv
       - wt Amoeboid rep1: GSM1825931
       - wt Early Culminant rep1: GSM1825932
       - wt Fruiting Body rep1: GSM1825933
       - wt Loose Aggregate rep1: GSM1825934
       - wt Mound rep1: GSM1825935
       - wt Slug rep1: GSM1825936
       - wt Streaming rep1: GSM1825937
       - wt T0 rep1: GSM1825938
    - Manually raised each log values expression to the power of 10 (as no log scale operations were provided) and averaged across all expression states for mean expression values.
3.	Created new python program to interpolate and save intron length data named Extracting_Intron_Length_Data.py saved in the Dictyostelium discoideum folder. All detailed steps are found within the program. Saved extracted data in Dictyostelium discoideum folder as Dictyostelium_discoideum_Extracted_Intron_Data.csv.
    - Manually added header replacing a full zero row appearing in the Dictyostelium_discoideum_Extracted_Intron_Data.csv data file appearing due to and extra nonessential data point in the raw annotations file
4.	Created new Python program named Graphing_Data_-_Figure_4_D.py in the Dictyostelium discoideum folder to merge the expression and intron data sets to graph the collected data. All detailed steps are found within the program.

