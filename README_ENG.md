## FELIS for C-CAT database <img src="source/FELIS.png" width=50>
Functions Especially for LIquid and Solid tumor clinical sequencing, for Japanese only [(Japanese version of this README file)](https://github.com/MANO-B/FELIS/blob/main/README.md).  
Copyright (c) 2024 Masachika Ikegami, Released under the [MIT license](https://opensource.org/license/mit).  


### Analysis software for C-CAT utilization data
[The Cancer Genome Information Center (C-CAT)](https://www.ncc.go.jp/jp/c_cat/use/index.html) at the National Cancer Center has a database of cancer gene panel tests (Comprehensive Genomic Profiling, CGP test) conducted by the insurance. Cancer Genome Information Management Center (C-CAT)]() collects the results of cancer gene panel tests (Comprehensive Genomic Profiling, CGP test) and clinical information. There is a system for secondary use of this information for the purpose of academic research and drug development. Currently, it can only be used for research that has undergone ethical review by the institution to which it belongs and by C-CAT, and organizations other than hospitals and academia are required to pay an annual fee of 7.8 million yen to use it, which is a high threshold, but a similar overseas database, [AACR project GENIE](https://www.aacr.org/professionals/research/aacr-project-genie/), a similar overseas database, in terms of detailed drug information and clinical information, and is expected to enable analysis from a new angle in the study of rare cancers and rare fractions. The C-CAT data will be used in the study of rare cancers and rare fractions.  
  
The GUI-based software lowers the barrier to analysis, allowing clinicians to conduct clinical analysis based on clinical questions in their daily practice. Felis is the scientific name of the cat, and C-CAT-related naming seems to be tied to the cat's name.

Please understand that only those who can obtain data from C-CAT can use this software.  

### The analysis method is based on the following papers
> 1) Tamura T et al, Selection bias due to delayed comprehensive genomic profiling in Japan, Cancer Sci, 114(3):1015-1025, 2023.  
      See also [this website](https://github.com/MANO-B/CCAT) for more information on left-lateral truncation bias.
> 2) Mochizuki T et al, Factors predictive of second-line chemotherapy in soft tissue sarcoma: An analysis of the National Genomic Profiling Database. Cancer Sci, 115(2):575-588, 2024.  

### System Requirements
#### Hardware Requirements
For analysis of a few thousand cases, there is no problem, but for analysis of tens of thousands of cases, 32 GB or more of memory is required.    
Survival analysis is performed using Monte Carlo simulations with Stan, so a CPU with at least 4 cores and as fast as possible is recommended.  
RAM: 4+ GB  
CPU: 4+ cores  
  
Survival analysis of 3000 cases and 30 genes with 64 GB RAM and M1MAX MacStudio will take approximately 1 hour.  

#### Software Requirements
##### R language
Please refer to the [website](https://syunsuke.github.io/r_install_guide_for_beginners/03_installation_of_R.html) to install R as appropriate.  
Although no specific version is specified, this software was created using v4.3.2.  
Below is a link to [Start R from the command line and work with it](http://kouritsu.biz/installing-r-on-mac/).
##### Rstan
Please refer to this [RStan Getting Started (Japanese)](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started-(Japanese))  
````
# If you have already installed rstan, please execute the following line
# remove.packages(c(“StanHeaders”, “rstan”))

install.packages(“rstan”, repos = c(“https://mc-stan.org/r-packages/”, getOption(“repos”)))
````
##### Shiny
[Shiny](https://shiny.posit.co) was used to make it a web app.
```
install.packages(“shiny”))
```
##### Package dependencies
Install the dependent packages from the ``R`` terminal.  
The first time you do this, it may take quite a while.  
```
install.packages(c('ggplot2', 'umap', 'tidyr', 'dbscan', 'shinyWidgets', 'readr', 'dplyr', 'stringr', 'RColorBrewer', 'gt', 'gtsummary', ' 'flextable', 'Rediscover', 'survival', 'gridExtra', 'survminer', 'tranSurv', 'DT', 'ggsci', 'scales', 'patchwork', 'sjPlot', 'sjlabelled', ' forcats', 'markdown', 'PropCIs', 'shinythemes', 'BiocManager'), dependencies = TRUE)
BiocManager::install(c(“ComplexHeatmap”), update=FALSE)
```



### Starting FELIS
- Obtaining Analysis Files
First, download the case information you wish to analyze from the C-CAT Utilization Search Portal.
Select the Japanese version, not the English version, and then select the case, as shown in the following image.  
Report CSV (all data output)  
Download the following two files: Report CSV (all data output) and Case CSV (all data output).  
The ZIP file can be extracted and converted back to a CSV file for use.  
　　

- Download FELIS  
Download the ZIP file of the version of FELIS you wish to use, and download and unzip it to an appropriate folder.  
In this example, the folder is “/User/C-C-CAT/Desktop/felis-cs”.  

- Starting FELIS
Start the web application with the following command.  
```
$ R

R version 4.3.2 (2023-10-31) -- “Eye Holes”
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: aarch64-apple-darwin20 (64-bit)
.
.
.
. 'help.start()' will give you HTML browser help. 
Type 'q()' to exit R.

> library(shiny)
> runApp('/User/C-C-CAT/Desktop/felis-cs')
```
  
### Reading analysis files.
**Input C-CAT files** tab.  
Select and load the downloaded case CSV and report CSV from the Browse... button in the upper left corner of the screen. button in the upper left corner of the screen.  
Multiple files can be imported by selecting them.  
In addition, you can optionally enter a correspondence table to change the drug or tissue type.  

### Specify the analysis target.  
Open the **Input C-CAT raw file** tab.  
Click the Start File Read and Analysis Settings button to display the setting items.  
Many items can be set.  
  

#### Filter by histology  
- Filter by histology  
　　Narrow down the tissue types to be analyzed.  
- Tissue type to be analyzed together (If none, do not select)  
　　Select a group of tissue types to be analyzed as a single tissue type.  
- Name of the tissue type to be analyzed together  
　　Select a name that is representative of the tissue types to be analyzed together.  
- Minimum patients for each histology  
　　Rare histologies can be renamed to the site of occurrence for analysis.  
　　Set the minimum number of patients for each histology to be analyzed.  
  
#### Filters for Clinical Matters  
- Filter by sex  
　　Filter by sex for analysis.  
- Filter by panel  
　　Filter by cancer gene panel test for analysis.  
- Age for analysis  
　　Filter by age for analysis.  
- Threshold age for oncoprint  
　　Set the threshold of Young/Old classification in Oncoprint.  
- Filter by performance status  
　　Filter by PS for analysis  
- Filter by smoking status  
　　Filter by smoking status  
  
#### Filter by gene  
- Genes to be focused on (if not, do not select)  
　　Select the genes to be prioritized for Oncoprint, survival analysis, etc.  
- Gene sets of interest 1 (If there is no gene set, it is not selected.)  
　　Select the gene sets of interest, if any.  
- Select cases based on mutations?  
　　Only cases with or without mutations can be selected for analysis.  
  
#### Filter on type of mutation  
- Genes of interest (if none, unselected)  
　　Select a gene that you would like to detail the site and pattern of mutation in particular.
　　Example: EGFR TKD mutation  
- Mutation type  
　　Select mutations to be analyzed together as a single mutation pattern.  
- Mutation type name  
　　Name the mutation pattern.  
- Treatment of pathological significance of the gene of interest  
　　The pathological significance of only this gene to be analyzed can be changed.  
- Treat the specified mutation independently or  
　　Only the specified mutation is treated as a single gene.  
　　Example: Rename EGFR TKD mutation to EGFR_TKD gene  
  
#### Other settings  
- Gene number for oncoprint  
　　This setting is used to narrow down the genes to be included in the oncoprint and survival analysis.  
　　This will especially affect the time required for survival analysis.  
- Display of Oncoprint  
　　Sets the sort order in Oncoprint.
- Variants for analysis  
　　Select whether to analyze only oncogenic mutations or all mutations regardless of their pathological significance.
- Handling of Fusion genes
　　If there are many partner genes, the number of each gene will be reduced.  
　　Select whether to analyze all the fusion genes together, such as NTRK fusion and ALK fusion.  
- Distance value for DBSCAN clustering  
　　Set the threshold value of the distance to discriminate in the clustering analysis.  
- Treatment lines to be analyzed  
　　Specify the line of drugs to be analyzed.  
　　If only 1st-line is specified, comparison with the previous treatment will not be performed.  
   
### Run Analysis  
Open the Analysis tab.  
A number of analyses can be performed.  
The results are displayed on the tab corresponding to each button.  
   
  
#### Display case summary  
Displays a summary of the selected cases in the **Case Summary** tab.  
- Classify by mutation pattern and display in the **Patients summary, by mutation pattern** tab.  
- Displays a summary of the selected cases categorized by histology in the **Patients summary, by histology** tab.  
  
#### View Oncoprint  
- Displays the genetic variants of the selected case in the **Oncoprint** tab.  
Displays a table of cases in the **Clinical and mutation information per patient** tab. You can download the file from the button in the upper left corner.  

  
#### Display Mutual Exclusion/Co-mutation
- The results of the mutual exclusivity analysis of gene mutation senses using the [Rediscover package](https://academic.oup.com/bioinformatics/article/38/3/844/6401995) are displayed in the **Mutual Exclusion/Co-mutation** tab.  
    Blue means that the relationship is mutually exclusive and red means that the relationship is co-mutation.  
    An asterisk is displayed if P<0.001.  
  
#### Displays the mutation rate of each gene per histological type  
- For genes with high mutation frequency, the frequency of gene mutation per tissue type is displayed in the **Mutation per tissue type** tab.  
  
#### Clustering based on gene mutation  
Clustering based on mutated genes is performed using [UMAP](https://arxiv.org/abs/1802.03426) and [DBSCAN](https://cdn.aaai.org/KDD/1996/KDD96-037.pdf).
Results are displayed under the **Clustering** tab.  
- Basic information for each histology is displayed under the **Basic data** tab.  
    - Driver: Percentage of cases with at least one carcinogenic mutation detected.  
    - option_and_treat: Frequency (%) of cases that received treatment recommended by the expert panel.  
    - time_before_CGP: median survival (days) from the start of palliative chemotherapy to CGP testing.  
    - time_after_CGP: median survival (days) from CGP testing to death.  
- The histological types and genetic variants clustered in each cluster are displayed in the **Clustered histological types and variants** tab.  
    The histotypes clustered at P<0.05 are displayed in order of their odds ratios compared to the other clusters, up to three.  
    Genetic variants with P<0.05 are displayed in descending order of odds ratio compared to other clusters, up to 3.  
- Age groups in each cluster are displayed in the **Cluster and Age Relationship Diagram** tab.  
- Organizational type for each cluster is displayed in the **Cluster and Organizational Type Relationship Diagram** tab.  
- The entropy of the clusters is displayed in the **Organizational Types Clustered in Fewer Clusters** tab to indicate whether the organizational types are clustered in a few clusters or in many clusters.  
    The entropy is calculated by Shannon entropy. The lower the value, the more clustered it tends to be.  
- A table of the relationship between clusters and organizational types is displayed in the **Clusters and Organizational Types Table** tab.  
    You can download it from the button on the upper left.  
- A table of the relationship between clusters and genetic variation is displayed in the **Table of Clusters and Genetic Variation** tab.  
    It can be downloaded by clicking the button in the upper left corner.  
- The frequency of mutations detected in each histology corresponding to drugs with Evidence level is displayed in the **Frequency of patients with targeted therapy per histology** tab.  
    Consistent with high Evidence level: Patients with both Evidence level A and B drugs are designated as A.  
- The relationship between mean survival from start of chemotherapy to CGP testing and from CGP testing to death for each histology is displayed on the **Pre-CGP and Post-CGP Duration Relationship** tab.  
- The relationship between the number of cases of each histology and information about treatment is displayed on the **Number of Patients per Histology and Treatment Attainment Rate** tab.  
    The number of patients and the percentage with Evidence level A, B or higher, and C or higher drugs are shown as a scatterplot.  
    Scatterplot of number of patients and percentage with recommended treatment, percentage with recommended treatment, and percentage of patients with recommended treatment who received the recommended treatment.  
- The relationship between time to CGP testing and information about treatment for each histology type is displayed in the **Time to CGP and Percent Treatment Reached** tab.  
    Scatterplot of time to CGP and percentage of patients with Evidence level A, B or higher, and C or higher drugs.  
    Scatter plots of time to CGP testing and percentage with recommended treatment, percentage with recommended treatment, and percentage of patients with recommended treatment who received the recommended treatment.  
- The relationship between survival after CGP testing and information about treatment for each histology type is displayed in the **Survival After CGP and Treatment Attainment** tab.  
    Scatterplot of survival after CGP and percentage of patients with Evidence level A, B or higher, and C or higher drugs.  
    Scatter plots of survival after CGP testing and percentage of patients with recommended treatment, percentage of patients who received recommended treatment, and percentage of patients with recommended treatment who received recommended treatment.    
  
#### Survival analysis after CGP testing  
Survival analysis after CGP testing focusing on gene mutation, treatment details, PS, etc.
A 95% confidence interval will be calculated using log-log transformation.  
It shows which patients have a poor prognosis and early CGP testing is recommended.  
The results are displayed under the tab **Survival after CGP**.  
- Survival analysis grouped by treatment recommendations is displayed under the **Survival and Treatment after CGP** tab.  
    UnrecomTreat(+): Patients who received treatment other than recommended treatment  
    RecomTreat(+): Patients who received the recommended treatment  
    Treat(-): Patients who did not receive any treatment after CGP testing  
- Survival analysis grouped by histology, PS, and presence/absence of genetic mutations is displayed in the **Survival and PS after CGP** tab.  
    For gene mutation analysis, if there is a gene of interest, the patients are grouped according to whether or not they have a mutation in one of the genes.  
    If there is no gene of interest, the mutations are grouped according to whether or not the gene with the highest mutation frequency is mutated.  
- The median survival (days) of the patients grouped by gene mutation is displayed in the **Survival and gene mutation after CGP, forest plot** tab.  
    If the mutation frequency is small, the 95% confidence interval is not displayed.  
- Kaplan-Meier survival curves grouped by gene mutation are displayed in the **Survival after CGP and Gene Mutation, KM-curve** tab.  
- A forest plot of the hazard ratio for death after CGP testing is displayed in the **Hazard ratio for survival after CGP, forest plot** tab.  
    Factors that match more than 95% are considered multicollinear and are excluded.
    Factors with more than 95% concordance are excluded as multicollinearity.  
    If a factor has no mortality events, the result is not displayed.  
- The results of the logistic regression analysis of the factors that led patients to the **recommended treatment** are displayed as a forest plot in the tab **Factors leading to the recommended treatment, forest plot**.  
    If there are factors that have not reached the recommended treatment at all, no values will be displayed.  
- Displays a summary table on the **Factors leading to recommended treatment, table** tab with a summary of the factors leading to **patients receiving recommended treatment**.  
  
#### Survival analysis of chemotherapy induction (takes time)  
Perform a survival analysis after palliative chemotherapy induction with left-sided cutting bias.
The analysis takes time in the order of tens of minutes due to the simulation using Stan.
Results are displayed under the **Survival after CTx** tab.  
- Survival analysis for the cases corrected for left-sided truncation bias, corrected for number at risk, and corrected for simulation is displayed under the tab **Overall Survival after CTx corrected for left-sided truncation bias**.  
- Survival analysis with grouping by gene mutation is displayed in the **Overall Survival After CTx Corrected for Left-Sided Cutting Bias** tab.  
    The gene mutation analysis is grouped by whether or not any of the genes of interest have a mutation.  
    If there is no gene of interest, the patients are grouped according to whether or not the gene with the highest mutation frequency has a mutation.  
- The difference (days) of median survival is calculated for each gene mutation and the result is displayed in the **Gene mutation and overall survival, forest plot** tab.  
    If there are few mortality events, the results will not be displayed.  
- Survival curves, grouped by gene mutation, are displayed in the **Gene mutation and overall survival, KM-curve** tab.  
  
#### List of drugs used in Palliative CTx (1st-4th line)  
- Regimens used in the 1st-4th line of palliative chemotherapy are extracted and displayed in the **Drug response and Drug table** tab.
The input may appear to be inaccurate, so it should be used only to see trends.  
**Select regimens of interest** panel will appear. Select regimens of interest for subsequent analysis.  
  
#### Response analysis of the drug selected by the above button  
Evaluate the relationship between the duration of drug response and genetic mutations and histology with a focus on Treatment on time (ToT).  
Results are displayed under the **Drug Response** tab.  
- Information on all drugs is summarized by line of therapy and displayed under the **Drug Usage, by Line of Therapy** tab.  
- All drug information is summarized by therapeutic effect and displayed under the **Drug Usage, By Therapeutic Effect** tab.  
- The information of drugs in the specified line is summarized by mutation pattern and displayed in the **Drug Usage in the Specified Line, by Mutation Pattern** tab.  
- The information of the drug for the specified line is summarized by tissue type and displayed in the **Drug Usage for the Specified Line, by Tissue Type** tab.  
- The drug information of the specified line is summarized by gene mutation of interest and displayed in the **Drug Usage of the Specified Line, by Gene Mutation** tab.  
- Information on the specified line and the specified drug with ToT information is summarized by mutation pattern and displayed in the **Specified line with ToT information, usage of the specified drug, by mutation pattern** tab.  
- The information on the designated line and the designated drug with ToT information is summarized by tissue type and displayed in the **Display of the designated line and the designated drug with ToT information, by tissue type** tab.  
- Information on designated lines and designated drugs with ToT information is summarized by gene mutation of interest and displayed in the **Display of designated lines and designated drugs with ToT information, by gene mutation** tab.  
- Information on the designated line and the designated drug with RESICT information is summarized by mutation pattern and displayed in the **Display of the designated line and the designated drug with RESICT information, by mutation pattern** tab.  
- Information on designated lines and designated drugs with RESICT information is summarized by tissue type and displayed in the **Display of designated lines and designated drugs with RESICT information, by tissue type** tab.  
- Information on designated lines and designated drugs with RESICT information is summarized by gene mutation of interest and displayed in the **Display of designated lines and designated drugs with RESICT information, by gene mutation** tab.  
- Displays waterfall plots and scatter plots of the relationship between the ToT of the drug of interest and the ToT of its previous treatment in the **Time on treatment, scatter plots** tab for the specified treatment and its previous treatment.  
- Scatter plots of the relationship between the ToT of the drug of interest and the ToT of its previous treatment are displayed in the **Time on treatment, scatter plots** tab for the specified treatment and its previous treatment.
    Censored cases are excluded.  
- A Kaplan-Meier survival curve of the ToT of the drug of interest compared to the ToT of its previous treatment and the ToT of other drugs, and the relationship between gene mutation and ToT is displayed in the **Time on treatment for the specified treatment and its previous treatment, KM-curve** tab.  
- Displays Kaplan-Meier survival curves for the relationship between ToT and genes and mutation patterns of interest in the **Time on treatment, KM-curve** tab for each gene and mutation of interest.  
- Displays Kaplan-Meier survival curves for the ToT for the drug and tissue type of interest in the **Time on treatment by tissue type** tab.  
- Displays the median OS forest plot of the ToT for the drug of interest and the gene mutation in the **Time on treatment, forest plot** tab for each gene mutation.  
- Displays the Kaplan-Meier survival curve of the ToT for the drug and gene mutation of interest in the **Time on treatment by gene mutation, KM-curve** tab.  
- Displays a table of Hazard ratios for factors leading to treatment discontinuation in the **Hazard ratio for time on treatment** tab.  
- Displays a forest plot of Odds ratio for factors leading to Objective response in the **Odds ratio for ORR, forest plot** tab.  
- Displays a table of Odds ratio of factors leading to Objective response in the **Odds ratio for ORR, table** tab.  
- Displays a forest plot of Odds ratio leading to disease control in the **Odds ratio for DCR, forest plot** tab.  
- Displays a table of Odds ratio of factors leading to disease control in the **Odds ratio for DCR, table** tab.  
- Displays a table of response by mutation pattern in the **Mutation pattern and RECIST** tab.  
    Calculate 95% confidence intervals using the Clopper-Pearson method.  
- A table of responses by histology is displayed in the **HISTOTYPE AND RECIST** tab.  
- A table of response by gene mutation is displayed on the **Gene Mutation and RECIST** tab.  
  
#### Displays instructions on how to use the software, etc. on the **DISCUSSION** tab.
　　
### Future Plans
- Add mutual exclusivity analysis between pathways  
- Add analysis to evaluate the association of RECIST response with clustering, histology, genetic mutations, etc.  
- Added an analysis to determine which is more predictive of drug response, genetic mutations from panel testing or results of tests performed prior to panel testing, such as HER2 immunostaining, MSI, etc.  
- Added analysis of the association between variant frequency and drug response in liquid sequencing.
- Evaluated differences in mutation clustering patterns between tissue types, including VUS.  

### Recommended FELIS versions for each version of the C-CAT database  
Since column names may be added or changed in each version of the C-CAT data, the appropriate version of FELIS is required.  
C-CAT database version 20240621: FELIS version 1.2.2  
  
### Version history
1.2.2: Added analysis by gene mutation pattern, corrected miscounting of total cases in survival analysis - 20240822  
1.2.1: Added support for analysis by gene mutation pattern, such as Exon19 mutation and TKD mutation - 20240821  
1.2.0: Improved user interface, corrected error in tissue types with only one gender - 20240821  
1.1.2: Added analysis of drug response - 20240820  
1.1.1: Fix error with only one diagnosis - 20240818  
1.1.0: Improved User Interface, addresses problem of poorly annotated oncogenic mutations, especially EGFR - 20240817  
1.0.0: Support C-CAT database version 20240621 - 20240815  
  
