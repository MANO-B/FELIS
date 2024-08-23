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
  
