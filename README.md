# FELIS for C-CAT database <img src="FELIS.png" width=50>
Functions Especially for LIquid and Solid tumor clinical sequencing, for Japanese only.

## C-CAT利活用データの解析ソフトウェア
国立がん研究センターに設置されている[がんゲノム情報管理センター(C-CAT)](https://www.ncc.go.jp/jp/c_cat/use/index.html)には保険診療で行われたがん遺伝子パネル検査(Comprehensive Genomic Profiling, CGP検査)の結果と臨床情報が集約されています。この情報を学術研究や医薬品等の開発を目的とした二次利活用する仕組みがあります。現状では所属施設の倫理審査とC-CATでの倫理審査を経た研究でのみ使用可能であり、また病院やアカデミア以外の組織では年間780万円の利用料金が必要と敷居が高いですが、類似した海外のデータベースである[AACR project GENIE](https://www.aacr.org/professionals/research/aacr-project-genie/)と比較して薬剤の情報や臨床情報が詳しい点で優れており、希少がん・希少フラクションの研究においてこれまでになかった切り口での解析が可能になると考えられています。  
  
C-CATのデータを用いるに当たってはビッグデータかつリアルワールドデータの解析には特有の問題があり、また一定程度のデータ処理を行うプログラミングの知識が必要になります。GUIを用いたソフトウェアにより解析の敷居を下げることで、臨床医の日常診療におけるクリニカルクエスチョンに基づいた探索的研究を容易とし、C-CAT利活用データの活用を促進するために本ソフトウェアを作成しました。

## 解析手法は以下の論文に基づきます
> 1) Tamura T et al., Selection bias due to delayed comprehensive genomic profiling in Japan, Cancer Sci, 114(3):1015-1025, 2023.  
      左側切断バイアスについては[こちらのwebsite](https://github.com/MANO-B/CCAT)も参照ください。
> 2) Mochizuki T et al., Factors predictive of second-line chemotherapy in soft tissue sarcoma: An analysis of the National Genomic Profiling Database, Cancer Sci, 115(2):575-588, 2024.  

## System Requirements

### Hardware Requirements

The scripts requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 4 GB of RAM (depending on the size of BAM file and the number of mutations). For optimal performance, we recommend a computer with the following specs:

RAM: 4+ GB  
CPU: 2+ cores, 2.6+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, M1 Macbook air) and internet of speed 40 Mbps.

### Software Requirements

### Samtools

Samtools is used for pre-processing to remove reads that are not related to mutations. Version 1.12 is what I am using, but I think older versions will work if they support multi-core processing.  

### R language

This script files runs on `R` for Windows, Mac, or Linux, which requires the R version 3.5.0 or later.
If you use version 3.4 or lower of R, you will have some difficulty installing the packages, but it is not impossible.  

### Package dependencies

Users should install the following packages prior to use the scripts, from an `R` terminal:

```
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
install.packages(c('stringr', 'dplyr', 'remotes'), dependencies = TRUE)
BiocManager::install(c("Rsamtools", "Biostrings", "GenomicAlignments", "GenomeInfoDb"), update=FALSE)

# install necessary genomes
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update=FALSE)
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", update=FALSE)
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", update=FALSE)
```

which will install in about 30 minutes on a recommended machine.

### Package Versions
The program does not use any of the functions specific to the following version of the packages, so there is no problem if you use the latest version of the package.  

```
> packageVersion("stringr")
[1] '1.4.0'
> packageVersion("dplyr")
[1] '1.0.2'
> packageVersion("Biostrings")
[1] '2.54.0'
> packageVersion("BSgenome.Hsapiens.UCSC.hg38")
[1] '1.4.1'
> packageVersion("GenomicAlignments")
[1] '1.22.1'
> packageVersion("Rsamtools")
[1] '2.0.3'
> packageVersion("remotes")
[1] '2.5.0'
> packageVersion("GenomeInfoDb")
[1] '1.22.1'
```

## Instructions for Use
See also https://rdrr.io/cran/MicroSEC/
- How to install
```
# Stable version (v2.1.5) from github (recommended)
remotes::install_github("MANO-B/MicroSEC", ref = 'v2.1.5')

# Developmental version from github
remotes::install_github("MANO-B/MicroSEC")

# Stable version (v2.1.3) from CRAN (not recommended)
install.packages('MicroSEC', dependencies = FALSE)

```
- How to use in command line
```
# download only once
wget https://raw.githubusercontent.com/MANO-B/MicroSEC/main/MicroSEC.R
# run the script
Rscript MicroSEC.R [working/output directory] [sample information tsv file] [progress bar Y/N]
```  
- Example
```
Rscript MicroSEC.R /mnt/HDD8TB/MicroSEC /mnt/HDD8TB/MicroSEC/source/Sample_list.txt Y  
```  
- How to use in R Console
```
```

# Version history
1.0.0: C-CAT database version 20240621に対応  
