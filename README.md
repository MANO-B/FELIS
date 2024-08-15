# FELIS for C-CAT database <img src="FELIS.png" width=50>
Functions Especially for LIquid and Solid tumor clinical sequencing, for Japanese only.

## C-CAT利活用データの解析ソフトウェア
国立がん研究センターに設置されている[がんゲノム情報管理センター(C-CAT)](https://www.ncc.go.jp/jp/c_cat/use/index.html)には保険診療で行われたがん遺伝子パネル検査(Comprehensive Genomic Profiling, CGP検査)の結果と臨床情報が集約されています。この情報を学術研究や医薬品等の開発を目的とした二次利活用する仕組みがあります。現状では所属施設の倫理審査とC-CATでの倫理審査を経た研究でのみ使用可能であり、また病院やアカデミア以外の組織では年間780万円の利用料金が必要と敷居が高いですが、類似した海外のデータベースである[AACR project GENIE](https://www.aacr.org/professionals/research/aacr-project-genie/)と比較して薬剤の情報や臨床情報が詳しい点で優れており、希少がん・希少フラクションの研究においてこれまでになかった切り口での解析が可能になると考えられています。  
  
C-CATのデータを用いるに当たってはビッグデータかつリアルワールドデータの解析には特有の問題があり、また一定程度のデータ処理を行うプログラミングの知識が必要になります。GUIを用いたソフトウェアにより解析の敷居を下げることで、臨床医の日常診療におけるクリニカルクエスチョンに基づいた探索的研究を容易とし、C-CAT利活用データの活用を促進するために本ソフトウェアを作成しました。命名にネコの名前縛りがあるようです。

## 解析手法は以下の論文に基づきます
> 1) Tamura T et al., Selection bias due to delayed comprehensive genomic profiling in Japan, Cancer Sci, 114(3):1015-1025, 2023.  
      左側切断バイアスについては[こちらのwebsite](https://github.com/MANO-B/CCAT)も参照ください。
> 2) Mochizuki T et al., Factors predictive of second-line chemotherapy in soft tissue sarcoma: An analysis of the National Genomic Profiling Database, Cancer Sci, 115(2):575-588, 2024.  

## System Requirements
### Hardware Requirements
数千例の解析であれば問題ありませんが、数万例の解析を行う場合は32GB以上のメモリが必要です。    
生存期間解析はStanを用いたモンテカルロ法でのシミュレーションを行います。4コア以上でできるだけ高速なCPUの使用が望まれます。  
RAM: 4+ GB  
CPU: 4+ cores  
  
3000例、30遺伝子についての生存期間解析を64 GB RAM, M1MAX MacStudioで行った場合、およそ1時間を要します。  

### Software Requirements
#### R language
適宜[ウェブサイト](https://syunsuke.github.io/r_install_guide_for_beginners/03_installation_of_R.html)を参照しRを導入ください。  
特にバージョンの指定はありませんが、本ソフトウェアはv4.3.2を使用して作成しました。  
以下、[コマンドラインからRを起動して作業を行います。](http://kouritsu.biz/installing-r-on-mac/)  
#### Rstan
こちらの[RStan Getting Started (Japanese)](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started-(Japanese))を参照ください。  
```
# もしすでにrstanをインストールしているならば次の行を実行してください
# remove.packages(c("StanHeaders", "rstan"))

install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```
#### Shiny
Webアプリとするために[Shiny](https://shiny.posit.co)を使用しました。
```
install.packages("shiny")
```
#### Package dependencies
依存しているパッケージ群を`R`ターミナルからインストールください。  
初めて実行する場合は相当に時間がかかると思われます。  
```
install.packages(c('ggplot2', 'umap', 'tidyr', 'dbscan', 'shinyWidgets', 'readr', 'dplyr', 'stringr', 'RColorBrewer', 'gt', 'gtsummary', 'flextable', 'Rediscover', 'survival', 'gridExtra', 'survminer', 'tranSurv', 'DT', 'ggsci', 'scales', 'patchwork', 'sjPlot', 'sjlabelled', 'forcats', 'BiocManager'), dependencies = TRUE)
BiocManager::install(c("ComplexHeatmap"), update=FALSE)
```

## 使用の流れ
- 解析ファイルの入手
まずは解析したい症例の情報をC-CAT利活用検索ポータルからダウンロードします。
症例を選択した上で、以下の画像の通り  
・レポートCSV（全データ出力）  
・症例CSV（全データ出力）  
の2つのファイルをダウンロードします。  
<img src="report.png"  height=300>      <img src="case.png" height=300>

- FELISのダウンロード
使用するバージョンのFELISのZIPファイルをダウンロードし、適当なフォルダにダウンロード・解凍してください。
ここでは"/User/C-CAT/Desktop/felis-cs"とします。  

- FELISの起動
以下のコマンドでWebアプリが起動します。  
```
$ R

R version 4.3.2 (2023-10-31) -- "Eye Holes"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: aarch64-apple-darwin20 (64-bit)
.
.
.
'help.start()' で HTML ブラウザによるヘルプがみられます。 
'q()' と入力すれば R を終了します。

> library(shiny)
> runApp('/User/C-CAT/Desktop/felis-cs')
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
