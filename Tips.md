## FELIS 解析例 <img src="source/FELIS.png" width=50>
必要性の高いであろう解析について方法をまとめてみました。  
  
  
### 薬剤奏効性解析
<details>
<summary>ある治療と、一つ前のコースの治療との効果比較をしたい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢、治療コースなどの絞り込みを行う<br>  
3. Analysis -> Drug response analysis -> List of drugs used in Palliative CTxボタンを押す<br>  
4. Results -> Drug response -> Tables -> Drug use, by line of treatmentから使用状況を確認する<br>  
5. Analysis -> Drug response analysis -> Choose drugs for treatment effect analysisで治療を選択する<br>
6. Analysis -> Drug response analysis -> Analyze with the setting selected aboveボタンで解析を行う<br>
7. Results -> Drug response -> Time on treatment -> Time on treatment and pre-treatment for the specified treatment, scatter plotで、同一患者群での前治療と指定治療のtime on treatmentの比較を行う<br>
8. Results -> Drug response -> Time on treatment -> Time on treatment and pre-treatment for the specified treatment, KM-curveで、同一患者群での前治療と指定治療のtime on treatmentの比較を行う<br>
<br>

</details>
  
<details>
<summary>組織型ごとの治療効果をみたい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢、治療コースなどの絞り込みを行う<br>  
3. Analysis -> Drug response analysis -> List of drugs used in Palliative CTxボタンを押す<br>  
4. Results -> Drug response -> Tables -> Drug use, by line of treatmentから使用状況を確認する<br>  
5. Analysis -> Drug response analysis -> Choose drugs for treatment effect analysisで治療を選択する<br>
6. Analysis -> Drug response analysis -> Analyze with the setting selected aboveボタンで解析を行う<br>
7. Results -> Drug response -> Time on treatment -> Time on treatment and pre-treatment for the specified treatment, KM-curveで、遺伝子変異の有無で群分けした指定治療のTime on treatmentをKaplan-Meier法で評価する<br>
<br>
全ての薬剤での治療期間と指定薬剤での治療期間に差がある場合、その遺伝子変異が指定薬剤のbiomarkerである可能性が示唆されます。
</details>
  
<details>
<summary>ある薬剤の効果における遺伝子変異の有無の意義をみたい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢、治療コースなどの絞り込み、探索したい遺伝子の指定を行う<br>  
3. Analysis -> Drug response analysis -> List of drugs used in Palliative CTxボタンを押す<br>  
4. Results -> Drug response -> Tables -> Drug use, by line of treatmentから使用状況を確認する<br>  
5. Analysis -> Drug response analysis -> Choose drugs for treatment effect analysisで治療を選択する<br>
6. Analysis -> Drug response analysis -> Analyze with the setting selected aboveボタンで解析を行う<br>
7. Results -> Drug response -> Time on treatment -> Time on treatment by tissue type, KM-curveで、全ての治療あるいは指定治療のTime on treatmentをKaplan-Meier法で評価する<br>
<br>
全ての薬剤での治療期間と指定薬剤での治療期間に差がある場合、その組織型に指定薬剤が有効ないし無効である可能性が示唆されます。
</details>
  
<details>
<summary>薬剤の効果と遺伝子変異の関係をvolcano plotで網羅的に確認したい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢、治療コースなどの絞り込みを行う<br>  
3. Analysis -> Drug response analysis -> List of drugs used in Palliative CTxボタンを押す<br>  
4. Results -> Drug response -> Tables -> Drug use, by line of treatmentから使用状況を確認する<br>  
5. Analysis -> Drug response analysis -> Choose drugs for treatment effect analysisで治療を選択する<br>
6. Analysis -> Drug response analysis -> Analyze with the setting selected aboveボタンで解析を行う<br>
7. Results -> Drug response -> Response rate -> Volcano plot for objective response rateで奏効性に関連する遺伝子変異を探索する<br>
<br>
右上の赤い遺伝子では変異があると奏効率が高く、左上の青い遺伝子では変異があると奏効率が低くなります。
</details>
  

  
### 生存期間解析
<details>
<summary>CGP検査後の予後と遺伝子変異の関係をみたい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢、治療コースなどの絞り込みを行う。とくにGenes of interestで注目する遺伝子セットを指定する。<br>  
3. Analysis -> Survival analysis after CGP test ボタンを押す<br>  
4. Results -> Survival after CGP -> Survival analysis -> Survival after CGP and performance statusから指定遺伝子セットないのいずれかに変異があるか否かで群分けした生存曲線を確認する。<br>
5. Results -> Survival after CGP -> Survival analysis -> Survival after CGP and mutations, forest plotで、変異頻度の高い遺伝子について、変異の有無での2群間での生存期間の比較を行う<br>
6. Results -> Survival after CGP -> Survival analysis -> Survival after CGP and mutations, KM-curveで、変異頻度の高い遺伝子について、変異の有無での2群間での生存曲線の比較を行う<br>
<br>
SettingのTiming for RMST measuring in survival analysis (years) で、forest plotで描画する生存期間(restricted mean survival time)の差を計算する時期を指定します。
</details>

  
<details>
<summary>緩和的化学療法導入後の予後と遺伝子変異の関係をみたい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢、治療コースなどの絞り込みを行う。とくにGenes of interestで注目する遺伝子セットを指定する。<br>  
3. Analysis -> Survival analysis after CTx induction ボタンを押す<br>  
4. Results -> Survival after CTx -> Genetic variants and survival, forest plotで、変異頻度の高い遺伝子について、変異の有無での2群間での生存期間の比較を行う<br>
5. Results -> Survival after CTx -> Genetic variants and survival, KM-curveで、変異頻度の高い遺伝子について、変異の有無での2群間での生存曲線の比較を行う<br>
<br>
左側切断バイアスを補正した場合としない場合で生存曲線が描かれます。<br>
</details>

  
<details>
<summary>ある薬剤を生存中に使用したか否かでの生存期間の差をみたい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢、治療コースなどの絞り込みを行う<br>  
3. Analysis -> Drug response analysis -> List of drugs used in Palliative CTxボタンを押す<br>  
4. Results -> Drug response -> Tables -> Drug use, by line of treatmentから使用状況を確認する<br>  
5. Analysis -> Drug response analysis -> Choose drugs for treatment effect analysisで治療を選択する<br>
6. Analysis -> Drug response analysis -> Analyze with the setting selected aboveボタンで解析を行う<br>
7. Results -> Drug response -> Survival after CGP -> Survival and drugで、緩和的化学療法導入後に指定したレジメンを使用したか否かでの2群で分けた生存曲線を確認する<br>
<br>
左側切断バイアスを補正した場合としない場合で生存曲線が描かれます。<br>
</details>

  
<details>
<summary>ある薬剤をCGP検査後に使用したか否かでの生存期間の差をみたい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢、治療コースなどの絞り込みを行う<br>  
3. Analysis -> Drug response analysis -> List of drugs used in Palliative CTxボタンを押す<br>  
4. Results -> Drug response -> Tables -> Drug use, by line of treatmentから使用状況を確認する<br>  
5. Analysis -> Drug response analysis -> Choose drugs for treatment effect analysisで治療を選択する<br>
6. Analysis -> Drug response analysis -> Analyze with the setting selected aboveボタンで解析を行う<br>
7. Results -> Drug response -> Survival after CGP -> Survival and drugで、CGP検査後に指定したレジメンを使用したか否か、そして治療を受けなかった群で2〜3群に分けた生存曲線を確認する<br>
<br>
通常のカプラン・マイアー生存曲線が描かれます。<br>
</details>

  
<details>
<summary>CGP検査後の死亡ハザードに関係する因子を抽出したい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢、治療コースなどの絞り込みを行う。<br>  
3. Analysis -> Survival analysis after CGP test ボタンを押す<br>  
4. Results -> Survival after CGP -> Survival analysis -> Hazard ratio for survival after CGP -genes から、単変量解析・多変量解析でのハザード比に関係する臨床情報や遺伝子変異を検討する。<br>
5. Results -> Survival after CGP -> Survival analysis -> Hazard ratio for survival after CGP -genes から、単変量解析・多変量解析でのハザード比に関係する臨床情報や遺伝子パターン（クラスタリング）を検討する。<br>
<br>
赤池情報量規準を用いて自動的に変数選択を行っています。
</details>
  

  
### 患者背景
<details>
<summary>患者背景のTableが欲しい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢の絞り込みを行う<br>  
3. AnalysisからCase summaryボタンを押す<br>  
4. Results -> Case summaryから結果を確認する<br>  
5. 必要があれば全体をコピーしてWordに保存する<br>
<br>
SettingのFilters on mutation types で選択した遺伝子変異の有無で群分けして表示されます。
</details>
  
<details>
<summary>患者背景のFigureが欲しい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢の絞り込みを行う<br>  
3. AnalysisからClustering based on variantsボタンを押す<br>  
4. Results -> Clustering analysis -> Basic dataから結果を確認する<br>  
5. 必要があれば全体をコピーしてWordに保存する<br>  
<br>
Driverの項目は、何らかのがん化変異(C-CAT evidence level "F")が検出された症例か否かを示します。<br>
Pts with recommended CTxはエキスパートパネルで推奨治療があった症例の割合を意味します。<br>
Pts received recommended CTxは推奨治療を実際に受けた症例の割合を意味します。<br>
Median time from CTx to CGPは緩和的化学療法開始日からCGP検査日までの期間の中央値を意味します。<br>
Median time from CGP to deathはCGP検査日から死亡までの期間のKaplan-Meier法での中央値を意味します。<br> 
</details>

<details>
<summary>オンコプリント（遺伝子変異の一覧表）が欲しい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢の絞り込みを行う<br>  
3. AnalysisからOncoprintボタンを押す<br>  
4. Results -> Oncoprint -> Figures -> Oncoprintから結果を確認する<br>  
5. AnalysisからMutation rate of each gene for each histologyボタンを押す<br>  
6. Results -> Variation by histologyから組織型ごとにどの遺伝子変異の頻度が高いのかを確認する<br>  
<br>
描画した元データはDownloadable tableからExcelファイルでダウンロード可能です。<br> 
</details>

<details>
<summary>ロリプロット（遺伝子のどのアミノ酸残基に変異が多いかの図）が欲しい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢の絞り込みを行う。3列目のGene for lolliplotから遺伝子を指定する。<br>  
3. AnalysisからOncoprintボタンを押す<br>  
4. Results -> Oncoprint -> Figures -> Lolliplot for the selected geneから結果を確認する<br>  
<br>
描画した元データはDownloadable tableからExcelファイルでダウンロード可能です。<br> 
</details>
  
<details>
<summary>遺伝子間の相互排他性・共変異の情報が欲しい</summary>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢の絞り込みを行う<br>  
3. AnalysisからMutually exclusive or co-occurring mutationボタンを押す<br>  
4. Results -> Mutually exclusivityから結果を確認する<br>  
<br>
X軸の遺伝子とY軸の遺伝子の交わるセルの色が青いと両者は相互排他的、赤いと共変異の関係です。<br> 
</details>


  
### データのキュレーション
<details>
<summary>組織型のキュレーションをしたい</summary>
2019年ころの症例を中心にして、詳細な組織型が登録されていない場合があります。<br>
各病院の担当者が入力した手入力の情報を基にして再分類することが可能です。<br>
1. Input C-CAT filesからcase/report CSVファイルを取り込む<br>  
2. Settingから組織型や年齢の絞り込みを行う<br>  
3. AnalysisからOncoprintボタンを押す<br>  
4. Oncoprint -> Downloadable tableの左上のExcelボタンから結果をダウンロードする<br>  
5. P列（病理診断名）、Q列（臨床診断名）、R列（提出検体の病理診断名）を参考に、S列（がん種.OncoTree.）を修正する<br>
6. Input C-CAT files -> Correspondence table between ID and histology (CSV) -> Download CSV file templateボタンを押し保存する<br>  
7. 5で作成した表のハッシュID列とがん種.OncoTree.列の内容をID列とHistology列に貼り付ける<br>  
8. Input C-CAT files -> Correspondence table between ID and histology (CSV)から作成したCSVファイルを取り込む<br>  
<br>
「がん種.OncoTree.」の記載と「がん種.OncoTree.LEVEL1.」の記載が同じ症例だけキュレーションすると労力が少ないと思います。<br>
</details>

  
### 左側切断バイアスのシミュレーション
<details>
<summary>左側切断バイアスについて確認したい</summary>
C-CATのデータのように、生存期間の測定開始日と検査日（観察開始日）が異なる場合、通常のカプラン・マイアー法では生存期間の推定が困難です。<br>
生存期間の測定開始日から検査日まで、全症例が生存している、Immortal biasが存在しているからです。<br>
生存期間の測定開始日から検査日までの生存期間と、検査日から最終観察日までの生存期間に分割すると、バイアスの一部が解消されます。<br>
ただし、「CGP検査を受けた患者は受けなかった患者と何が違うのか」は究極的にはわからず、ある程度の選択バイアスの解消は不可能と考えます。<br>
1. Results -> Bias correction simulation -> An example of bias adjustmentを開く。<br>  
2. お好みに応じてパラメタを調整する<br>  
3. Left-truncation bias adjustment simulationボタンを押す<br>  
4. 真の生存曲線、通常のカプラン・マイアー生存曲線、CGP検査前後で分割した生存曲線、ベイズ推定でのバイアス解消手法による生存曲線が描画されます。<br>  
<br>
概ね良好なバイアスの補正ができているのではないでしょうか。<br>
</details>
  


