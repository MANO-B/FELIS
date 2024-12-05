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
7. Results -> Drug response -> Time on treatment -> Time on treatment by tissue type, KM-curveで、全ての治療あるいは指定治療のTime on treatmentをKaplan-Meier法で評価する<br>
<br>
全ての薬剤での治療期間と指定薬剤での治療期間に差がある場合、その組織型に指定薬剤が有効ないし無効である可能性が示唆されます。
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


