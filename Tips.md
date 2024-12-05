## FELIS 解析例 <img src="source/FELIS.png" width=50>
必要性の高いであろう解析について方法をまとめてみました。

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

