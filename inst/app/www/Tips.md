# 💡 FELIS 解析Tips集
必要性の高いであろう解析について方法をまとめてみました。各項目をクリックして詳細な手順をご確認ください。

---

## 💊 薬剤奏効性解析

<details>
<summary><b>ある治療と、一つ前のコースの治療との効果比較をしたい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う
3. **Analysis** -> **Drug response analysis** -> **List of drugs used in Palliative CTx** ボタンを押す
4. **Results** -> **Drug response** -> **Tables** -> **Drug use, by line of treatment** から使用状況を確認する
5. **Analysis** -> **Drug response analysis** -> **Choose drugs for treatment effect analysis** で治療を選択する
6. **Analysis** -> **Drug response analysis** -> **Analyze with the setting selected above** ボタンで解析を行う
7. **Results** -> **Drug response** -> **Time on treatment** -> **Time on treatment and pre-treatment for the specified treatment, scatter plot** で、同一患者群での前治療と指定治療の time on treatment の比較を行う
8. **Results** -> **Drug response** -> **Time on treatment** -> **Time on treatment and pre-treatment for the specified treatment, KM-curve** で、同一患者群での前治療と指定治療の time on treatment の比較を行う

</details>

<details>
<summary><b>二つの治療の治療効果比較をしたい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う
3. **Analysis** -> **Drug response analysis** -> **List of drugs used in Palliative CTx** ボタンを押す
4. **Results** -> **Drug response** -> **Tables** -> **Drug use, by line of treatment** から使用状況を確認する
5. **Analysis** -> **Drug response analysis** -> **Choose drugs for treatment effect analysis** で二つ以上の治療を選択する。比較したい2レジメンを、それぞれ **Drug set 1** と **Drug set 2** に入力する
6. **Analysis** -> **Drug response analysis** -> **Analyze with the setting selected above** ボタンで解析を行う
7. **Results** -> **Drug response** -> **Survival after CGP** -> **Survival and drug** で、緩和的化学療法導入後に指定したレジメンで2群で分けた生存曲線を確認する

> 💡 **Tip:**
> とくに、1st lineの治療に限定すると、例えば膵がんの1st line治療としてFOLFIRINOXとGEM + nab-PTXのどちらが生存期間が優れるかの比較ができたりします。
> 発売が新しい薬剤では左側切断バイアスが強く出るため、治療開始日でマッチングさせた解析も加えました。

</details>

<details>
<summary><b>組織型ごとの治療効果をみたい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う
3. **Analysis** -> **Drug response analysis** -> **List of drugs used in Palliative CTx** ボタンを押す
4. **Results** -> **Drug response** -> **Tables** -> **Drug use, by line of treatment** から使用状況を確認する
5. **Analysis** -> **Drug response analysis** -> **Choose drugs for treatment effect analysis** で治療を選択する
6. **Analysis** -> **Drug response analysis** -> **Analyze with the setting selected above** ボタンで解析を行う
7. **Results** -> **Drug response** -> **Time on treatment** -> **Time on treatment and pre-treatment for the specified treatment, KM-curve** で、遺伝子変異の有無で群分けした指定治療の Time on treatment を Kaplan-Meier 法で評価する

> 💡 **解釈のポイント:**
> 全ての薬剤での治療期間と指定薬剤での治療期間に差がある場合、その遺伝子変異が指定薬剤の biomarker である可能性が示唆されます。

</details>

<details>
<summary><b>ある薬剤の効果における遺伝子変異の有無の意義をみたい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込み、探索したい遺伝子の指定を行う
3. **Analysis** -> **Drug response analysis** -> **List of drugs used in Palliative CTx** ボタンを押す
4. **Results** -> **Drug response** -> **Tables** -> **Drug use, by line of treatment** から使用状況を確認する
5. **Analysis** -> **Drug response analysis** -> **Choose drugs for treatment effect analysis** で治療を選択する
6. **Analysis** -> **Drug response analysis** -> **Analyze with the setting selected above** ボタンで解析を行う
7. **Results** -> **Drug response** -> **Time on treatment** -> **Time on treatment by tissue type, KM-curve** で、全ての治療あるいは指定治療の Time on treatment を Kaplan-Meier 法で評価する

> 💡 **解釈のポイント:**
> 全ての薬剤での治療期間と指定薬剤での治療期間に差がある場合、その組織型に指定薬剤が有効ないし無効である可能性が示唆されます。

</details>

<details>
<summary><b>薬剤の効果と遺伝子変異の関係を volcano plot で網羅的に確認したい</b></summary>



1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う
3. **Analysis** -> **Drug response analysis** -> **List of drugs used in Palliative CTx** ボタンを押す
4. **Results** -> **Drug response** -> **Tables** -> **Drug use, by line of treatment** から使用状況を確認する
5. **Analysis** -> **Drug response analysis** -> **Choose drugs for treatment effect analysis** で治療を選択する
6. **Analysis** -> **Drug response analysis** -> **Analyze with the setting selected above** ボタンで解析を行う
7. **Results** -> **Drug response** -> **Response rate** -> **Volcano plot for objective response rate** で奏効性に関連する遺伝子変異を探索する

> 💡 **見方:**
> 右上の赤い遺伝子では変異があると奏効率が高く、左上の青い遺伝子では変異があると奏効率が低くなります。

</details>

---

## 📈 生存期間解析

<details>
<summary><b>CGP検査後の予後と遺伝子変異の関係をみたい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う。とくに **Genes of interest** で注目する遺伝子セットを指定する。
3. **Analysis** -> **Survival analysis after CGP test** ボタンを押す
4. **Results** -> **Survival after CGP** -> **Survival analysis** -> **Survival after CGP and performance status** から指定遺伝子セットないのいずれかに変異があるか否かで群分けした生存曲線を確認する。
5. **Results** -> **Survival after CGP** -> **Survival analysis** -> **Survival after CGP and mutations, forest plot** で、変異頻度の高い遺伝子について、変異の有無での2群間での生存期間の比較を行う
6. **Results** -> **Survival after CGP** -> **Survival analysis** -> **Survival after CGP and mutations, KM-curve** で、変異頻度の高い遺伝子について、変異の有無での2群間での生存曲線の比較を行う

> 📝 **補足:**
> **Setting** の `Timing for RMST measuring in survival analysis (years)` で、forest plot で描画する生存期間 (restricted mean survival time) の差を計算する時期を指定します。

</details>

<details>
<summary><b>緩和的化学療法導入後の予後と遺伝子変異の関係をみたい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う。とくに **Genes of interest** で注目する遺伝子セットを指定する。
3. **Analysis** -> **Survival analysis after CTx induction** ボタンを押す
4. **Results** -> **Survival after CTx** -> **Genetic variants and survival, forest plot** で、変異頻度の高い遺伝子について、変異の有無での2群間での生存期間の比較を行う
5. **Results** -> **Survival after CTx** -> **Genetic variants and survival, KM-curve** で、変異頻度の高い遺伝子について、変異の有無での2群間での生存曲線の比較を行う

> 📝 **補足:**
> 左側切断バイアスを補正した場合としない場合で生存曲線が描かれます。

</details>

<details>
<summary><b>ある薬剤を生存中に使用したか否かでの生存期間の差をみたい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う
3. **Analysis** -> **Drug response analysis** -> **List of drugs used in Palliative CTx** ボタンを押す
4. **Results** -> **Drug response** -> **Tables** -> **Drug use, by line of treatment** から使用状況を確認する
5. **Analysis** -> **Drug response analysis** -> **Choose drugs for treatment effect analysis** で治療を選択する
6. **Analysis** -> **Drug response analysis** -> **Analyze with the setting selected above** ボタンで解析を行う
7. **Results** -> **Drug response** -> **Survival after CGP** -> **Survival and drug** で、緩和的化学療法導入後に指定したレジメンを使用したか否かでの2群で分けた生存曲線を確認する

> 📝 **補足:**
> 左側切断バイアスを補正した場合としない場合で生存曲線が描かれます。

</details>

<details>
<summary><b>ある薬剤をCGP検査後に使用したか否かでの生存期間の差をみたい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う
3. **Analysis** -> **Drug response analysis** -> **List of drugs used in Palliative CTx** ボタンを押す
4. **Results** -> **Drug response** -> **Tables** -> **Drug use, by line of treatment** から使用状況を確認する
5. **Analysis** -> **Drug response analysis** -> **Choose drugs for treatment effect analysis** で治療を選択する
6. **Analysis** -> **Drug response analysis** -> **Analyze with the setting selected above** ボタンで解析を行う
7. **Results** -> **Drug response** -> **Survival after CGP** -> **Survival and drug** で、CGP検査後に指定したレジメンを使用したか否か、そして治療を受けなかった群で2〜3群に分けた生存曲線を確認する

> 📝 **補足:**
> 通常のカプラン・マイアー生存曲線が描かれます。

</details>

<details>
<summary><b>CGP検査後の死亡ハザードに関係する因子を抽出したい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う。
3. **Analysis** -> **Survival analysis after CGP test** ボタンを押す
4. **Results** -> **Survival after CGP** -> **Survival analysis** -> **Hazard ratio for survival after CGP -genes** から、単変量解析・多変量解析でのハザード比に関係する臨床情報や遺伝子変異を検討する。
5. **Results** -> **Survival after CGP** -> **Survival analysis** -> **Hazard ratio for survival after CGP -clusters** （※UIパス名に準拠）から、単変量解析・多変量解析でのハザード比に関係する臨床情報や遺伝子パターン（クラスタリング）を検討する。

> 📝 **補足:**
> 赤池情報量規準（AIC）を用いて自動的に変数選択を行っています。

</details>

---

## 👥 患者背景

<details>
<summary><b>患者背景のTableが欲しい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢の絞り込みを行う
3. **Analysis** から **Case summary** ボタンを押す
4. **Results** -> **Case summary** から結果を確認する
5. 必要があれば全体をコピーしてWordに保存する

> 📝 **補足:**
> **Setting** の `Filters on mutation types` で選択した遺伝子変異の有無で群分けして表示されます。

</details>

<details>
<summary><b>患者背景のFigureが欲しい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢の絞り込みを行う
3. **Analysis** から **Clustering based on variants** ボタンを押す
4. **Results** -> **Clustering analysis** -> **Basic data** から結果を確認する
5. 必要があれば全体をコピーしてWordに保存する

> 📖 **項目の意味:**
> * **Driver**: 何らかのがん化変異 (C-CAT evidence level "F") が検出された症例か否かを示します。
> * **Pts with recommended CTx**: エキスパートパネルで推奨治療があった症例の割合。
> * **Pts received recommended CTx**: 推奨治療を実際に受けた症例の割合。
> * **Median time from CTx to CGP**: 緩和的化学療法開始日からCGP検査日までの期間の中央値。
> * **Median time from CGP to death**: CGP検査日から死亡までの期間のKaplan-Meier法での中央値。

</details>

<details>
<summary><b>オンコプリント（遺伝子変異の一覧表）が欲しい</b></summary>



1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢の絞り込みを行う
3. **Analysis** から **Oncoprint** ボタンを押す
4. **Results** -> **Oncoprint** -> **Figures** -> **Oncoprint** から結果を確認する
5. **Analysis** から **Mutation rate of each gene for each histology** ボタンを押す
6. **Results** -> **Variation by histology** から組織型ごとにどの遺伝子変異の頻度が高いのかを確認する

> 📝 **補足:**
> 描画した元データは **Downloadable table** からExcelファイルでダウンロード可能です。

</details>

<details>
<summary><b>ロリプロット（遺伝子のどのアミノ酸残基に変異が多いかの図）が欲しい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢の絞り込みを行う。3列目の `Gene for lolliplot` から遺伝子を指定する。
3. **Analysis** から **Oncoprint** ボタンを押す
4. **Results** -> **Oncoprint** -> **Figures** -> **Lolliplot for the selected gene** から結果を確認する

> 📝 **補足:**
> 描画した元データは **Downloadable table** からExcelファイルでダウンロード可能です。
> 現状ではエキソンスキッピングやイントロンの変異には対応していません。

</details>

<details>
<summary><b>遺伝子間の相互排他性・共変異の情報が欲しい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢の絞り込みを行う
3. **Analysis** から **Mutually exclusive or co-occurring mutation** ボタンを押す
4. **Results** -> **Mutually exclusivity** から結果を確認する

> 💡 **見方:**
> X軸の遺伝子とY軸の遺伝子の交わるセルの色が**青い**と両者は相互排他的、**赤い**と共変異の関係です。

</details>

---

## 🎯 治療到達性解析

<details>
<summary><b>どのような患者にCGP検査を行うと治療到達率が高いのかを知りたい</b></summary>



1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う。とくに **Genes of interest** で注目する遺伝子セットを指定する。
3. **Analysis** -> **CGP benefit prediction analysis** ボタンを押す
4. **Results** -> **CGP benefit prediction** -> **Factors leading to treatment** -> **Factors leading to treatment, pre-CGP, Nomogram** から、検査前に得られる患者の臨床情報に基づいて治療到達率を予測するノモグラムを得る。
5. **Results** -> **CGP benefit prediction** -> **Factors leading to treatment** -> **Factors leading to treatment, pre-CGP, Odds ratio** から、検査前に得られる患者の臨床情報が治療到達率に与える影響を単変量・多変量で解析する。
6. **Results** -> **CGP benefit prediction** -> **Factors leading to treatment** -> **ROC curve of nomogram** から、ノモグラムによる予測の精度をROC曲線で確認する。
7. **Results** -> **CGP benefit prediction** -> **Factors leading to treatment** -> **Factors leading to treatment, decision curve** から、ノモグラムによる予測の臨床的有用性をdecision curve analysisで評価した結果を確認する。
8. **Results** -> **CGP benefit prediction** -> **Factors leading to treatment** -> **Analyze your data** から、特定の患者さんの情報を入力すると治療到達率が予想される。

> 📖 **参考:**
> Decision curve analysis については <a href="https://mskcc-epi-bio.github.io/decisioncurveanalysis/index.html" target="_blank">こちら</a> や <a href="https://github.com/MANO-B/FELIS/blob/main/decision_curve_analysis.md" target="_blank">こちら</a> を参照下さい。

</details>

<details>
<summary><b>CGP検査後の生存期間と治療到達率の関係性を知りたい</b></summary>

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢、治療コースなどの絞り込みを行う。とくに **Genes of interest** で注目する遺伝子セットを指定する。
3. **Analysis** -> **Survival analysis after CGP test** ボタンを押す
4. **Results** -> **Survival after CGP** -> **Survival analysis** -> **Survival period and treatment reach rate** から、CGP検査後の生存期間と治療到達率について移動平均を取ったグラフを確認する。

> 💡 **考察:**
> CGP検査後に短期で死亡する患者は治療到達率が低いため、いかに予後が悪そうな患者さんに早めに検査を行うか、そしてPSが保たれ一定程度の生存期間がある患者さんに検査を行うかが重要と思われます。

</details>

---

## 🔧 データのキュレーション

<details>
<summary><b>組織型のキュレーションをしたい</b></summary>

2019年ころの症例を中心にして、詳細な組織型が登録されていない場合があります。各病院の担当者が入力した手入力の情報を基にして再分類することが可能です。

1. **Input C-CAT files** から case/report CSVファイルを取り込む
2. **Setting** から組織型や年齢の絞り込みを行う
3. **Analysis** から **Oncoprint** ボタンを押す
4. **Oncoprint** -> **Downloadable table** の左上のExcelボタンから結果をダウンロードする
5. P列（病理診断名）、Q列（臨床診断名）、R列（提出検体の病理診断名）を参考に、S列（がん種.OncoTree.）を修正する
6. **Input C-CAT files** -> **Correspondence table between ID and histology (CSV)** -> **Download CSV file template** ボタンを押し保存する
7. 5で作成した表のハッシュID列とがん種.OncoTree.列の内容をID列とHistology列に貼り付ける
8. **Input C-CAT files** -> **Correspondence table between ID and histology (CSV)** から作成したCSVファイルを取り込む

> 💡 **Tip:**
> 「がん種.OncoTree.」の記載と「がん種.OncoTree.LEVEL1.」の記載が同じ症例だけキュレーションすると労力が少ないと思います。

</details>

---

## ⏳ 左側切断バイアスのシミュレーション

<details>
<summary><b>左側切断バイアスについて確認したい</b></summary>



C-CATのデータのように、生存期間の測定開始日と検査日（観察開始日）が異なる場合、通常のカプラン・マイアー法では生存期間の推定が困難です。生存期間の測定開始日から検査日まで、全症例が生存している、Immortal biasが存在しているからです。
生存期間の測定開始日から検査日までの生存期間と、検査日から最終観察日までの生存期間に分割すると、バイアスの一部が解消されます。
ただし、「CGP検査を受けた患者は受けなかった患者と何が違うのか」は究極的にはわからず、ある程度の選択バイアスの解消は不可能と考えます。

1. **Results** -> **Bias correction simulation** -> **An example of bias adjustment** を開く。
2. お好みに応じてパラメタを調整する
3. **Left-truncation bias adjustment simulation** ボタンを押す
4. 真の生存曲線、通常のカプラン・マイアー生存曲線、CGP検査前後で分割した生存曲線、ベイズ推定でのバイアス解消手法による生存曲線が描画されます。

> 💡 **解釈:**
> 概ね良好なバイアスの補正ができているのではないでしょうか。

</details>

---

## 💾 その他

<details>
<summary><b>図を保存したい</b></summary>

図を右クリックし、拡張子を `.png` として名前をつけて保存して下さい。

</details>

<details>
<summary><b>サマリーの表を保存したい</b></summary>

画面上でテキストを選択し、コピーしてWordかExcelに貼り付けて保存して下さい。

</details>

<details>
<summary><b>生データの表を保存したい</b></summary>

左上にあるボタンでエクセルファイルあるいはCSVファイルなどで名前をつけて保存して下さい。

</details>
