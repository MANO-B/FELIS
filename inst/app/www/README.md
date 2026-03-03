# FELIS: Flexible Exploration for LIquid and Solid tumor clinical sequencing data
### — C-CAT Secondary Use Data Analysis Platform —

FELIS（Flexible Exploration for LIquid and Solid tumor clinical sequencing data）は、C-CAT（二次利用）データを対象に、GUI（R/Shiny）を通じてコホート定義から可視化、そして臨床的に重要なバイアスを考慮したアウトカム解析までを一通り行うためのローカル実行型Webアプリです。

> **注意**: 本ソフトは C-CAT二次利用データを適法に取得できる方のみが対象です。各施設の倫理審査・データ利用規約等に従って利用してください。

---

## 1. 概要 (Overview)
C-CATは、日本のがんゲノム医療で実施されるCGP検査結果と臨床情報が集約されるナショナルDBです。
FELISは、その二次利用データを対象に、以下のプロセスをノーコード/低コードで反復探索できることを狙っています。
1. **コホート構築**
2. **変異の要約と可視化**
3. **治療・予後の解析**

特に、CGPのリアルワールドデータで問題になりやすい **遅延到達（delayed entry）** や **左側切断（left truncation）** を意識した解析導線をGUIとして統合している点が特徴です。

---

## 2. GUIの全体構造 (Menu Mapping)

- **Settings**: 解析対象（症例・組織型・遺伝子・治療など）の指定、しきい値、モデル設定など。
- **Results**: 各解析の図表（保存/ダウンロード含む）。
- **Instruction / Tips**: アプリ内で README / Tips を表示（Instruction タブは `README.md` を表示）。

---

## 3. Results：出力の意味（網羅版リファレンス）

以下は、FELISのメニュー構造に沿った各アウトカムの解説です。

### ■ Case summary
* **Summarized by mutation pattern**
  患者の基本属性や臨床変数を、指定した変異パターン等で要約したテーブル。
* **Summarized by histology**
  組織型（OncoTree等）単位での症例数、臨床情報、変異要約。

### ■ Oncoprint（図・テーブル）
* **Oncoprint**
  選択コホートで頻度の高い遺伝子の変異景観（患者×遺伝子）を可視化。
* **Lolliplot for the selected gene**
  指定遺伝子のアミノ酸変化（ホットスポット等）の頻度分布。
* **Table of clinical and mutation information per patient**
  患者単位の臨床情報と変異情報のテーブル。解析結果の元データ確認・二次解析用。

### ■ Mutual exclusivity（共起/排他）
* **Figure for probability / odds ratio**
  変異ペアの共起/排他傾向を、確率・オッズ比などで可視化（青＝排他、赤＝共起）。
* **table_mutually_exclusive**
  ペアごとの統計量（p値等）を一覧。

### ■ Variant rate by histology（組織型別頻度）
* **figure_mut_subtype_plot**
  組織型サブタイプ別の頻度上位遺伝子の変異率を比較。

### ■ Mutation and treatment option（変異クラスタリング等）
* **Basic data**
  組織型ごとの基本統計（年齢、性別、変異、TMB、治療オプション/実治療、CGP前後の期間など）を図表化。
* **UMAP clustering based on mutations / Cluster and histology relationship**
  変異パターンに基づく患者クラスタ（UMAP + DBSCAN等）を作成し、クラスタと組織型・変異の富化を表示。
* **Heterogeneity within histologic types**
  組織型ごとのクラスタ分布（集中/分散）をエントロピー等で評価。
* **Frequency of patients with targeted therapy**
  各患者におけるエビデンスレベル別の治療オプション頻度を集計。

### ■ Survival after CGP（CGP後生存）
* **Survival and clinical information**
  CGP検査日を起点にした生存曲線（KM等）と、層別（組織型・PS・治療など）。
* **Custom survival analysis**
  2群比較をGUIで定義し、KM/RMST等を出力。PSM/IPWなどの背景調整（Love plot、重み分布、PS分布等）も提供。
* **Survival and mutations, forest plot**
  変異有無等の層別でRMST差/推定効果をフォレストプロットで表示。
* **Hazard ratio**
  Coxモデル等での単変量/多変量の推定結果（gtsummary形式）。
* **Survival period and treatment reach rate**
  “治療開始”をアウトカムとし、死亡を競合リスクとして扱うCIF（累積発生関数）等の到達解析。

### ■ CGP benefit prediction（到達予測・DCA/ROC）
* **Nomogram**
  CGP前情報から“治療到達”を予測するノモグラムを表示。
* **Odds ratio**
  到達に関連する因子のオッズ比（単変量/多変量）。
* **Decision curve**
  DCA（Decision Curve Analysis）で、予測モデルを使う臨床的価値（Net benefit）を評価。
* **ROC curve of nomogram**
  ROC（AUC等）による予測精度の評価。
* **Input data（Input_data）**
  GUIから臨床変数を手入力し、上記モデルで予測値を返すシミュレーション。

### ■ Overall survival with risk-set adjustment（左側切断補正手法：risk-set）
* **Survival and clinical information（SurvivalandtreatmentafterCTx）**
  生存の起点を「化学療法開始」等に置いた時の、左側切断（delayed entry）を補正した推定。
* **Custom survival analysis**
  2群比較におけるRisk-set補正後の比較。
* **Frequent variants and survival / Diagnosis and survival / Mutational cluster and survival**
  変異、診断、クラスタなどを説明変数にした生存解析の図（フォレスト/曲線）。
* **Hazard ratio for survival after CTx**
  補正枠組み下でのCoxモデル等の推定結果テーブル。

### ■ Survival after CTx with Bayesian inference（左側切断補正手法：Bayesian）
* **Survival corrected for left-truncation bias**
  Stanを用いたベイズ推定により、左側切断を考慮した生存曲線の中央値と信頼区間を出力。
* **Custom survival analysis**
  ベイズ枠組みでの2群間比較。
* **Genetic variants and survival / Diagnosis and survival**
  変異や診断で層別した、補正後の生存曲線や効果推定。

### ■ Survival after CTx with control cohort data（experimental）
* **Custom survival analysis**
  院内がん登録を対照コホート情報として用いた実験的解析手法。現状はSettingsでStage4患者のみを選択する必要があり、また診断後の生存期間解析に限定。

### ■ Bias correction simulation（シミュレーション）
* **Left-truncation bias adjustment simulation**
  パラメータを与えて、左側切断補正の挙動を視覚的に理解するためのシミュレーター。

### ■ Drug response（ToT/ORR/AE/薬剤別生存）
* **Settings**
  解析対象の投与ライン、解析期間（Time on treatment (ToT)など）、対象薬剤群（まとめ上げ）などを指定。
* **Drug usage data**
  患者×薬剤使用の元データテーブル（ダウンロード可能）。
* **Treatment time and clinical information**
  ToTのKM曲線、事前治療との関係（散布図）、フォレストプロット。
* **Treatment time comparison**
  2群定義をGUIで指定し、ToTを比較。
* **ToT by gene mutation cluster / by mutated genes / by mutation pattern**
  クラスタ・遺伝子変異・変異パターンでToTを層別化。
* **Hazard ratio on time on treatment（genes / clusters）**
  ToTをアウトカムとしたCox等の推定結果テーブル。
* **Volcano plot（ToT / ORR / AE）**
  レジメン×遺伝子等の関連を、効果量と有意性でvolcano plotを表示。
* **Cumulative incidence of adverse effect**
  AEをイベントとして競合リスクを意識した累積発生を表示。
* **Survival and drug**
  “投薬開始日”起点の生存を薬剤別に表示。

---

## 4. 統計・疫学的バイアスへの対応 (Methodology)
FELISの設計思想:
1. **セキュア環境**: オフライン環境で動作し、クラウドを使わずに高度な可視化を実現。
2. **選択バイアスの克服**: CGP RWD特有の **delayed entry / left truncation** による生存の過大評価（不死時間バイアス）を、Risk-set補正やベイズシミュレーションで解決。

---

## 5. 引用 (References)
研究成果に使用する際は、以下の文献を引用してください。

```bibtex
@article{Mano2024FELIS,
  title={FELIS: Flexible Exploration for LIquid and Solid tumor clinical sequencing data},
  author={IKEGAMI, Masachika},
  journal={GitHub Repository},
  url={[https://github.com/MANO-B/FELIS](https://github.com/MANO-B/FELIS)},
  year={2024}
}
```
