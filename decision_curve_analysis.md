## CGP検査に対するDecision curve analysisについて考えたこと
Decision curve analysisは、結果を予測するアルゴリズムの有用性を評価する一つの手法として開発されました。  
[Wikipedia](https://en.wikipedia.org/wiki/Decision_curve_analysis)  
そのR言語をはじめとしたプログラミング言語での実装および理論の詳細な説明がMSKCCからなされています。  
[Decision curve analysis software tutorial](https://mskcc-epi-bio.github.io/decisioncurveanalysis/index.html)
患者さんや主治医の意向、すなわちどの程度の


### 解析の限界
C-CATに登録された患者は当然ながら全例「CGP検査を受けた」というバイアスがあります。主治医や患者さんの漠然とした選択基準に従いCGP検査が行われています。
一方で、その基準に従ってCGP検査を受けなかった患者さんの情報はC−CATから得ることはできず、畢竟新たに開発したノモグラムの外的妥当性がどの程度あるのかは不明確です。
現状のノモグラムおよびDecision curveの活用としては、「臨床的にCGP検査の適応を考慮した患者さんに本当に検査を実施すべきか判断する」というシチュエーションになろうかと思われます。  
