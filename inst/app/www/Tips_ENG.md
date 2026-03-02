# 💡 FELIS Analysis Tips
We have summarized the methods for analyses that are likely to be highly necessary. Please click on each item to check the detailed procedures.

---

## 💊 Drug response analysis

<details>
<summary><b>I want to compare the efficacy of a certain treatment with the treatment of the previous line</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Select drugs from **Results** -> **Drug response** -> **Settings**
3. Check the usage status from **Results** -> **Drug response** -> **Summary Tables**
4. Compare the time on treatment of the previous treatment and the specified treatment in the same patient group using **Results** -> **Drug response** -> **Time on treatment** -> **Comparison of the treatment duration of the selected drug and the preceding treatment line**
5. Specify **CTx lines in which drugs used**, etc., and compare the time on treatment of the previous treatment and the specified treatment using **Results** -> **Drug response** -> **Time on treatment** -> **Treatment time comparison**

> 💡 **Tip:**
> If the ToT is longer than the treatment of one line prior, that treatment may have high efficacy for the specified patient group.

</details>

<details>
<summary><b>I want to compare the treatment efficacy of two treatments</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Select drugs from **Results** -> **Drug response** -> **Settings**
3. Check the usage status from **Results** -> **Drug response** -> **Summary Tables**
4. Compare the time on treatment of the treatments selected in Settings using **Results** -> **Drug response** -> **Time on treatment** -> **Treatment duration in the selected regimen, grouped by regimen**
5. Specify **CTx lines in which drugs used**, etc., and compare the time on treatment of the two treatments using **Results** -> **Drug response** -> **Time on treatment** -> **Treatment time comparison**
6. Select two or more treatments in **Analysis** -> **Survival after drug initiation date** -> **Choose drugs for treatment effect analysis**. Enter the 2 regimens you want to compare into **Drug set 1** and **Drug set 2**, respectively.
7. In **Analysis** -> **Drug response analysis** -> **Survival and drug**, press the Analyze button after setting, and check the survival curves divided into 2 groups by the specified regimens after the introduction of palliative chemotherapy.

> 💡 **Tip:**
> For example, it is possible to compare whether FOLFIRINOX or GEM + nab-PTX has a superior survival time as a 1st line treatment for pancreatic cancer.

</details>

<details>
<summary><b>I want to see the treatment efficacy for each histology</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Select drugs from **Results** -> **Drug response** -> **Settings**
3. Check the usage status from **Results** -> **Drug response** -> **Summary Tables**
4. Compare the time on treatment of the treatments selected in Settings using **Results** -> **Drug response** -> **Time on treatment** -> **Treatment duration in the selected regimen, grouped by regimen**
5. Evaluate the Time on treatment of all treatments grouped by the presence or absence of gene mutations using the Kaplan-Meier method in **Results** -> **Drug response** -> **Time on treatment** -> **Treatment duration of all drugs in the selected line, grouped by diagnosis**
6. Evaluate the Time on treatment of the specified treatment grouped by the presence or absence of gene mutations using the Kaplan-Meier method in **Results** -> **Drug response** -> **Time on treatment** -> **Time on treatment and pre-treatment for the specified treatment, KM-curve**

> 💡 **Interpretation Point:**
> If the treatment period does not change between histologies, and the treatment period for the specified drug is long in a specific histology, it suggests the possibility that the specified drug is highly effective for that histology.

</details>

<details>
<summary><b>I want to see the significance of the presence or absence of a gene mutation on the effect of a certain drug</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Select drugs from **Results** -> **Drug response** -> **Settings**
3. Check the usage status from **Results** -> **Drug response** -> **Summary Tables**
4. Compare the time on treatment of the treatments selected in Settings using **Results** -> **Drug response** -> **Time on treatment** -> **Treatment duration in the selected regimen, grouped by regimen**
5. Evaluate the Time on treatment of all treatments or the specified treatment using the Kaplan-Meier method in **Results** -> **Drug response** -> **Time on treatment** -> **Treatment duration of all drugs in the selected line, grouped by diagnosis**
6. Evaluate the Time on treatment of the specified treatment grouped by the presence or absence of gene mutations using the Kaplan-Meier method in **Results** -> **Drug response** -> **Time on treatment** -> **Treatment duration in the selected regimen, grouped by detailed mutation**
7. Evaluate the Time on treatment using the Kaplan-Meier method after dividing into 2 groups based on the presence or absence of gene mutations, etc., in **Results** -> **Drug response** -> **Time on treatment** -> **Time on treatment by tissue type, KM-curve**

> 💡 **Interpretation Point:**
> If the treatment period for all drugs does not change depending on the presence or absence of gene mutations, and the treatment period for the specified drug differs depending on the presence or absence of gene mutations, it suggests the possibility that the gene mutation is a biomarker for the specified drug.

</details>

<details>
<summary><b>I want to comprehensively check the relationship between drug effects and gene mutations using a volcano plot</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Select drugs from **Results** -> **Drug response** -> **Settings**
3. Check the usage status from **Results** -> **Drug response** -> **Summary Tables**
4. Explore gene mutations related to Time on treatment in **Results** -> **Drug response** -> **Time on treatment** -> **Volcano plot for treatment time, hazard ratio**
5. Explore gene mutations related to Objective response in **Results** -> **Drug response** -> **Response rate** -> **Volcano plot for objective response rate**

> 💡 **How to read:**
> For light blue genes, the response rate is high when there is a mutation, and for red genes, the response rate is low when there is a mutation. The table can be viewed from the bottom.

</details>

---

## 📈 Survival analysis

<details>
<summary><b>I want to see the relationship between prognosis after CGP testing and gene mutations</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Check the survival curves grouped by whether or not there is a mutation in any of the specified gene set from **Results** -> **Survival after CGP** -> **Survival analysis** -> **Survival after CGP and performance status**.
3. Compare the survival time between the 2 groups with and without mutations for genes with high mutation frequencies in **Results** -> **Survival after CGP** -> **Survival analysis** -> **Survival after CGP and mutations, forest plot**
4. Compare the survival curves between the 2 groups with and without mutations for genes with high mutation frequencies in **Results** -> **Survival after CGP** -> **Survival analysis** -> **Survival after CGP and mutations, KM-curve**
5. Compare the survival curves between 2 groups divided by the presence or absence of mutations, etc., in **Results** -> **Survival after CGP** -> **Survival analysis** -> **Custom survival analysis**

> 📝 **Note:**
> In `Timing for RMST measuring in survival analysis (years)` of **Setting**, specify the time to calculate the difference in survival time (restricted mean survival time) drawn in the forest plot.

</details>

<details>
<summary><b>I want to see the relationship between prognosis after the introduction of palliative chemotherapy and gene mutations</b></summary>

1. Perform filtering such as histology, age, and treatment course from **Setting**. In particular, specify the gene set of interest in **Genes of interest**.
2. Check whether a reasonable survival time is predicted with cKendall tau around 0 in **Results** -> **Overall survival with risk-set adjustment** -> **Survival and clinical information** -> **Entire cohort**
3. Compare survival times by grouping according to clinical information and gene mutation information in **Results** -> **Overall survival with risk-set adjustment** -> **Custom survival analysis**
4. Compare the survival time between the 2 groups with and without mutations for genes with high mutation frequencies in **Results** -> **Survival after CTx with Bayesian inference** -> **Genetic variants and survival, forest plot**
5. Compare the survival curves between the 2 groups with and without mutations for genes with high mutation frequencies in **Results** -> **Survival after CTx with Bayesian inference** -> **Genetic variants and survival, KM-curve**
6. Compare survival times by grouping according to clinical information and gene mutation information in **Results** -> **Survival after CTx with control cohort data (experimental)** -> **Custom survival analysis**

> 📝 **Note:**
> Survival analysis is performed after correcting for left-truncation bias. There are 3 types of correction methods, but the method corrected with external information called hospital-based cancer registries seems to be the most valid (<a href="https://github.com/MANO-B/FELIS/blob/main/paper/新規手法説明書.pdf" target="_blank">here</a>).

</details>

<details>
<summary><b>I want to see the difference in survival time depending on whether a certain drug was used during survival</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Select drugs from **Results** -> **Drug response** -> **Settings**
3. Check the usage status from **Results** -> **Drug response** -> **Summary Tables**
4. Check the survival curves divided into 2 groups by whether or not the specified regimen was used after the introduction of palliative chemotherapy in **Results** -> **Drug response** -> **Survival after CGP** -> **Survival and drug**

> 📝 **Note:**
> Survival curves are drawn for cases where the left-truncation bias is corrected and where it is not.

</details>

<details>
<summary><b>I want to see the difference in survival time depending on whether a certain drug was used after the CGP test</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Compare the survival curves between 2 groups divided by drug usage history, etc., in **Results** -> **Survival after CGP** -> **Survival analysis** -> **Custom survival analysis**

> 📝 **Note:**
> Normal Kaplan-Meier survival curves are drawn.

</details>

<details>
<summary><b>I want to extract factors related to the mortality hazard after the CGP test</b></summary>

1. Perform filtering such as histology and age from **Setting**.
2. Examine clinical information and gene mutations related to the hazard ratio in univariate and multivariate analysis from **Results** -> **Survival after CGP** -> **Survival analysis** -> **Hazard ratio**.

> 📝 **Note:**
> Variable selection is automatically performed by backward elimination using the Akaike Information Criterion (AIC). For multicollinearity, factors with VIF > 10 are excluded.

</details>

---

## 👥 Patient Background

<details>
<summary><b>I want a Table of the patient background</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Check the results from **Results** -> **Case summary**
5. If necessary, copy the whole thing and save it in Word

> 📝 **Note:**
> It is displayed grouped by the presence or absence of gene mutations selected in `Filters on mutation types` of **Setting**.

</details>

<details>
<summary><b>I want a Figure of the patient background</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Check the results from **Results** -> **Clustering analysis** -> **Basic data**
5. If necessary, copy the whole thing and save it in Word

> 📖 **Meaning of items:**
> * **Driver**: Indicates whether the case has any oncogenic mutation (C-CAT evidence level "F") detected.
> * **Pts with recommended CTx**: Means the percentage of cases that had recommended treatments by the expert panel.
> * **Pts received recommended CTx**: Means the percentage of cases that actually received the recommended treatments.
> * **Median time from CTx to CGP**: Means the median time from the start date of palliative chemotherapy to the CGP test date.
> * **Median time from CGP to death**: Means the median time from the CGP test date to death by the Kaplan-Meier method.

</details>

<details>
<summary><b>I want an oncoprint (a list of gene mutations)</b></summary>

1. Perform filtering such as histology and age from **Setting**
2. Check the results from **Results** -> **Oncoprint** -> **Figures** -> **Oncoprint**
3. Download the patient's clinical information and mutation information from **Results** -> **Oncoprint** -> **Downloadable table**
3. Check which gene mutations are highly frequent for each histology from **Results** -> **Variation by histology**

> 📝 **Note:**
> The drawn raw data can be downloaded as an Excel file from **Downloadable table**.

</details>

<details>
<summary><b>I want a lolliplot (a diagram showing which amino acid residues of a gene have many mutations)</b></summary>

1. Perform filtering such as histology and age from **Setting**
4. Select a gene and check the results from **Results** -> **Oncoprint** -> **Figures** -> **Lolliplot for the selected gene**

> 📝 **Note:**
> Currently, it does not support exon skipping or intron mutations.

</details>

<details>
<summary><b>I want information on mutual exclusivity and co-occurring mutations between genes</b></summary>

1. Perform filtering such as histology and age from **Setting**
3. Check the results from **Results** -> **Mutually exclusivity**

> 💡 **How to read:**
> If the color of the cell where the gene on the X-axis and the gene on the Y-axis intersect is **blue**, the two are mutually exclusive, and if it is **red**, they have a co-occurring mutation relationship.

</details>

---

## 🎯 Treatment reachability analysis

<details>
<summary><b>I want to know what kind of patients have a high treatment reach rate when a CGP test is performed</b></summary>

1. Perform filtering such as histology, age, and treatment course from **Setting**. In particular, specify the gene set of interest in **Genes of interest**.
2. Obtain a nomogram that predicts the treatment reach rate based on clinical information of the patient obtained before the test from **Results** -> **CGP benefit prediction** -> **Factors leading to treatment** -> **Factors leading to treatment, pre-CGP, Nomogram**.
3. Analyze the effect of clinical information of the patient obtained before the test on the treatment reach rate in univariate and multivariate analysis from **Results** -> **CGP benefit prediction** -> **Factors leading to treatment** -> **Factors leading to treatment, pre-CGP, Odds ratio**.
4. Check the accuracy of the prediction by the nomogram with an ROC curve from **Results** -> **CGP benefit prediction** -> **Factors leading to treatment** -> **ROC curve of nomogram**.
5. Check the results of evaluating the clinical utility of the prediction by the nomogram with decision curve analysis from **Results** -> **CGP benefit prediction** -> **Factors leading to treatment** -> **Factors leading to treatment, decision curve**.
6. The treatment reach rate will be predicted when you enter the information of a specific patient from **Results** -> **CGP benefit prediction** -> **Factors leading to treatment** -> **Input your data**.

> 📖 **Reference:**
> For decision curve analysis, please refer <a href="https://mskcc-epi-bio.github.io/decisioncurveanalysis/index.html" target="_blank">here</a> or <a href="https://github.com/MANO-B/FELIS/blob/main/decision_curve_analysis.md" target="_blank">here</a>.

</details>

<details>
<summary><b>I want to know the relationship between survival time after CGP testing and treatment reach rate</b></summary>

1. Perform filtering such as histology, age, and treatment course from **Setting**. In particular, specify the gene set of interest in **Genes of interest**.
2. Check the graph taking the moving average of the survival time after CGP testing and the treatment reach rate from **Results** -> **Survival after CGP** -> **Survival analysis** -> **Survival period and treatment reach rate**.

> 💡 **Consideration:**
> Since the treatment reach rate is low for patients who die shortly after the CGP test, it seems important how early to perform the test on patients with a likely poor prognosis, and to perform the test on patients whose PS is maintained and have a certain degree of survival time. The recommended treatment reach rate after CGP when death is treated as a competing risk can also be considered as the proportion of patients who can reach the recommended treatment by Upfront CGP.

</details>

---

## ⏳ Left-truncation bias simulation

<details>
<summary><b>I want to check about left-truncation bias</b></summary>

As in the data of C-CAT, when the start date of measurement of survival time and the test date (start date of observation) are different, it is difficult to estimate the survival time using the normal Kaplan-Meier method. This is because all cases are alive from the start date of measurement of survival time to the test date, and immortal bias exists.
When the survival time is divided into the survival time from the start date of measurement to the test date, and the survival time from the test date to the final observation date, part of the bias is eliminated.
However, it is ultimately unknown "how the patients who received the CGP test differ from those who did not," and we consider that eliminating selection bias to a certain extent is impossible.

1. Open **Results** -> **Bias correction simulation** -> **An example of bias adjustment**.
2. Adjust the parameters to your liking
3. Press the **Left-truncation bias adjustment simulation** button
4. The true survival curve, the normal Kaplan-Meier survival curve, the survival curve divided before and after the CGP test, and the survival curve by the bias elimination method using Bayesian estimation will be drawn.

> 💡 **Interpretation:**
> You can probably see that a generally good bias correction has been achieved.

</details>

---

## 💾 Others

<details>
<summary><b>I want to save the figure</b></summary>

Right-click the figure and save it with a name with the extension `.png`.

</details>

<details>
<summary><b>I want to save the summary table</b></summary>

Select the text on the screen, copy it, paste it into Word or Excel, and save it.

</details>

<details>
<summary><b>I want to save the raw data table</b></summary>

Save it with a name as an Excel file or CSV file using the button on the top left.

</details>
