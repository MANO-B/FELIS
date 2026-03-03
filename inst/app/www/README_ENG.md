# FELIS: Flexible Exploration for LIquid and Solid tumor clinical sequencing data
### — C-CAT Secondary Use Data Analysis Platform —

FELIS is a local-execution Web application (R/Shiny) designed for researchers to perform cohort definition, visualization, and bias-aware outcome analysis using C-CAT secondary use data.

> **Note**: This software is intended only for those who can legally obtain C-CAT secondary use data. Please comply with ethical reviews and data usage terms.

---

## 1. Overview
C-CAT (Center for Cancer Genomics and Data Management) is the national database for cancer genomic medicine in Japan. 
FELIS aims to enable "no-code/low-code" iterative exploration of:
1. **Cohort construction**
2. **Mutation summarization and visualization**
3. **Treatment and prognosis analysis**

Specifically, it integrates analysis workflows that account for **delayed entry** and **left truncation**, which are common issues in real-world CGP (Cancer Genomic Profiling) data.

---

## 2. GUI Structure (Menu Mapping)

- **Settings**: Specify analysis targets (cases, histology, genes, treatments), thresholds, and model settings.
- **Results**: Display and download charts/tables for each analysis.
- **Instruction / Tips**: Display this README and operational tips within the app.

---

## 3. Results: Detailed Reference of Outputs

Below is a comprehensive list of what each output in the FELIS menu signifies.

### ■ Case summary
- **Summarized by mutation pattern**: A table summarizing patient attributes and clinical variables based on specified mutation patterns.
- **Summarized by histology**: Case counts, clinical info, and mutation summaries grouped by histology (e.g., OncoTree).

### ■ Oncoprint
- **Oncoprint**: Visualization of the mutational landscape (Patients x Genes) for high-frequency genes in the selected cohort.
- **Lolliplot for the selected gene**: Frequency distribution of amino acid changes (hotspots) for a specific gene.
- **Table of clinical and mutation information per patient**: Raw patient-level data for verification and secondary analysis.

### ■ Mutual exclusivity
- **Figure for probability / odds ratio**: Visualization of co-occurrence or mutual exclusivity trends between gene pairs.
- **table_mutually_exclusive**: Statistical table (p-values, etc.) for each pair. 
  - *Interpretation*: Red typically indicates co-occurrence, Blue indicates exclusivity.

### ■ Variant rate by histology
- **figure_mut_subtype_plot**: Comparison of mutation rates of top genes across histological subtypes to identify characteristic mutations.

### ■ Mutation and treatment option
- **Basic data**: Statistics (Age, Sex, TMB, Treatment options, etc.) by histology.
- **UMAP clustering based on mutations**: Dimensionality reduction and clustering (DBSCAN/HDBSCAN) to identify patient clusters based on mutation patterns.
- **Cluster and histology relationship**: Enrichment analysis of clusters across histology and mutations.
- **Heterogeneity within histologic types**: Evaluation of cluster distribution (entropy) within each histological type.
- **Frequency of patients with targeted therapy**: Aggregated frequency of treatment options by evidence level.

### ■ Survival after CGP
- **Survival and clinical information**: Survival curves (KM, etc.) starting from the date of CGP testing.
- **Custom survival analysis**: GUI-defined two-group comparison with KM/RMST and adjustments like PSM/IPW (including Love plots and weight distributions).
- **Survival and mutations, forest plot**: Forest plot showing RMST differences or estimated effects stratified by mutation status.
- **Hazard ratio**: Multivariate analysis results (Cox model) in table format.
- **Survival period and treatment reach rate**: Cumulative Incidence Function (CIF) for "treatment reach," treating death as a competing risk.

### ■ CGP benefit prediction
- **Nomogram**: Predictive model (logistic regression) for reaching recommended treatment based on pre-CGP info.
- **Odds ratio**: ORs for factors related to treatment reach.
- **Decision curve**: Decision Curve Analysis (DCA) to evaluate the clinical net benefit of the model.
- **ROC curve of nomogram**: Evaluation of predictive accuracy (AUC).
- **Input data**: Simulator to return predicted values based on manually entered clinical variables.

### ■ Overall survival with risk-set adjustment (Left-truncation correction)
- **Survival and clinical information (after CTx)**: Survival estimation starting from "Chemotherapy initiation" with risk-set adjustment to handle the bias of delayed CGP entry.
- **Custom survival analysis**: Two-group comparison using the risk-set adjustment framework.
- **Frequent variants / Diagnosis / Mutational cluster and survival**: Survival analysis (Forest/Curves) with mutations, histology, or clusters as explanatory variables.
- **Hazard ratio for survival after CTx**: Cox model estimation results for corrected workflows.

### ■ Survival after CTx with Bayesian inference
- **Survival corrected for left-truncation bias**: Bayesian estimation using Stan to output survival curves with Credible Intervals (CI).
- **Custom survival analysis (BayesCustom)**: Two-group comparison within the Bayesian framework.
- **Genetic variants / Diagnosis and survival**: Corrected survival curves stratified by variants or histology.

### ■ Survival after CTx with control cohort data (Experimental)
- **Custom survival analysis (ControlCustom)**: Experimental analysis framework utilizing external control cohort information.

### ■ Bias correction simulation
- **Left-truncation bias adjustment simulation**: Simulation output to visualize how correction methods behave under different parameters. Used for sensitivity checks.

### ■ Drug response
- **Settings (Drugusebylineoftreatment)**: Define lines of therapy, ToT/TTF definitions, and drug groupings. Provides an overview of commonly used regimens.
- **Drug usage data (Drugperpatient)**: Raw data table of patient-drug interactions.
- **Treatment time and clinical information**: ToT KM curves, correlation with pre-treatment (scatter plots), and forest plots.
- **Treatment time comparison (ToT_interactive)**: Interactive GUI for comparing ToT between two defined groups.
- **ToT by gene mutation cluster / by mutated genes / by mutation pattern**: Stratification of ToT by clusters, specific mutations, or patterns.
- **Hazard ratio on time on treatment**: Cox estimation for ToT outcomes.
- **Volcano plot (ToT / ORR / AE)**: Exploration of associations between regimens/genes and outcomes (effect size vs. significance).
- **Cumulative incidence of adverse effect**: Visualization of AE occurrences considering competing risks.
- **Survival and drug (Survival_drug)**: Survival analysis starting from the "Drug initiation date."

---

## 4. Methodological Background (Bias Correction)
As detailed in the FELIS paper:
- C-CAT data requires secure/offline processing, making cloud-based platforms difficult to use.
- Real-world CGP data suffers from **selection bias** due to delayed entry/left truncation.
- FELIS addresses this by separating **Post-CGP analysis** (Naive KM/CIF) and **Post-CTx analysis** (Risk-set adjustment/Bayesian simulation).

---

## 5. References
If you use FELIS in your research, please cite:

```bibtex
@article{Mano2024FELIS,
  title={FELIS: Flexible Exploration for LIquid and Solid tumor clinical sequencing data},
  author={Mano, B.},
  journal={GitHub Repository},
  url={https://github.com/MANO-B/FELIS},
  year={2024}
}
```

---
**License**: MIT License
**Developer**: [MANO-B](https://github.com/MANO-B)
