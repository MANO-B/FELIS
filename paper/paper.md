---
title: "FELIS: An R package for reproducible analysis of real-world clinicogenomic sequencing data"
tags:
  - R
  - Shiny
  - clinical genomics
  - survival analysis
  - real-world data
  - cancer genomics
authors:
  - name: Masachika Ikegami
    orcid: 0000-0003-0491-7142   # ← ORCIDを入れてください
    affiliation: 1
affiliations:
  - name: National Cancer Center Japan
    index: 1
date: 2025-12-25
bibliography: paper.bib
---

## Summary

Large-scale real-world clinicogenomic datasets provide critical opportunities to study rare cancers, uncommon molecular alterations, and treatment outcomes in routine clinical practice. However, such datasets pose substantial analytical challenges, including data heterogeneity, selection bias, left truncation, and complex confounding structures.

FELIS (Flexible Exploration for LIquid and Solid tumor clinical sequencing data) is an open-source R package that provides an interactive, reproducible, and transparent framework for analyzing real-world clinical sequencing data from liquid and solid tumors. FELIS integrates exploratory visualization, survival analysis, and propensity score–based causal inference methods within a Shiny-based graphical user interface, enabling researchers to perform statistically rigorous analyses without sacrificing transparency or reproducibility.

The software was developed to support analyses of data from the Center for Cancer Genomics and Advanced Therapeutics (C-CAT), Japan’s national clinicogenomic database, but is applicable to other large-scale real-world sequencing datasets with similar structure.

## Statement of Need

The increasing adoption of comprehensive genomic profiling in routine oncology practice has generated large volumes of real-world sequencing data. While numerous tools exist for variant annotation and molecular interpretation, there is a lack of open-source software specifically designed to address downstream analytical challenges unique to real-world clinicogenomic data, such as delayed testing, left truncation, selection bias, and treatment heterogeneity.

FELIS addresses this gap by providing a unified analytical environment that combines bias-aware statistical modeling with interactive data exploration. By packaging these methods in a reproducible R package with a graphical user interface, FELIS lowers the barrier to rigorous analysis for clinical researchers and facilitates transparent, auditable research workflows suitable for regulatory and academic settings.

## Software Description

### Design and Implementation

FELIS is implemented as an R package and distributed via GitHub. The analytical core is written in R and leverages widely used statistical and visualization libraries. The user interface is implemented using Shiny, allowing interactive cohort selection, visualization, and model fitting.

The software can be executed locally or deployed on a Shiny Server. Application logic is bundled within the package, enabling version-controlled deployment and simplifying maintenance.

### Functionality

Key functionalities include:

- Interactive exploration of clinicogenomic cohorts  
- Time-to-event analyses using Kaplan–Meier estimation and Cox proportional hazards models  
- Multivariable survival modeling with covariate adjustment  
- Propensity score matching and weighting with diagnostic visualization  
- Visualization of cohort characteristics and treatment outcomes  
- Support for tumor-only and tumor–normal sequencing data  

### Reproducibility

FELIS follows standard R package development practices. All analytical steps are explicitly defined in code, outputs are written only to user-writable directories, and version history is tracked via GitHub. This design facilitates reproducible research and independent verification.

## Related Work

The methodological framework implemented in FELIS builds on prior peer-reviewed studies addressing biases in real-world clinicogenomic data, including delayed testing and left truncation. In particular, FELIS extends methods described by Tamura et al. (2023) for bias-aware analysis of comprehensive genomic profiling data.

## Availability

- **Source code:** https://github.com/MANO-B/FELIS  
- **License:** MIT  
- **Operating systems:** Platform independent  
- **Programming language:** R  

## References
