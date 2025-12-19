# Reviewer FAQ

## What problem does FELIS solve?
FELIS enables reproducible analysis of large-scale real-world clinicogenomic
data.

## How does FELIS differ from existing tools?
It is specifically designed for C-CAT data and integrates bias-aware statistical
methods.

## Is FELIS reproducible?
Yes. FELIS follows standard R package conventions and version control practices.

## Is there related publications?
Yes. Some analytical methods implemented in FELIS have been previously
published.

- Mochizuki T et al., Cancer Science, 2024: Variant-based clustering analysis
  implemented in the genomic stratification module.
- Tamura T et al., Cancer Science, 2023: Survival analysis accounting for
  left truncation and delayed entry implemented in the outcome analysis module.

These publications describe the underlying methods, while the present JOSS
submission focuses on their software implementation, integration, and
practical usability in clinicogenomic research workflows.

## Why does FELIS require rstan / cmdstanr and such a large number of dependencies?
FELIS targets real-world clinicogenomic analyses where survival modeling must
explicitly account for delayed entry and left truncation, which are common in
comprehensive genomic profiling (CGP) practice. Several core analyses
implemented in FELIS rely on Bayesian or likelihood-based survival models that
are naturally and robustly expressed using the Stan probabilistic programming
framework.

As a result, rstan and cmdstanr are required dependencies rather than optional
components. The large dependency footprint reflects the breadth of clinically
oriented analyses provided by FELIS—including bias-aware survival analysis,
genomic stratification, clustering, visualization, and reporting—within a single
integrated application intended for end users without programming expertise.

