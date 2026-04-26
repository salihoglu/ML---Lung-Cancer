# Ferroptosis-Linked Prognostic Modeling in Lung Adenocarcinoma

> **Status:** Manuscript submitted to a peer-reviewed journal; not yet published.

## Overview

This repository accompanies an unpublished scientific study entitled:

**Multi-Cohort Machine Learning Identifies an Externally Validated Ferroptosis-Linked Prognostic Signature in Lung Adenocarcinoma**

The study presents a multi-cohort machine-learning framework for identifying a reproducible ferroptosis-linked prognostic signature in lung adenocarcinoma (LUAD). The analysis integrates ferroptosis annotation resources, co-expression network biology, survival machine learning, external validation, tumour microenvironment profiling, model interpretability, and structure-informed druggability assessment.

LUAD remains clinically heterogeneous, and conventional clinicopathological staging does not fully explain patient-level survival variability. This work investigates whether ferroptosis-related genes can be systematically prioritized into a stable prognostic model that generalizes across independent transcriptomic cohorts.

## Study Summary

The candidate gene space was defined by intersecting **FerrDb V3** ferroptosis annotations with a ferroptosis-associated **WGCNA** module derived from the independent **GSE81089** cohort.

A final **25-gene prognostic signature** was locked using:

- bagged meta-univariate survival screening
- leave-one-cohort-out stability analysis
- multi-cohort model benchmarking

Training was performed across four LUAD cohorts:

- TCGA-LUAD
- GSE31210
- GSE72094
- GSE136961

Together, these training cohorts included **1,149 patients** and **344 survival events**.

A total of **150 configurations** from **nine survival-learning model families** were benchmarked. Final model performance was assessed in three strictly held-out external validation cohorts:

- GSE50081
- GSE68465
- GSE30219-ADC

These external cohorts included **912 patients** and **507 survival events**.

## Key Findings

The best internal model achieved a leave-one-cohort-out **Uno C-index of 0.723**.

External validation demonstrated consistent prognostic performance, with C-indices of:

- **0.611**
- **0.691**
- **0.699**

Across held-out cohorts, the continuous-risk score was associated with overall survival, with a random-effects pooled hazard ratio of:

**HR = 1.42, 95% CI 1.15–1.75, p = 1.0 × 10⁻³, I² = 25%**

In TCGA-LUAD, the signature remained independently prognostic after adjustment for age, sex, and pathological stage:

**HR = 1.84, 95% CI 1.38–2.47**

Biological interpretation supported an immune-inflamed, interferon- and complement-enriched low-risk phenotype. Additional analyses nominated structure-supported candidate targets for future experimental follow-up.

## Analytical Components

This repository may include code and analysis modules for:

- ferroptosis-related gene set curation
- WGCNA-based module prioritization
- multi-cohort survival modeling
- leave-one-cohort-out validation
- external validation
- risk-score generation
- pathway enrichment analysis
- tumour microenvironment deconvolution
- STRING network analysis
- permutation and SHAP-based interpretability
- immunotherapy response-propensity analysis
- signature-reversal connectivity analysis
- structure-informed druggability prioritization

## Manuscript Information

**Title:**  
*Multi-Cohort Machine Learning Identifies an Externally Validated Ferroptosis-Linked Prognostic Signature in Lung Adenocarcinoma*

**Running title:**  
Ferroptosis-linked immune signature and targets in LUAD

**Author:**  
Rana Salihoglu

**Affiliations:**  
Department of Computer Engineering, Recep Tayyip Erdogan University, Rize, Türkiye  
Chair of Bioinformatics, Julius-Maximilians-Universität Würzburg, Würzburg, Germany

**Keywords:**  
ferroptosis; lung adenocarcinoma; machine learning; drug target identification; WGCNA; tumour microenvironment; external validation; structural biology; druggability

## Citation

The associated manuscript has been submitted to a peer-reviewed journal.

Citation details will be added after publication.

## Contact

Rana Salihoglu  
ranasalihoglu@gmail.com
