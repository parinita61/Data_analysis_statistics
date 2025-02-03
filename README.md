## CpG Methylation Analysis - Bioinformatics Coding Challenge

## Overview

This repository contains code for analyzing CpG methylation data to assess the reliability of phased methylation patterns (PMPs) as biomarkers for tissue differentiation. The tasks are divided into three sub-tasks:

### Sub-task 1: Coverage Analysis
- **Goal**: Calculate the median and coefficient of variation (CV) for single CpG coverage in each tissue type and generate visual plots summarizing coverage statistics.
- **Key Steps**: 
  - Compute the median and CV for coverage.
  - Visualize coverage statistics using plots.

### Sub-task 2: Biomarker Identification
- **Goal**: Identify PMPs with high specificity for tissue differentiation, minimize false positives for Tissue #1, and assign confidence to each PMP. Also, calculate the mean variant read fraction (VRF) for each PMP.
- **Key Steps**:
  - Statistical or machine learning approach to assign confidence (e.g., p-values).
  - Calculate mean VRF for each PMP in both tissues.

### Sub-task 3: Addressing Key Questions
- **Goal**: Analyze the impact of sequencing depth on specificity, estimate the required reads for confident classification, and validate the hypothesis by comparing the specificity of top PMPs with individual CpG sites.
- **Key Steps**:
  - Assess the effect of sequencing depth on specificity.
  - Estimate the read threshold for confident classification of Tissue #2 at a sequencing depth of 1 million reads.
  - Compare the specificity of the top PMPs with individual CpG sites.

## Usage

- Clone this repository and run the respective Python scripts for each sub-task.
- Make sure the dataset is placed in the correct directory as mentioned in the scripts.

