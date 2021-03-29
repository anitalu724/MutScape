# MutScape: analytical toolkit for mutational landscape in cancer genomics

## Introduction
MutScape provides a comprehensive pipeline of filtering, combination, transformation, analysis and visualization. MutScape can not only preprocess millions of mutation records in a few minutes, but offers various analyses, including driver gene detection, mutational signature, large-scale alteration identification, and actionable biomarker annotation simultaneously. Furthermore, MutScape supports both somatic variants in Variant Call Format (VCF) and Mutation Annotation Format (MAF), and leverages caller combination strategies to quickly eliminate false-positives. With only two simple commands, robust results and publication-quality images are generated automatically. Herein, we demonstrate the performance of MutScape using breast cancer samples from The Cancer Genome Atlas (TCGA) that correctly reproduce known results. More significantly, it enables discovery of novel results for cancer genomics studies through the advanced features in MutScape.

## Implementation
MutScape has simply separated into two main modules: data preprocessing and analysis and visualization. 
![](https://github.com/anitalu724/MutScape/blob/main/mutscape/examples/images/fig1.pdf)
### Data Preprocessing

