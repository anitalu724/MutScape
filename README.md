# MutScape: analytical toolkit for mutational landscape in cancer genomics

## Introduction
MutScape provides a comprehensive pipeline of filtering, combination, transformation, analysis and visualization. MutScape can not only preprocess millions of mutation records in a few minutes, but offers various analyses, including driver gene detection, mutational signature, large-scale alteration identification, and actionable biomarker annotation simultaneously. Furthermore, MutScape supports both somatic variants in Variant Call Format (VCF) and Mutation Annotation Format (MAF), and leverages caller combination strategies to quickly eliminate false-positives. With only two simple commands, robust results and publication-quality images are generated automatically. Herein, we demonstrate the performance of MutScape using breast cancer samples from The Cancer Genome Atlas (TCGA) that correctly reproduce known results. More significantly, it enables discovery of novel results for cancer genomics studies through the advanced features in MutScape.
## Installation
### Prerequisite

This toolkit is recommended to implement in a new environment of conda. If you have never install conda, please refer to [Miniconda website](https://docs.conda.io/en/latest/miniconda.html).
After installed Miniconda, please create a brand-new environment and activate it:
```shell
conda create --name [env_name]
source activate [env_name]
```

### Install MutScape

MutScape is provided on Github website, please download it.

    git clone https://github.com/anitalu724/MutScape.git

### Install python modules
If you have never install `pip`, install it by `conda`.

    conda install -c anaconda pip

To make sure our code smoothly implement, we need to install several modules that are used in MutScape:

    bash MutScape/mutscape/installation/install_module.sh


### Install vcf2maf module
For the VCF transforms to the MAF, this procedure is implemented by `vcf2maf` utility, which processes variant annotation and transcript prioritization. You can refer to  [this script]((https://github.com/mskcc/vcf2maf)) or just follow the commands below. (Before this step, you must be sure that you have install [ Ensembl's VEP](https://gist.github.com/ckandoth/61c65ba96b011f286220fa4832ad2bc0))

    bash MutScape/mutscape/installation/install_vcf2maf.sh
    mv ../MutScape ./
## Implementation
MutScape has simply separated into two main modules: data preprocessing and analysis and visualization. 

### Data Preprocessing



