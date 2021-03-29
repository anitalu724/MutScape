# MutScape: analytical toolkit for mutational landscape in cancer genomics

## Introduction
MutScape provides a comprehensive pipeline of filtering, combination, transformation, analysis and visualization. MutScape can not only preprocess millions of mutation records in a few minutes, but offers various analyses, including driver gene detection, mutational signature, large-scale alteration identification, and actionable biomarker annotation simultaneously. Furthermore, MutScape supports both somatic variants in Variant Call Format (VCF) and Mutation Annotation Format (MAF), and leverages caller combination strategies to quickly eliminate false-positives. With only two simple commands, robust results and publication-quality images are generated automatically. Herein, we demonstrate the performance of MutScape using breast cancer samples from The Cancer Genome Atlas (TCGA) that correctly reproduce known results. More significantly, it enables discovery of novel results for cancer genomics studies through the advanced features in MutScape.
## Installation
* **Prerequisite**
This toolkit is recommended to implement in a new environment of conda. If you have never install conda, please refer to [Miniconda website](https://docs.conda.io/en/latest/miniconda.html).
After installed Miniconda, please create a brand-new environment and activate it:
```shell=
conda create --name [env_name]
source activate [env_name]
```

* **Install useful modules**
To make sure our code smoothly implement, we need to install several modules that are used in MutScape:
```shell=
pip install termcolor
pip install tqdm
pip install numpy
pip install pandas
pip install vcfpy
pip install seaborn
pip install oncodriveclust
pip3 install comut
pip install sklearn
pip install fonttools
pip install PyVCF
```


## Implementation
MutScape has simply separated into two main modules: data preprocessing and analysis and visualization. 

### Data Preprocessing



