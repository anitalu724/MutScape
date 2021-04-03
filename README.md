# MutScape: analytical toolkit for mutational landscape in cancer genomics

<blockquote class="imgur-embed-pub" lang="en" data-id="a/oC2EmEr" data-context="false" ><a href="//imgur.com/a/oC2EmEr"></a></blockquote><script async src="//s.imgur.com/min/embed.js" charset="utf-8"></script>

## Introduction
MutScape provides a comprehensive pipeline of filtering, combination, transformation, analysis and visualization. MutScape can not only preprocess millions of mutation records in a few minutes, but offers various analyses, including driver gene detection, mutational signature, large-scale alteration identification, and actionable biomarker annotation simultaneously. Furthermore, MutScape supports both somatic variants in Variant Call Format (VCF) and Mutation Annotation Format (MAF), and leverages caller combination strategies to quickly eliminate false-positives. With only two simple commands, robust results and publication-quality images are generated automatically. Herein, we demonstrate the performance of MutScape using breast cancer samples from The Cancer Genome Atlas (TCGA) that correctly reproduce known results. More significantly, it enables discovery of novel results for cancer genomics studies through the advanced features in MutScape.

## Prerequisite installation

### Install Miniconda3
Many modules for this toolkit will be installed by `conda`.
If you have never install conda, please refer to [Miniconda website](https://docs.conda.io/en/latest/miniconda.html). For high compatibility, we recommended users to install `Miniconda3-py37_4.9.2`. (SHA256 hash  79510c6e7bd9e012856e25dcb21b3e093aa4ac8113d9aa7e82a86987eabe1c31)

There is a script for users to install Miniconda quickly.

    wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.9.2-Linux-x86_64.sh
    sha256sum Miniconda3-py37_4.9.2-Linux-x86_64.sh
    bash Miniconda3-py37_4.9.2-Linux-x86_64.sh
    export PATH="$HOME/miniconda3/bin:$PATH"

MutScape is prefered to be implementing under a brand-new conda environment.

    conda create --name [env_name]
    conda activate [env_name]
### Install Ensembl's VEP
If you have already install Ensembl's VEP, you may skip this part and directly into the next part for install `vcf2maf`. (However, you must confirm that your VEP version is compatible to vcf2maf. Here, we recommended to install `ensembl-vep=102.0`. )

    conda install -c bioconda -c conda-forge samtools=1.10 ucsc-liftover=377 bcftools=1.10.2
    conda install -c bioconda -c conda-forge -c defaults ensembl-vep=102.0 

### Install vcf2maf
For the VCF transforms to the MAF, this procedure is implemented by `vcf2maf` utility, which processes variant annotation and transcript prioritization. You can refer to  [this script]((https://github.com/mskcc/vcf2maf)) or just follow the commands below. (Before this step, you must be sure that you have install [ Ensembl's VEP](https://gist.github.com/ckandoth/61c65ba96b011f286220fa4832ad2bc0))

    export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
    curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-*
    perl vcf2maf.pl --man
    perl maf2maf.pl --man

Before we start to use vcf2maf, we need some references.<br/>
***:warning:  Since these files are quite large, it may take a long time to download them!*** <br/>
:information_source:  Here we recommended to download **102_GRCh37**

    mkdir -p $HOME/.vep/homo_sapiens/102_GRCh37/
    wget ftp://ftp.ensembl.org/pub/grch37/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
    mv Homo_sapiens.GRCh37.dna.toplevel.fa.gz  $HOME/.vep/homo_sapiens/102_GRCh37/
    gzip -d $HOME/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
    bgzip -i $HOME/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa
    samtools faidx $HOME/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-102/variation/indexed_vep_cache/homo_sapiens_vep_102_GRCh37.tar.gz
    mv homo_sapiens_vep_102_GRCh37.tar.gz $HOME/.vep/
    tar -zxf $HOME/.vep/homo_sapiens_vep_102_GRCh37.tar.gz -C $HOME/.vep/
    
### Install MutScape
MutScape is provided on Github website, please download it.

    git clone https://github.com/anitalu724/MutScape.git

## Implementation

If you have never install `pip`, install it by `conda`.

    conda install -c anaconda pip

To make sure our code smoothly implement, we need to install several modules that are used in MutScape:

    cd MutScape/mutscape
    bash installation/install_module.sh

MutScape has simply separated into two main modules: data preprocessing and analysis and visualization. 

### Data Preprocessing

    python3 dataPreprocess.py \
    -f examples/tsv/testData.tsv \
    -vf GI [1,3] \
    -c \
    -v2m 8 \
    -o examples/output \
    -m examples/meta 

