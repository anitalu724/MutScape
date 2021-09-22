# MutScape: an analytical toolkit for probing the mutational landscape in cancer genomics

![all_figure](https://github.com/anitalu724/MutScape/blob/main/mutscape/data/MutScape_movie.gif?raw=true)


## Introduction
Cancer genomics has been evolving rapidly, fueled by the emergence of numerous studies and public databases through next-generation sequencing technologies. However, the downstream programs used to preprocess and analyze data on somatic mutations are scattered in different tools, most of which need a long time for computation and require specific input formats. Here, we developed a user-friendly Python toolkit, MutScape, which provides a comprehensive pipeline of filtering, combination, transformation, analysis, and visualization. MutScape can not only preprocess millions of mutation records in a few minutes, but offers various analyses simultaneously, including driver gene detection, mutational signature, large-scale alteration identification, and actionable biomarker annotation. Furthermore, MutScape supports somatic variant data in both Variant Call Format (VCF) and Mutation Annotation Format (MAF), and leverages caller combination strategies to quickly eliminate false-positives. With only two simple commands, robust results and publication- quality images are generated automatically. Herein, we demonstrate the ability of MutScape to correctly reproduce known results using breast cancer samples from The Cancer Genome Atlas. More significantly, discovery of novel results in cancer genomics studies is enabled through the advanced features in MutScape.


## Quick Installation
Before implement quick installation, please be sure that you have installed MiniConda3, created a new conda environment and activate it. Also, to make this implementation run smoothly, please confirm that the Internet is connected always and the server/computer has enough storage memory.

    git clone https://github.com/anitalu724/MutScape.git
    bash MutScape/mutscape/installation/quickInstall_1.sh
    bash vcf2maf-1.6.20/MutScape/mutscape/installation/quickInstall_2.sh

## Prerequisite installation

### Requirements ####
The latest tested version in parentheses:
1. Using Miniconda (py37_4.9.2) to install:

    samtools (v1.10), ucsc-liftover (v377), bcftools (v1.10.2), htslib (v1.10.2) and ensembl-vep (v102.0)

2. Download vcf2maf (v1.6.20) and git clone MutScape (v1.0)

3. Download VEP cache data of GRCh37 and the reference FASTA (v102.0)


### Install Miniconda3
Numerous modules for this toolkit will be installed by `conda`.
If you have never installed conda, please refer to [Miniconda website](https://docs.conda.io/en/latest/miniconda.html). For high compatibility, we recommended users install `Miniconda3-py37_4.9.2`. (SHA256 hash  79510c6e7bd9e012856e25dcb21b3e093aa4ac8113d9aa7e82a86987eabe1c31)

There is a script for users to install Miniconda quickly.

    wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.9.2-Linux-x86_64.sh
    sha256sum Miniconda3-py37_4.9.2-Linux-x86_64.sh
    bash Miniconda3-py37_4.9.2-Linux-x86_64.sh
    export PATH="$HOME/miniconda3/bin:$PATH"

MutScape is preferred to be implementing under a brand-new conda environment.

    conda create --name MutScape
    conda activate MutScape

### Install Ensembl's VEP
If you have already install Ensembl's VEP, you may skip this part and directly into the next part to install `vcf2maf`. (However, you must confirm that your VEP version is compatible to vcf2maf. Here, we recommended installing `ensembl-vep=102.0`. )

    conda install -c bioconda -c conda-forge samtools=1.10 ucsc-liftover=377 bcftools=1.10.2 htslib==1.10.2
    conda install -c bioconda -c conda-forge -c defaults ensembl-vep=102.0 

### Install vcf2maf
For transforming the VCF into the MAF, this procedure is implemented by `vcf2maf` utility, which processes variant annotation and transcript prioritization. You can refer to  [this script]((https://github.com/mskcc/vcf2maf)) or just follow the commands below. (Before this step, you must be sure that you have installed [Ensembl's VEP](https://gist.github.com/ckandoth/61c65ba96b011f286220fa4832ad2bc0))

    wget https://github.com/mskcc/vcf2maf/archive/refs/tags/v1.6.20.tar.gz
    tar -zxf v1.6.20.tar.gz
    cd vcf2maf-1.6.20
    perl vcf2maf.pl --man
    perl maf2maf.pl --man

Before we start to use vcf2maf, we need to download VEP cache data and the reference FASTA.<br/>
***:warning:  Since these files are quite large, it may take a long time to download them!*** <br/>
***:warning:  Be sure that your available memory is at least 30GB!*** <br/>
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
MutScape is provided on the Github website, please download it.

    git clone https://github.com/anitalu724/MutScape.git

If you have never installed `pip`, install it by `conda`.

    conda install -c anaconda pip

To make sure all code smoothly implement, you need to install several modules that are used in MutScape:

    cd MutScape/mutscape
    bash installation/install_module.sh

## Implementation

MutScape has simply separated into two main modules: data preprocessing and analysis and visualization. Detailed structure please refer to [**Fig. 1**](https://github.com/anitalu724/MutScape/blob/main/mutscape/data/Fig1.pdf).

### Data Preprocessing
MutScape accepts both VCF and MAF files as input data. 
For multiple VCF/MAF files will be implemented simultaneously, MutScape requires a limited-format TSV file as input. For the detailed format please refer to example files such as `examples/tsv/testData_vcf.tsv` and `examples/tsv/testData_maf.tsv` or just see [Wiki](https://github.com/anitalu724/MutScape/wiki/Column-information-of-input-TSV-file).

#### Quick start from VCFs
For VCFs as input data, `-f`, `-o` and `-m` are required while `-vf`, `-ra`, `-v2m` and `-mf` are optional. 
Some simple test commands are displayed below.</br>
See [Wiki](https://github.com/anitalu724/MutScape/wiki/Specific-arguments-in-data-preprocessing) for detailed information.

![S1](https://github.com/anitalu724/MutScape/blob/main/mutscape/data/S1.gif?raw=true)

    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -o examples/output \
    -m examples/meta \
    -vf CI "*,*,*,6,*,*,*,*"


    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -o examples/output \
    -m examples/meta \
    -vf GI [1,3] \
    -v2m 8 
    
    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -o examples/output \
    -m examples/meta \
    -vf GI "{1: [*,*], 2 : [1, 300000]}" CI "15,15,0,6,0,0.05,8,8" PA 0 AV 0.9 \
    -v2m 
    
    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -o examples/output \
    -m examples/meta \
    -v2m 8 \
    -mf GI [1,3]


    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -ra examples/test_data/vcf/reject.vcf examples/test_data/vcf/accept.vcf \
    -o examples/output \
    -m examples/meta \
    -vf CI "*,*,*,6,*,*,*,*" \
    -v2m 8 \
    -mf GI [1,3]


#### Quick start from MAFs
For MAFs as input data, `-f`, `-o` and `-m` are required while `-mf` are optional. 
Some simple test commands are displayed below.

    python3 dataPreprocess.py \
    -f examples/tsv/testData_maf.tsv \
    -mf GI [1:3] \
    -o examples/output \
    -m examples/meta 
    
    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_maf.tsv \
    -mf GI [1:3] CI "15,15,0,0,0,0.05,8,8" TE [BLCA,5] PAC 1 HY 500 \
    -o examples/output \
    -m examples/meta


### Analysis and Visualization
MutScape provides 9 different analyses and some of them generate plots after analysis.</br>
See [Wiki](https://github.com/anitalu724/MutScape/wiki/Table-of-arguments-for-analysis-and-visualization) for detailed information.

#### Quick start
Some simple test commands are displayed below.

1. Significantly mutated gene detection
    
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -smg \
    -o examples/output \
    -p examples/pic/


2. Known cancer gene annotation
    
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -kcga \
    -o examples/output \
    -p examples/pic/
    

3. Mutation burden statistics
    
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -tmb 60456963 \
    -o examples/output \
    -p examples/pic/
    
4. CoMut plot analysis
    > Output figure is shown like [**Fig. 2**](https://github.com/anitalu724/MutScape/blob/main/mutscape/examples/images/Fig2.pdf).</br>
    > See [Wiki](https://github.com/anitalu724/MutScape/wiki/Column-information-of-input-TSV-file-for-CoMut-plot) for detailed information.

    
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -cm 60456963 \
    -o examples/output \
    -p examples/pic/


    python3 mafAnalysis.py \
    -cmp examples/tsv/comut.tsv examples/tsv/comut_info.tsv 0 comut.pdf \
    -o examples/output \
    -p examples/pic/
    
5. Mutational signature
   > Output figure is shown like [**Fig. 3**](https://github.com/anitalu724/MutScape/blob/main/mutscape/examples/images/Fig3.pdf).

    python3 mafAnalysis.py \
    -f examples/test_data/maf/ms.maf \
    -ms 0 "[SBS1, SBS5, SBS40, SBS87]" \
    -o examples/output \
    -p examples/pic/


    python3 mafAnalysis.py \
    -f examples/test_data/maf/ms.maf \
    -ms 1 "[2,9,10]" \
    -o examples/output \
    -p examples/pic/


    python3 mafAnalysis.py \
    -f examples/test_data/maf/ms.maf \
    -ms 2 "[3]" \
    -o examples/output \
    -p examples/pic/


6. HRD Score
   > Output figure is shown like [**Fig. 4A, B**](https://github.com/anitalu724/MutScape/blob/main/mutscape/examples/images/Fig4.pdf).
    
    python3 mafAnalysis.py \
    -hrd examples/tsv/hrd.tsv grch37 \
    -o examples/output \
    -p examples/pic/
    

7. Whole-genome doubling (WGD) and Chromosome instability (CIN)
   > Output figure is shown like [**Fig. 4C, D**](https://github.com/anitalu724/MutScape/blob/main/mutscape/examples/images/Fig4.pdf).
    
    python3 mafAnalysis.py \
    -wgdcin examples/tsv/hrd.tsv \
    -o examples/output \
    -p examples/pic/
    

8. HRD, CIN and WGD Comparison
    
    python3 mafAnalysis.py \
    -hcwc examples/tsv/hcw_comparison.tsv grch37 \
    -o examples/output \
    -p examples/pic/
    

9. Actionable mutation (drug) annotation
   `[your_oncokb_token]` is gotten from [OncoKB Website](https://www.oncokb.org/). You must create  your own account and get your personal API token.
   > Output figure is shown like [**Fig. 5**](https://github.com/anitalu724/MutScape/blob/main/mutscape/examples/images/Fig5.pdf).
    
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -oncokb ../oncokb-annotator/ [your_oncokb_token] 4 examples/test_data/oncokb/clinical_input.txt \
    -o examples/output \
    -p examples/pic/
    

## Reference
### If you use MutScape in your work, please cite
> Cheng-Hua Lu*, Chia-Shin Wu*, Mong-Hsun Tsai, Liang-Chuan Lai, Eric Y. Chuang (2021) MutScape: an analytical toolkit for probing the mutational landscape in cancer genomics (Submitted)



