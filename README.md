# MutScape: analytical toolkit for mutational landscape in cancer genomics

![](https://github.com/anitalu724/MutScape/blob/main/mutscape/data/movie.gif?raw=true)

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
MutScape accepts both VCF and MAF files as input data. 
For multiple files will be implemented simultaneously, the user should enter a format-limited TSV file. The detailed format please refer to `examples/tsv/testData_vcf.tsv` and `examples/tsv/testData_maf.tsv`

Commands are listed below:
* `--file`, `-f` : The relative path of the input TSV file.
* `--vcf_filter`, `-vf` : Some parameters for VCF filtering.
    - GI: Genome Interval (ex: `GI [1,3]` , `GI [2:4]` or `GI "{1: [*,*], 2 : [1, 300000]}"`)
    - CI: Caller Information (ex: `CI "15,15,0,0,0,0.05,8,8"`)
    - PA: Keep or exclude non-PASS tag (ex: `PA 1`)
    - AV: Artifact variant filter: FFPE filter (ex: `AV 0.9`)
* `--combine`, `-c` : No parameter is required. VCFs combination will be implemented if this command is ordered.
* `--vcf2maf`, `-v2m` : The parameter is `max_filter_ac` which is an integer when transforming files to MAF. (Refer to [vcf2maf](https://github.com/mskcc/vcf2maf))
* `--output`, `-o` : The path for storing output files. This path must end with a folder.
* `--meta`, `-m` : The path for storing metafiles. This path must end with a folder.
* `--maf_filter`, `-mf` : Some parameters for MAF filtering.
    - GI: Genome Interval (ex: `GI [1,3]` , `GI [2:4]` or `GI "{1: [*,*], 2 : [1, 300000]}"`)
    - CI: Caller Information (ex: `CI "15,15,0,0,0,0.05,8,8"`)
    - TE: Tissue Expression (ex: `TE [breast,5]`)
    - PF: Population Frequency (ex: `PF 1`)
    - HY: Hypermutation or Sample Exclusion (ex: `HY 500`)

#### VCF
For VCFs as input data, `-f`, `-c`, `-v2m`, `-o` and `-m` are required while `-vf` and `-mf` are optional. 
Some simple test commands are displayed below.
    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -vf GI [1,3] \
    -c \
    -v2m 8 \
    -o examples/output \
    -m examples/meta 
    
    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -vf GI "{1: [*,*], 2 : [1, 300000]}" CI "15,15,0,0,0,0.05,8,8" PA 0 AV 0.9 \
    -c \
    -v2m 8 \
    -o examples/output \
    -m examples/meta
    
    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -c \
    -v2m 8 \
    -o examples/output \
    -m examples/meta \
    -mf GI [1,3]


#### MAF
For MAFs as input data, `-f`, `-o` and `-m` are required while `-mf` are optional. 
Some simple test commands are displayed below.

    python3 dataPreprocess.py \
    -f examples/tsv/testData_maf.tsv \
    -mf GI [1:3] \
    -o examples/output \
    -m examples/meta 
    
    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_maf.tsv \
    -mf GI [1:3] CI "15,15,0,0,0,0.05,8,8" TE [breast,5] PF 1 HY 500 \
    -o examples/output \
    -m examples/meta 
    
### MAF Analysis and Visualization
MutScape provides 8 different analyses and some of them generate plots after analysis.

1. Significantly mutated gene detection
	```
	python3 mafAnalysis.py \
	-f examples/test_data/maf/TCGA_test.maf \
	-smg \
	-o examples/output \
	-p examples/pic/
	```
2. Known cancer gene annotation
    ```
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -kcga \
    -o examples/output \
	-p examples/pic/
    ```
3. Mutation burden statistics
    ```
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -tmb 60456963 \
    -o examples/output \
	-p examples/pic/
    ```
4. CoMut plot analysis
    ```
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -cm \
    -o examples/output \
	-p examples/pic/


    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -cmp examples/tsv/comut.tsv examples/tsv/comut_info.tsv 0 comut.pdf \
    -o examples/output \
	-p examples/pic/
    ```
5. Mutational signature
    ```
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
    ```
6. HRD Score
    ```
    python3 mafAnalysis.py \
    -f examples/test_data/maf/hrd.maf \
    -hrd examples/tsv/hrd.tsv grch37 \
    -o examples/output \
	-p examples/pic/
    ```
7. Whole-genome doubling (WGD) and Chromosome instability (CIN)
    ```
    python3 mafAnalysis.py \
    -f examples/test_data/maf/hrd.maf \
    -wgdcin examples/tsv/hrd.tsv \
    -o examples/output \
	-p examples/pic/
    ```
8. Actionable mutation(drug) annotation
    ```
    ```
