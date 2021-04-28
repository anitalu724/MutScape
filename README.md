# MutScape: analytical toolkit for mutational landscape in cancer genomics

![](https://github.com/anitalu724/MutScape/blob/main/mutscape/data/movie.gif?raw=true)

## Introduction
MutScape provides a comprehensive pipeline of filtering, combination, transformation, analysis and visualization. MutScape can not only preprocess millions of mutation records in a few minutes, but offers various analyses, including driver gene detection, mutational signature, large-scale alteration identification, and actionable biomarker annotation simultaneously. Furthermore, MutScape supports both somatic variants in Variant Call Format (VCF) and Mutation Annotation Format (MAF), and leverages caller combination strategies to quickly eliminate false-positives. With only two simple commands, robust results and publication-quality images are generated automatically. Herein, we demonstrate the performance of MutScape using breast cancer samples from The Cancer Genome Atlas (TCGA) that correctly reproduce known results. More significantly, it enables discovery of novel results for cancer genomics studies through the advanced features in MutScape.

## Prerequisite installation

### Requirements ####
Latest tested version in parentheses:
1. Using Miniconda (py37_4.9.2) to install:

    samtools (v1.10), ucsc-liftover (v377), bcftools (v1.10.2), ensembl-vep (v102.0)

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

    conda install -c bioconda -c conda-forge samtools=1.10 ucsc-liftover=377 bcftools=1.10.2
    conda install -c bioconda -c conda-forge -c defaults ensembl-vep=102.0 

### Install vcf2maf
For transforming the VCF into the MAF, this procedure is implemented by `vcf2maf` utility, which processes variant annotation and transcript prioritization. You can refer to  [this script]((https://github.com/mskcc/vcf2maf)) or just follow the commands below. (Before this step, you must be sure that you have installed [ Ensembl's VEP](https://gist.github.com/ckandoth/61c65ba96b011f286220fa4832ad2bc0))

    export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
    curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-*
    perl vcf2maf.pl --man
    perl maf2maf.pl --man

Before we start to use vcf2maf, we need to download VEP cache data and the reference FASTA.<br/>
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
MutScape is provided on the Github website, please download it.

    git clone https://github.com/anitalu724/MutScape.git

If you have never installed `pip`, install it by `conda`.

    conda install -c anaconda pip

To make sure all code smoothly implement, you need to install several modules that are used in MutScape:

    cd MutScape/mutscape
    bash installation/install_module.sh

## Implementation

MutScape has simply separated into two main modules: data preprocessing and analysis and visualization. Detailed structure please refer to **Fig1**.

### Data Preprocessing
MutScape accepts both VCF and MAF files as input data. 
For multiple VCF/MAF files will be implemented simultaneously, MutScape requires a limited-format TSV file as input. In the detailed format please refer to `examples/tsv/testData_vcf.tsv` and `examples/tsv/testData_maf.tsv`

#### Quick start from VCFs
For VCFs as input data, `-f`, `-c`, `-v2m`, `-o` and `-m` are required while `-vf` and `-mf` are optional. 
Some simple test commands are displayed below.

![S1](https://github.com/anitalu724/MutScape/blob/main/mutscape/data/S1.gif?raw=true)

    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -vf GI [1,3] \
    -v2m 8 \
    -o examples/output \
    -m examples/meta 
    
    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -vf GI "{1: [*,*], 2 : [1, 300000]}" CI "15,15,0,0,0,0.05,8,8" PA 0 AV 0.9 \
    -v2m 8 \
    -o examples/output \
    -m examples/meta
    
    
    python3 dataPreprocess.py \
    -f examples/tsv/testData_vcf.tsv \
    -v2m 8 \
    -o examples/output \
    -m examples/meta \
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
    -mf GI [1:3] CI "15,15,0,0,0,0.05,8,8" TE [breast,5] PAC 1 HY 500 \
    -o examples/output \
    -m examples/meta
    
#### MutScape specific arguments for data preprocessing
<table>
   <tr>
      <td>Argument name(s)</td>
      <td>Default value</td>
      <td>Descriptions</td>
      <td>Example</td>
   </tr>
   <tr>
      <td>Required arguments</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>--file, -f</td>
      <td>NULL</td>
      <td>The relative path of the input TSV file</td>
      <td>ex: -f examples/tsv/testData_vcf.tsv</td>
   </tr>
   <tr>
      <td>--output, -o</td>
      <td>NULL</td>
      <td>The path for storing output files. This path must end with a folder.</td>
      <td>ex: -o examples/output</td>
   </tr>
   <tr>
      <td>--meta, -m</td>
      <td>NULL</td>
      <td>The path for storing metafiles. This path must end with a folder.</td>
      <td>ex: -m examples/meta</td>
   </tr>
   <tr>
      <td>--combine, -c</td>
      <td>TRUE</td>
      <td>VCFs combination will be implemented if this command is ordered</td>
      <td>ex: -c</td>
   </tr>
   <tr>
      <td>--vcf2maf, -v2m</td>
      <td>10</td>
      <td>Annotated mutations as common in "FILTER" column of MAFs, which population allele counts across at least one ExAC subpopulation are > defult or a given threshold</td>
      <td>ex: -v2m 8</td>
   </tr>
   <tr>
      <td>Optional arguments</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>--vcf_filter, -vf</td>
      <td></td>
      <td>VCF filtering steps will only be implemented when one or more of the following parameters are given</td>
      <td></td>
   </tr>
   <tr>
      <td>GI</td>
      <td>[]</td>
      <td>Select mutations from one or more genomic positions or specific chromosomal regions</td>
      <td>ex: -vf  GI [1,3]  : chromosome 1 and 3; GI [2:4]  : chromosome 2-4; GI "{1: [*,*], 2 : [1, 300000]}"  : chromosome 1 and position 1-300000 in chromosome 2</td>
   </tr>
   <tr>
      <td>CI</td>
      <td>[]</td>
      <td>Sift valid mutations based on the information in the VCF file</td>
      <td>ex: -vf  CI "15,15,0,6,0,0.05,8,8"  : Sift mutations with the criteria (total depth of reads in normal>15, total depth of reads in tumor>15, mutant allelic depth in normal=0, mutant allelic depth in tumor≧6,  mutant allelic fraction in normal=0, mutant allelic fraction in tumor≧0.05, normal LOD score≧8, tumor LOD score≧8)</td>
   </tr>
   <tr>
      <td>PA</td>
      <td>0</td>
      <td>Keep or exclude mutations with non-PASS tag</td>
      <td>ex: -vf  PA 1  : exclude; -vf  PA 0  : keep </td>
   </tr>
   <tr>
      <td>AV</td>
      <td>[]</td>
      <td>Discard false-positive calls from FFPE samples with values > a given threshold, which is only for Mutect2 VCF files</td>
      <td>ex: -vf  AV 0.9  : Discard mutations with |(F1R2alt-F2R1alt)/(F1R2alt+F2R1alt)| > 0.9</td>
   </tr>
   <tr>
      <td>--maf_filter, -mf</td>
      <td></td>
      <td>MAF filtering steps will only be implemented when one or more of the following parameters are given</td>
      <td></td>
   </tr>
   <tr>
      <td>CI</td>
      <td>[]</td>
      <td>Sift valid muations based on the information in the MAF file</td>
      <td>ex: -mf  CI "15,15,0,6,"  ; Sift mutations with the criteria (total depth of reads in normal>15, total depth of reads in tumor>15, mutant allelic depth in normal=0, mutant allelic depth in tumor≧6) </td>
   </tr>
   <tr>
      <td>TE</td>
      <td>[]</td>
      <td>Exclude non- or low-expressed gene mutations in the specific tissue</td>
      <td>ex: -mf  TE [breast,5]  : exclude low-expressed gene mutations with average expression < 5 in the breast tissue from TCGA samples</td>
   </tr>
   <tr>
      <td>PF</td>
      <td>0</td>
      <td>Filter ‘common_variant’, which annotated by vcf2maf (-v2m)</td>
      <td>ex: -mf  PAC 1  : exclude; -mf  PF 0  : keep</td>
   </tr>
   <tr>
      <td>HY</td>
      <td>[]</td>
      <td>Exclude hypermutators to avoids statistical bias</td>
      <td>ex: -mf  HY 500  ; exclude samples with >500 mutations</td>
   </tr>
</table>

#### Column information of input TSV file (for VCFs):
<table>
   <tr>
      <td>Columns</td>
      <td>Descriptions</td>
   </tr>
   <tr>
      <td>NORMAL </td>
      <td>Each matched normal sample name </td>
   </tr>
   <tr>
      <td>TUMOR</td>
      <td>Each tumor sample name </td>
   </tr>
   <tr>
      <td>MuSE </td>
      <td>Paths of MuSE VCF for each sample </td>
   </tr>
   <tr>
      <td>Mutect2</td>
      <td>Paths of Mutect2 VCF for each sample </td>
   </tr>
   <tr>
      <td>SomaticSniper</td>
      <td>Paths of SomaticSniper VCF for each sample </td>
   </tr>
   <tr>
      <td>Strelka2</td>
      <td>Paths of Strelka2 VCF for each sample </td>
   </tr>
   <tr>
      <td>VarScan2</td>
      <td>Paths of VarScan2 VCF for each sample </td>
   </tr>
   <tr>
      <td>At Least # CALLS</td>
      <td>Self-defined criteria for further filtering (ex: Given that users enter 2 here, MutScape retains mutations that are identified by at least 2 variant callers) </td>
   </tr>
   <tr>
      <td>At Most # REJECT</td>
      <td>Self-defined criteria for further filtering (ex: Given that users enter 1 here, MutScape discards mutations that are rejected by at most 1 caller) </td>
   </tr>
</table>

#### Column information of input TSV file (for MAFs):
<table>
   <tr>
      <td>Columns</td>
      <td>Descriptions</td>
   </tr>
   <tr>
      <td>MAF</td>
      <td>Paths of MAF for each sample</td>
   </tr>
</table>

    
### MAF Analysis and Visualization
MutScape provides 8 different analyses and some of them generate plots after analysis.

#### Quick start
Some simple test commands are displayed below.

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
    > Output figure is shown like **Fig2**.

    ```
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -cm 60456963 \
    -o examples/output \
    -p examples/pic/


    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -cmp examples/tsv/comut.tsv examples/tsv/comut_info.tsv 0 comut.pdf \
    -o examples/output \
    -p examples/pic/
    ```
5. Mutational signature
   > Output figure is shown like **Fig3**.
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
   > Output figure is shown like **Fig4.A, B**.
    ```
    python3 mafAnalysis.py \
    -f examples/test_data/maf/hrd.maf \
    -hrd examples/tsv/hrd.tsv grch37 \
    -o examples/output \
    -p examples/pic/
    ```
7. Whole-genome doubling (WGD) and Chromosome instability (CIN)
   > Output figure is shown like **Fig4.C, D**.
    ```
    python3 mafAnalysis.py \
    -f examples/test_data/maf/hrd.maf \
    -wgdcin examples/tsv/hrd.tsv \
    -o examples/output \
    -p examples/pic/
    ```
8. Actionable mutation (drug) annotation
   `[your_oncokb_token]` is gotten from [OncoKB Website](https://www.oncokb.org/). You must create  your own account and get your personal API token.
   > Output figure is shown like **Fig5**.
    ```
    python3 mafAnalysis.py \
    -f examples/test_data/maf/TCGA_test.maf \
    -oncokb ../oncokb-annotator/ [your_oncokb_token] 4 examples/test_data/oncokb/clinical_input.txt \
    -o examples/output \
    -p examples/pic/
    ```

#### MutScape specific arguments for analysis and visualization
<table>
   <tr>
      <td>Argument name(s)</td>
      <td>Default value</td>
      <td>Descriptions</td>
      <td>Example</td>
   </tr>
   <tr>
      <td>Required arguments</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>--file, -f</td>
      <td>NULL</td>
      <td>The relative path of the input TSV file</td>
      <td>ex: -f examples/test_data/maf/TCGA_test.maf</td>
   </tr>
   <tr>
      <td>--output, -o</td>
      <td>NULL</td>
      <td>The path for storing output files. This path must end with a folder.</td>
      <td>ex: -o examples/output</td>
   </tr>
   <tr>
      <td>--picture, -p</td>
      <td>NULL</td>
      <td>The path for storing output images. This path must end with a folder.</td>
      <td>ex: -p examples/pic</td>
   </tr>
   <tr>
      <td>Optional arguments</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>--significantly_mutated_gene, -smg</td>
      <td>False </td>
      <td>Implement significantly mutated gene detection for input cohort</td>
      <td>ex: -smg</td>
   </tr>
   <tr>
      <td>--known_cancer_gene_annotaiton, -kcga</td>
      <td>False </td>
      <td>Annotate known cancer gene in MAFs</td>
      <td>ex: -kcga</td>
   </tr>
   <tr>
      <td>--tumor_mutation_burden, -tmb</td>
      <td>[]</td>
      <td>Calculate tumor mutation burden for each sample and generate a summary table</td>
      <td>ex: -tmb 60456963  ;  the value indicates the sequencing taget region</td>
   </tr>
   <tr>
      <td>--comut_analysis, -cm</td>
      <td>False </td>
      <td>Summarize mutation counts and mutation types of genes for each sample</td>
      <td>ex: -cm</td>
   </tr>
   <tr>
      <td>--comut_plot, -cmp</td>
      <td>[]</td>
      <td>Generate CoMut plot for input cohort</td>
      <td>ex: -cmp examples/tsv/comut.tsv examples/tsv/comut_info.tsv 0 comut.pdf  :  input of CoMut plot and ouput of the figure with selected format </td>
   </tr>
   <tr>
      <td>--mutational_signature, -ms</td>
      <td>[]</td>
      <td>Estimate optimal number of sinatures, extract mutational signatures for input cohort, and perform visualization</td>
      <td>ex: two steps;  -ms 1 "[2,9,10]"  :  (step 1) estimate optimal number of mutational signatures (2-9) with 10 epoch; -ms 2 "[3]"  : (step 2) Generate plots  according to extracted 3 mutational signatures</td>
   </tr>
   <tr>
      <td>--hrd_score, -hrd :</td>
      <td>[]</td>
      <td>Calculate HRD score for each sample, generate a summary table,  and perform visualization</td>
      <td>ex: -hrd examples/tsv/hrd.tsv grch37  : enter the summarized CNAs of input cohort with information of  genome reference that users employed</td>
   </tr>
   <tr>
      <td>--wgd_cin, -wgdcin</td>
      <td>NULL</td>
      <td>Calculate CIN level for each sample, identify WGD cohort, generate a summary table,  and perform visualization</td>
      <td>ex: -wgdcin examples/tsv/hrd.tsv  : enter the summarized CNAs of input cohort </td>
   </tr>
   <tr>
      <td>--oncokb_annotator, -oncokb</td>
      <td>[]</td>
      <td>Annotate actionable mutation(drug) and perform visualization</td>
      <td>ex: -oncokb ../oncokb-annotator/  [your_oncokb_token]   4   examples/test_data/oncokb/clinical_input.txt  : enter the path of oncokb-annotator, your personal API token of  OncoKB,  choose evidence levels (4 for Level 1-4, 3 for Level 1-3, etc.) for visualization, and the path of  clinical data of input cohort </td>
   </tr>
</table>


#### Column information of input TSV file (for CoMut plot; -cmp command):
##### Input 1
<table>
   <tr>
      <td>Columns</td>
      <td>Descriptions</td>
      <td>Data preparation</td>
      <td>In the detailed format please refer to examples below</td>
   </tr>
   <tr>
      <td>Same Patient</td>
      <td>The path of information for each paired sample</td>
      <td>User-defined TSV file </td>
      <td>examples/test_data/comut/sp.tsv</td>
   </tr>
   <tr>
      <td>Copy Number Alteration </td>
      <td>The path of CNA data for each sample</td>
      <td>User-defined TSV file from the output of CNA calling</td>
      <td>examples/test_data/comut/cna.tsv</td>
   </tr>
   <tr>
      <td>Mutation Type</td>
      <td>The path of mutation types of genes for each sample</td>
      <td>Mutscape output of CoMut plot analysis (-cm)</td>
      <td>examples/test_data/comut/mutation_data.tsv</td>
   </tr>
   <tr>
      <td>Purity </td>
      <td>The path of tumor purity for each sample</td>
      <td>User-defined TSV file from the ouput of purity estimation </td>
      <td>examples/test_data/comut/purity.tsv</td>
   </tr>
   <tr>
      <td>Mutation Signature</td>
      <td>The path of singature contributions for each sample</td>
      <td>Mutscape output of mutational signature analysis (-ms)</td>
      <td>examples/test_data/comut/sig_contribution.tsv</td>
   </tr>
   <tr>
      <td>Mutation Classification</td>
      <td>The path of mutation counts for each sample</td>
      <td>Mutscape output of CoMut plot analysis (-cm)</td>
      <td>examples/test_data/comut/mutation_count.tsv</td>
   </tr>
   <tr>
      <td>Frequency</td>
      <td>The path of mutation frequency of genes for each sample</td>
      <td>User-defined TSV file </td>
      <td></td>
   </tr>
   <tr>
      <td>Whole Genome Doubling</td>
      <td>The path of information of identified WGD samples</td>
      <td>Mutscape output of WGD analysis (-wgdcin)</td>
      <td>examples/test_data/comut/wgd.tsv</td>
   </tr>
</table>

##### Input 2
<table>
   <tr>
      <td>Columns</td>
      <td>Descriptions</td>
      <td>Data preparation</td>
   </tr>
   <tr>
      <td>Copy Number Alteration</td>
      <td>Selected genes of CNAs for Visualization</td>
      <td>Genes according to GISTIC2 or or known cancer gene annotation from other tools</td>
   </tr>
   <tr>
      <td>Mutation Type</td>
      <td>Selected genes of short variants for Visualization</td>
      <td>Genes according to significantly mutated gene detection or known cancer gene annotation from MutScape</td>
   </tr>
   <tr>
      <td>Mutation Signature</td>
      <td>Name of signatures</td>
      <td>User-defined names according to cosine similarity between identified signatures and COSMIC signatures</td>
   </tr>
</table>