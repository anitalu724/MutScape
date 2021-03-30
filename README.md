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

### Install useful modules

To make sure our code smoothly implement, we need to install several modules that are used in MutScape:

    pip3 install termcolor
    pip3 install tabulate
    pip3 install tqdm
    pip3 install numpy
    pip3 install pandas
    pip3 install vcfpy
    pip3 install pysam
    pip3 install seaborn
    pip3 install oncodriveclust
    pip3 install sklearn
    pip3 install fonttools
    pip3 install PyVCF
    pip3 install comut


### Install vcf2maf module
For the VCF transforms to the MAF, this procedure is implemented by `vcf2maf` utility, which processes variant annotation and transcript prioritization. You can refer to  [this script]((https://github.com/mskcc/vcf2maf)) or just follow the commands below. 
```shell
export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-*
perl vcf2maf.pl --man
perl maf2maf.pl --man

curl -sL https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh -o /tmp/miniconda.sh
sh /tmp/miniconda.sh -bfp $HOME/miniconda3

export PATH=$HOME/miniconda3/bin:$PATH

conda install -qy -c conda-forge -c bioconda -c defaults ensembl-vep==100.2 samtools==1.9 bcftools==1.9 ucsc-liftover==377

vep_install --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CACHEDIR $HOME/.vep

curl -sLO https://raw.githubusercontent.com/Ensembl/ensembl-vep/release/100/examples/homo_sapiens_GRCh38.vcf
vep --species homo_sapiens --assembly GRCh38 --offline --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir $HOME/.vep --fasta $HOME/.vep/homo_sapiens/100_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --input_file homo_sapiens_GRCh38.vcf --output_file homo_sapiens_GRCh38.vep.vcf --polyphen b --af --af_1kg --af_esp --regulatory
```

### Install MutScape
For the necessity of definite file path, you should install MutScape in `mskcc-vcf2maf-958809e` folder.
```shell
git clone https://github.com/anitalu724/MutScape.git
```


## Implementation
MutScape has simply separated into two main modules: data preprocessing and analysis and visualization. 

### Data Preprocessing



