
echo -e "\e[1;33m \nWARNING: Before you implement this bash script, please confirm you have already installed Miniconda3 (best in version py37_4.9.2), created a new environment and also activate it.\n \e[0m"
echo -e "\e[1;33m Press enter to continue... \e[0m" 

read -n 1 input
if [ "$input" = "" ]
then
    echo -e "\e[1;35m \nInstall Ensembl's VEP \n \e[0m"
    conda install -c bioconda -c conda-forge samtools=1.10 ucsc-liftover=377 bcftools=1.10.2
    conda install -c bioconda -c conda-forge -c defaults ensembl-vep=102.0 

    echo -e "\e[1;35m \nInstall vcf2maf \n \e[0m"
    wget https://github.com/mskcc/vcf2maf/archive/refs/tags/v1.6.20.tar.gz
    tar -zxf v1.6.20.tar.gz

    echo -e "\e[1;35m \nMove MutScape tool into vcf2maf\n \e[0m"
    mv MutScape vcf2maf-1.6.20/
    cd vcf2maf-1.6.20
    perl vcf2maf.pl --man

    conda install -c anaconda pip
    cd MutScape/mutscape
    bash installation/install_module.sh

else
    echo "Quit"
fi


