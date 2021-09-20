
echo -e "\e[1;33m \nWARNING: Before you implement this bash script, please confirm you have already installed Miniconda3 (best in version py37_4.9.2), created a new environment and also activate it.\n \e[0m"
echo -e "\e[1;33m Press enter to continue... \e[0m" 

read -n 1 input
if [ "$input" = "" ]
then
    echo -e "\e[1;35m Install Ensembl's VEP \n \e[0m"
    conda install -c bioconda -c conda-forge samtools=1.10 ucsc-liftover=377 bcftools=1.10.2
    conda install -c bioconda -c conda-forge -c defaults ensembl-vep=102.0 

    echo -e "Install vcf2maf and move MutScape tool into vcf2maf"
    wget https://github.com/mskcc/vcf2maf/archive/refs/tags/v1.6.20.tar.gz
    tar -zxf v1.6.20.tar.gz
    mv MutScape vcf2maf-1.6.20/
    cd vcf2maf-1.6.20

else
    echo "Quit"
fi


