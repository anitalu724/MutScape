
echo -e "\e[1;33m \nWARNING: Before you implement this bash script, please confirm you have already installed Miniconda3 (best in version py37_4.9.2), created a new environment and also activate it.\n \e[0m"
echo -e "\e[1;33mPress enter to continue... \e[0m" 

read -n 1 input
if [ "$input" = "" ]
then
    
    echo -e "\e[1;35m \nInstall vcf2maf \n \e[0m"
    wget https://github.com/mskcc/vcf2maf/archive/refs/tags/v1.6.20.tar.gz
    tar -zxf v1.6.20.tar.gz
    rm -rf v1.6.20.tar.gz

    echo -e "\e[1;35m \nMove MutScape tool into vcf2maf\n \e[0m"
    mv MutScape vcf2maf-1.6.20/
    
    echo -e "\e[1;35m \nInstall pip...\n \e[0m"
    conda install -c anaconda pip
    cd vcf2maf-1.6.20/MutScape/mutscape

    echo -e "\e[1;35m \nInstall modules needed for MutScape...\n \e[0m"
    echo -e "\e[1;35m \nInstall termcolor...\n \e[0m"
    pip3 install termcolor
    echo -e "\e[1;35m \nInstall tabulate...\n \e[0m"
    pip3 install tabulate
    echo -e "\e[1;35m \nInstall tqdm...\n \e[0m"
    pip3 install tqdm
    echo -e "\e[1;35m \nInstall numpy...\n \e[0m"
    pip3 install numpy
    echo -e "\e[1;35m \nInstall pandas...\n \e[0m"
    pip3 install pandas
    echo -e "\e[1;35m \nInstall vcfpy...\n \e[0m"
    pip3 install vcfpy
    echo -e "\e[1;35m \nInstall pysam...\n \e[0m"
    pip3 install pysam
    conda config --add channels r
    conda config --add channels bioconda
    conda install pysam
    echo -e "\e[1;35m \nInstall seaborn...\n \e[0m"
    pip3 install seaborn
    echo -e "\e[1;35m \nInstall oncodriveclust...\n \e[0m"
    pip3 install oncodriveclust
    echo -e "\e[1;35m \nInstall sklearn...\n \e[0m"
    pip3 install sklearn
    echo -e "\e[1;35m \nInstall fonttools...\n \e[0m"
    pip3 install fonttools
    echo -e "\e[1;35m \nInstall PyVCF...\n \e[0m"
    pip3 install PyVCF
    echo -e "\e[1;35m \nInstall comut...\n \e[0m"
    pip3 install comut
    echo -e "\e[1;35m \nInstall scarHRD...\n \e[0m"
    conda install -c conda-forge r-devtools
    Rscript installation/scarHRD.r
    echo -e "\e[1;32m \nDone!\n \e[0m"

else
    echo "Quit"
fi


