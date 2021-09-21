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
