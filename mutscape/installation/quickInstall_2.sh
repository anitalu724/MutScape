echo -e "\e[1;33m \nWARNING: Before you implement this bash script, please confirm you have already run quickInstall_1.sh. Also, this installation may last for hours, please confirm the Internet connected always.\n \e[0m"
echo -e "\e[1;33mPress enter to continue... \e[0m" 

read -n 1 input
if [ "$input" = "" ]
then
    export PATH="$HOME/miniconda3/bin:$PATH"
    mkdir -p $HOME/.vep/homo_sapiens/102_GRCh37/
    echo -e "\e[1;35m \nDownload Homo_sapiens.GRCh37.dna.toplevel.fa.gz...\n \e[0m"
    wget ftp://ftp.ensembl.org/pub/grch37/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
    mv Homo_sapiens.GRCh37.dna.toplevel.fa.gz  $HOME/.vep/homo_sapiens/102_GRCh37/
    echo -e "\e[1;35m \nDecompress Homo_sapiens.GRCh37.dna.toplevel.fa.gz...\n \e[0m"
    gzip -d $HOME/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
    bgzip -i $HOME/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa
    samtools faidx $HOME/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
    echo -e "\e[1;35m \nDownload homo_sapiens_vep_102_GRCh37.tar.gz...\n \e[0m"
    wget ftp://ftp.ensembl.org/pub/release-102/variation/indexed_vep_cache/homo_sapiens_vep_102_GRCh37.tar.gz
    mv homo_sapiens_vep_102_GRCh37.tar.gz $HOME/.vep/
    echo -e "\e[1;35m \nDecompress homo_sapiens_vep_102_GRCh37.tar.gz...\n \e[0m"
    tar -zxf $HOME/.vep/homo_sapiens_vep_102_GRCh37.tar.gz -C $HOME/.vep/
    echo -e "\e[1;32m \nDone!\n \e[0m"
else
    echo "Quit!"
fi