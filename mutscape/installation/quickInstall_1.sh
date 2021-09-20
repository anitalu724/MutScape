
echo -e "\e[1;33m \nWARNING: Before you implement this bash script, please confirm you have already install Miniconda3 (best in version py37_4.9.2).\n \e[0m"
echo -e "\e[1;33m Press enter to continue... \e[0m" 

read -n 1 input
if [ "$input" = "" ]
then
    echo -e "\e[1;35m Create a new Conda environment \n \e[0m"
    conda create --name MutScape
    activate MutScape
else
    echo "Quit"
fi


