echo -e "\e[1;33m \nThis is quickInstall_1.sh\n \e[0m"

echo -e "\e[1;33m \nWARNING: Before you implement this bash script, please confirm you have already install Miniconda3 (best in version py37_4.9.2).\n \e[0m"

read -s -n $'Press enter to continue...' input

if [ "$input" = "" ]
then
    echo "Enter"
else
    echo "quit"
fi


