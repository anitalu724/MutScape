echo -e "\e[1;33m \nThis is quickInstall_1.sh\n \e[0m"

echo -e "\e[1;33m \nWARNING: Before you implement this bash script, please confirm you have already install Miniconda3 (best in version py37_4.9.2).\n \e[0m"

read -t 10 -r -s -p $'Press enter to continue...' input

if [ "input" = "\n" ]
then
    echo "Enter"
else
    echo "quit"
fi


