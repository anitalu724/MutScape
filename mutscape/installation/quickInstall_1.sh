echo -e "\e[1;33m \nThis is quickInstall_1.sh\n \e[0m"

echo -e "\e[1;33m \nBefore you implement this bash script, please confirm you have already install Miniconda3 (best in version py37_4.9.2).\n \e[0m"

if [[ $key = "" ]]; then 
    echo 'You pressed enter!'
else
    echo "You pressed '$key'"
fi

echo -i "Press enter to continue..."

