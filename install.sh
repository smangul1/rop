dir=$PWD 
db_url="https://googledrive.com/host/0B_NUyiE86yDwaUxoVjhlSjN5SkE/database.tar"

echo "Databases will be renewed"
echo "Do you wish to install the database?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) wget $db_url; tar-xvf $dir/database.tar; break;;
        No ) exit 1;;
    esac
done