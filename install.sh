dir=$PWD 
db_url="https://googledrive.com/host/0B_NUyiE86yDwaUxoVjhlSjN5SkE/database.tar"

cd db/
wget $db_url
tar -xvf $dir/db/database.tar
