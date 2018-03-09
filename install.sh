#!/usr/bin/env bash

echo '--------------------------------------------------------------------------------'
echo 'Read Origin Protocol: Installer'
echo '--------------------------------------------------------------------------------'
cat README.md
echo '--------------------------------------------------------------------------------'
DIR=`dirname $(readlink -f $0)`
cd $DIR/tools

# Uninstall previous versions. Restore shebangs and exit if the --clean option
# is selected.
rm -fr imrep metaphlan2 MiniConda
if [ $# -eq 1 ] && [ $1 = '--clean' ]; then
    sed -i '1c #!/usr/bin/env python2.7' ../rop.py
    sed -i '1c #!/usr/bin/env python2.7' ../getDB.py
    echo 'Cleaning complete. Please use install.sh to reinstall.'
    exit 0
fi

# Install ImReP.
git clone https://github.com/mandricigor/imrep.git
cd imrep
./install.sh
cd ..

# Install MetaPhlAn 2.
hg clone https://bitbucket.org/biobakery/metaphlan2
cd metaphlan2
ln -s ../../db_human/databases
cd ..

# Install MiniConda and add shebangs.
if [ $# -ne 1 ] || [ $1 != '--no-miniconda' ]; then
    ./install-MiniConda.sh
    MiniConda="$PWD/MiniConda/bin/python"
    sed -i "1c #!$MiniConda" metaphlan2/metaphlan2.py
    sed -i "1c #!$MiniConda" metaphlan2/strainphlan.py
    sed -i "1c #!$MiniConda" metaphlan2/utils/read_fastx.py
    sed -i "1c #!$MiniConda" ../rop.py
    sed -i "1c #!$MiniConda" ../getDB.py
else
    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/metaphlan2.py
    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/strainphlan.py
    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/utils/read_fastx.py
    sed -i '1c #!/usr/bin/env python2.7' ../rop.py
    sed -i '1c #!/usr/bin/env python2.7' ../getDB.py
fi

echo 'Installation complete. Please use getDB.py to download databases.'
