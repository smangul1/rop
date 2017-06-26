Please make sure that the basic unix commands (tar or zip) are available on the cluster. The package  can be downloaded as zip or tar.gz compressed archives. In a Unix environment, you can obtain and uncompress it from the command line:

```
wget https://sourceforge.net/projects/rop2/files/latest/download/rop.tar.gz
tar -zxvf rop.tar.gz
```

Now we require installation using `INSTALL.sh`. (No installation is required in versions previous to v1.0.6)



ROP requires Python 2.7. Please ensure that Python 2.7 is available on your system. This `INSTALL.sh` script will install all other prerequisite software (no administrator permissions required). ROP is distributed with several open source components that were developed by other groups (see section Third Party software from [here](https://sergheimangul.wordpress.com//rop/)).


In addition, please use `getDB.py` to download prepared reference database(s), which can be download as shown bellow. 

```
cd rop
python getDB.py ~/
```


If you planning to use ROP for mouse, you need to download `db_mouse`, using the following command:

```
python getDB.py --organism mouse ~/
```
Make sure you have enough (50G) space in `~/` to store the database. 

Where `~/` means that you are downloading the `db` in your home directory. After running `python getDB.py ~/`, please check the `~/db_human/`. You expect to see the following subdirectories in the `db_human` directory:

```
drwxr-xr-x 3 serghei eeskin 20480 Aug 10 13:46 antibody
drwxr-xr-x 2 serghei eeskin  8192 Aug 10 13:50 bacteria
drwxr-xr-x 2 serghei eeskin  4096 Aug 10 13:45 bowtie2Index
drwxr-xr-x 2 serghei eeskin  4096 Aug 10 13:46 BWAIndex
-rw-r--r-- 1 serghei eeskin   377 Sep  5 06:42 cookies.txt
drwxr-xr-x 2 serghei eeskin 16384 Aug 10 13:48 eupathdb
drwxr-xr-x 2 serghei eeskin  4096 Aug 10 13:46 metaphlan
drwxr-xr-x 2 serghei eeskin  4096 Aug 20 21:44 repeats
drwxr-xr-x 2 serghei eeskin  4096 Aug 10 13:44 rRNA
drwxr-xr-x 2 serghei eeskin  4096 Aug 10 13:46 virus
```


To be sure that the db download was successful, please check the rop directory. Given that you are at `rop` directory, run `ls -l`. Make sure you see `db_human -> /u/home/s/serghei/project/db_human/`. This is the link pointing to the `db` you downlaoded. 



In the versions prior to 1.0.3 the `getDB.py` was named as `installation.py`. Upgrading to the latest version is strongly encouraged. 
 
Downloading reference databases (refDB) might take up to  45 minutes on average and requires 50Gb of available space to download the references.

The 'getDB.py' provides the possibility to connect the current version of ROP to the existing database. Please use `--link2db` option. 

Prior to the release 1.0.3 refDB was downloaded to the ROP directory. Now the location to download refDB is a choice of the user. This allows to store the refDB independently from ROP and connect the current version of ROP to the refDB. As you get a new version of ROP there is no need to download refDB again, instead you can just connect the new version of ROP to the existing database (via `--link2db` option).


The full list of parameters of 'getDB.py' can be accessed from the help message:

```
usage: python rop.py [-h] [--organism ORGANISM] [--repeat] [--immune]
                     [--metaphlan] [--circRNA] [--microbiome] [--link2db]
                     [--f]
                     dirDB

optional arguments:
  -h, --help           show this help message and exit

Necessary Inputs:
  dirDB                directory where the reference databases will be
                       downloaded

Select Database (multi-selection is possible):
  --organism ORGANISM  downloads the (mouse) version of the database files
  --repeat             Set up database for repeat sequences ONLY (lost repeat
                       reads)
  --immune             Set up database for VDJ gene segments of B cell and T
                       cell receptors ONLY (immune reads)
  --metaphlan          Set up database for Metaphlan2 ONLY
  --circRNA            Set up database for circular RNAs ONLY
  --microbiome         Set up database for microbiome ONLY

Connect database with the new release:
  --link2db            Connect the existing reference database with ROP
  --f                  Reconnects the ROP to the new database. Please note the
                       existing link will removed
```


