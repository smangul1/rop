#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# INTRO
# ------------------------------------------------------------------------------

set -e
echo '--------------------------------------------------------------------------------'
echo 'Read Origin Protocol: Installer'
echo '--------------------------------------------------------------------------------'
DIR=`dirname $(readlink -f "$0")`
sed '/##/ q' "$DIR/README.md" | head -n -2 | tail -n +3
echo '--------------------------------------------------------------------------------'

# ------------------------------------------------------------------------------
# CONSTANTS
# ------------------------------------------------------------------------------

declare -A DB_ID_HUMAN=(
    ['viral_vipr']='1fIxhnwNSPj6NL2R44bqYkYu2T8OLfqpk'
    ['fungi']='1yBeBjnrnHtxZliruu3oC8NjZ3wHg-WQ3'
    ['BWAindex']='19Uscw8KrPyUiuPcErrXpyOxN0PqUtbOZ'
    ['protozoa']='1_dPn8kk3I--Icy0gwTorFneV1sor1dU2'
    ['metaphlan']='15UGuZ4klBjIEYV-tv6t1nYa2GdyadZAm'
    ['repeats']='1rK1m7sWbiG2cTahLY5SDuDe26KvCWSPL'
    ['ribosomal.DNA']='1kFI5waihEpoZE8DBNYhlfCwJXGsOf4NS'
    ['viral']='1HfnoEhoYzlvo4f6Ap_LZvIrxqzIiKrv-'
)

declare -A DB_ID_MOUSE=(
    ['viral_vipr']='1fIxhnwNSPj6NL2R44bqYkYu2T8OLfqpk'  # same as human
    ['fungi']='1yBeBjnrnHtxZliruu3oC8NjZ3wHg-WQ3'  # same as human
    ['BWAindex']='17g26AzgYhvpY5J5Spofs_rVwXVgVxtrL'
    ['protozoa']='1_dPn8kk3I--Icy0gwTorFneV1sor1dU2'  # same as human
    ['metaphlan']='15UGuZ4klBjIEYV-tv6t1nYa2GdyadZAm'  # same as human
    ['repeats']='1f4xhA8Ku0bMQagY78qDpjb5-6H2AkxIp'
    ['ribosomal.DNA']='1AnzuhTtEGU8_QnHu9ueofSAJoXK-Sg7B'
    ['viral']='1HfnoEhoYzlvo4f6Ap_LZvIrxqzIiKrv-'  # same as human
)

declare -A DB_MD5_HUMAN=(
    ['viral_vipr']='9dce447328dfbc3a62cc7dd5b052242f'
    ['fungi']='9f2d304fd5c49981682b2bb7a900a30e'
    ['BWAindex']='4f009e3732d9f513e7b19b58edc41c13'
    ['protozoa']='23e12115a5e9d526553c901e772731f5'
    ['metaphlan']='3c9b9d6414d86a0c3d5018aefa5aaec4'
    ['repeats']='109a97423f505b73a7e55995b827e2fd'
    ['ribosomal.DNA']='9663a0e1121a0b122c8e23f41c558083'
    ['viral']='7ce95144827603a64dc5996aa0112cc0'
)

declare -A DB_MD5_MOUSE=(
    ['viral_vipr']='9dce447328dfbc3a62cc7dd5b052242f'  # same as human
    ['fungi']='9f2d304fd5c49981682b2bb7a900a30e'  # same as human
    ['BWAindex']='a5166f9e134582a66030878b0b91aed2'
    ['protozoa']='23e12115a5e9d526553c901e772731f5'  # same as human
    ['metaphlan']='3c9b9d6414d86a0c3d5018aefa5aaec4'  # same as human
    ['repeats']='223777c32bfdf8cfe43ca8a01499d4ee'
    ['ribosomal.DNA']='971401e7eae7a721097ae91cb7237b46'
    ['viral']='7ce95144827603a64dc5996aa0112cc0'  # same as human
)

# ------------------------------------------------------------------------------
# PARSE OPTIONS
# ------------------------------------------------------------------------------

# Test for getopt availability.
set +e
getopt --test
if [ $? -ne 4 ]; then
    echo "Error: Environment doesn't support getopt." >&2
    exit 1
fi
set -e

# Call getopt.
SHORT_OPTIONS='cfnl:d:o:rs:h'
LONG_OPTIONS='clean,force,native,link:,db-dest:,organism:,reinstall,select-db:,help'
set +e
PARSED=`getopt --options="$SHORT_OPTIONS" --longoptions="$LONG_OPTIONS" --name "$0" -- "$@"`
if [ $? -ne 0 ]; then
    exit 1  # getopt will have printed the error message
fi
set -e
eval set -- "$PARSED"

# Set default options.
CLEAN_ONLY=false
FORCE=false
NATIVE=false
LINK=''
DB_DEST="$DIR"
ORGANISM='human'
REINSTALL=false
SELECT_DB='all'

# Review parsed options. -c and -l override conflicting options.
while true; do
    case "$1" in
        -c|--clean-only)
            # Just remove installed tools.
            CLEAN_ONLY=true
            shift
            ;;
        -f|--force)
            # Unlink databases.
            FORCE=true
            shift
            ;;
        -n|--native)
            # Use native python.
            NATIVE=true
            shift
            ;;
        -l|--link)
            # Link databases instead of downloading.
            LINK="$2"
            shift 2
            ;;
        -d|--db-dest)
            # Change database download location.
            DB_DEST="$2"
            shift 2
            ;;
        -o|--organism)
            # Organism to download databases for.
            ORGANISM="$2"
            shift 2
            ;;
        -r|--reinstall)
            # Reinstall tools, even if they're already present.
            REINSTALL=true
            shift
            ;;
        -s|--select-db)
            # Database(s) to download for the specified organism.
            SELECT_DB=`tr ',' ' ' <<<"$2"`
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [-cfnrh] [-l LINK] [-d DB_DEST] [-o ORGANISM]"\
                '[-s SELECT_DB]' >&2
            exit 0
            ;;
        --)
            # End of options.
            shift
            break
            ;;
        *)
            echo "Error parsing options." >&2
            exit 1
            ;;
    esac
done

# ------------------------------------------------------------------------------
# DOWNLOAD TOOLS
# ------------------------------------------------------------------------------

cd "$DIR/tools"

# Skip this section if neither -c nor -r are selected and there is a previous
# installation (as indicated by the presence of the imrep directory).
echo '----- Checking for existing installations --------------------------------------'
if [ $CLEAN_ONLY = false ] && [ $REINSTALL = false ] && [ -d 'imrep' ]; then
    echo 'Existing installation found. Skipping tools download. To reinstall,' \
        'please use the -r option.'
else
    echo '----- Removing previous versions -----------------------------------------------'
    rm -fr imrep metaphlan2 MiniConda
    if [ $CLEAN_ONLY = true ]; then
        echo 'Done: Cleaning complete.'
        exit 0
    fi

    # Download ImReP.
    echo '----- Downloading ImRep --------------------------------------------------------'
    git clone https://github.com/mandricigor/imrep.git
    cd imrep
    ./install.sh
    cd ..

    # Download MetaPhlAn 2.
    #echo '----- Downloading MetaPhlAn 2 --------------------------------------------------'
    #hg clone https://bitbucket.org/biobakery/metaphlan2
    #cd metaphlan2
    #ln -s ../../db_human/databases
    #cd ..

    # Download MiniConda and add shebangs.
    echo '----- Setting up Python environment --------------------------------------------'
    if [ $NATIVE = false ]; then
        ./install-MiniConda.sh
        cd MiniConda/lib
        ln -s libncursesw.so.5 libtinfow.so.5
        cd ../..
        MiniConda="$PWD/MiniConda/bin/python"
    #    sed -i "1c #!$MiniConda" metaphlan2/metaphlan2.py
    #    sed -i "1c #!$MiniConda" metaphlan2/strainphlan.py
    #    sed -i "1c #!$MiniConda" metaphlan2/utils/read_fastx.py
    #else
    #    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/metaphlan2.py
    #    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/strainphlan.py
    #    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/utils/read_fastx.py
    fi
fi

# ------------------------------------------------------------------------------
# LINK/UNLINK DATABASES
# ------------------------------------------------------------------------------

cd "$DIR"

echo '----- Checking for existing databases ------------------------------------------'
if [ -h "db_$ORGANISM" ] || [ -d "db_$ORGANISM" ]; then
    if [ $FORCE = true ]; then
        echo 'Unlinking existing database.'
        if [ -h "db_$ORGANISM" ]; then
            rm "db_$ORGANISM"
        else
            rm -r "db_$ORGANISM"
        fi
    else
        echo 'Existing database found. Skipping database download.' \
            'To unlink the current database, please use the -f option.'
        exit 0
    fi
fi

if [ "$LINK" != '' ]; then
    echo '----- Linking database -----------------------------------------------------'
    if [ -d "$LINK" ]; then
        ln -s "$LINK"
        echo 'Done: Database linked.'
        exit 0
    else
        echo "Error: Link target doesn't exist." >&2
        exit 1
    fi
fi

# ------------------------------------------------------------------------------
# DOWNLOAD DATABASES
# ------------------------------------------------------------------------------

cd "$DB_DEST"
mkdir "db_$ORGANISM"
cd "db_$ORGANISM"

download_list=$'ribosomal.DNA\nBWAindex'
for database in $SELECT_DB; do
    case "$database" in
        basic)
            ;;
        repeats)
            download_list+=$'\nrepeats'
            ;;
        microbiome)
            download_list+=$'\nmetaphlan\nviral\nviral_vipr\nfungi\nprotozoa'
            ;;
        metaphlan)
            #download_list+=$'\nmetaphlan'
            ;;
        viral)
            download_list+=$'\nviral\nviral_vipr'
            ;;
        fungi)
            download_list+=$'\nfungi'
            ;;
        protozoa)
            download_list+=$'\nprotozoa'
            ;;
        all)
            download_list+=$'\nrepeats'
            download_list+=$'\nmetaphlan\nviral\nviral_vipr\nfungi\nprotozoa'
            ;;
        *)
            echo 'Error: Unknown database.' >&2
            exit 1
            ;;
    esac
done
download_list=`echo "$download_list" | sort -u`

echo '----- Downloading databases ----------------------------------------------------'
for download in $download_list; do
    echo "Downloading item: $download for $ORGANISM"
    success=false
    while [ $success = false ]; do
        case "$ORGANISM" in
            human)
                db_id="${DB_ID_HUMAN[$download]}"
                db_md5="${DB_MD5_HUMAN[$download]}"
                ;;
            mouse)
                db_id="${DB_ID_MOUSE[$download]}"
                db_md5="${DB_MD5_MOUSE[$download]}"
                ;;
            *)
                echo 'Error: Unknown ORGANISM.' >&2
                exit 1
                ;;
        esac
        confirm_code=`curl --silent --insecure --cookie-jar cookies.txt \
            "https://docs.google.com/uc?export=download&id=$db_id" \
            | sed -rn 's .*confirm=([0-9A-Za-z_]+).* \1\n p'`
        curl --location --insecure --cookie cookies.txt -o "$download.tar.gz" \
            "https://docs.google.com/uc?export=download&confirm=$confirm_code&id=$db_id"
        rm cookies.txt
        if [ `md5sum "$download.tar.gz" | sed 's \(.*\)\ .* \1 '` = "$db_md5" ]; then
            tar -zxvf "$download.tar.gz"
            rm "$download.tar.gz"
            success=true
        else
            echo "Download of $download for $ORGANISM failed (checksum" \
                'mismatch. Retrying.'
        fi
    done
done

cd "$DIR"
if [ `readlink -e "$DB_DEST"` != "$DIR" ]; then
    ln -s "$DB_DEST/db_$ORGANISM"
fi
echo "Done: Reference databases are ready"
