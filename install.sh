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
    ['viral_vipr']='1ISkC9Z4oYgtxlitMbB-C3T9tZgErjR8_'
    ['fungi']='1NH4bx1ZpDpiv_wqYcxwR02kBalVdMF2f'
    ['BWAindex']='13H20BSdp6JMG3BzDSL0JE3m9X12Jaj6n'
    ['protozoa']='1X5ElXNkkTBxd9Wui6ttvi3LhdcfDzTMc'
    ['metaphlan']='1RnhchW875-DqED7xkWV5W1oFFESKLE1r'
    ['repeats']='1W4wMCUSTgVgnUxVocYy1CqCSmEYWJyuF'
    ['ribosomal.DNA']='1T4JdkSJUZVvOZdZcWgmk_u0XSfLdxqdF'
    ['viral']='1ZzCuwW-38yU6igoIJYKwi8KFPwSxcvsg'
)

declare -A DB_ID_MOUSE=(
    ['viral_vipr']='1ISkC9Z4oYgtxlitMbB-C3T9tZgErjR8_'  # same as human
    ['fungi']='1NH4bx1ZpDpiv_wqYcxwR02kBalVdMF2f'  # same as human
    ['BWAindex']='1L9fMlRjGcB85_orar1ObK_VVE3IfbUsd'
    ['protozoa']='1X5ElXNkkTBxd9Wui6ttvi3LhdcfDzTMc'  # same as human
    ['metaphlan']='1RnhchW875-DqED7xkWV5W1oFFESKLE1r'  # same as human
    ['repeats']='1z888jreCH3HSUfbn-DvkBdNdeZ34tLLe'
    ['ribosomal.DNA']='1YtPRdUciVl2K0GybH4VPs_mf-l-LvUsX'
    ['viral']='1ZzCuwW-38yU6igoIJYKwi8KFPwSxcvsg'  # same as human
)

declare -A DB_MD5_HUMAN=(
    ['viral_vipr']='3b5d1e88e8fd64db2fd997283e06c7ea'
    ['fungi']='99dff18fae92549b1b7e92fa3fede4ff'
    ['BWAindex']='1bdaba6c5987a0bdb37657b6fb33a2fc'
    ['protozoa']='33751a0312882dbbd9699cb53d78459d'
    ['metaphlan']='a07a407ce7071df373af9d6795704bc9'
    ['repeats']='c6939256a21dd166c9c38380afbd9f7c'
    ['ribosomal.DNA']='f339943ddfff31533e749dc28b70065c'
    ['viral']='6a16fc9726ebad359a6b3719703dba83'
)

declare -A DB_MD5_MOUSE=(
    ['viral_vipr']='3b5d1e88e8fd64db2fd997283e06c7ea'  # same as human
    ['fungi']='99dff18fae92549b1b7e92fa3fede4ff'  # same as human
    ['BWAindex']='6acd8bb77687dcd139a27d1bcab72eca'
    ['protozoa']='33751a0312882dbbd9699cb53d78459d'  # same as human
    ['metaphlan']='a07a407ce7071df373af9d6795704bc9'  # same as human
    ['repeats']='8988a4b1c671bcf2bbce8a4b08d564ef'
    ['ribosomal.DNA']='78b73b3548feb7113cf56a5fd61d3ce2'
    ['viral']='6a16fc9726ebad359a6b3719703dba83'  # same as human
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
    echo '----- Downloading MetaPhlAn 2 --------------------------------------------------'
    git clone https://github.com/biobakery/metaphlan2
    cd metaphlan2
    ln -s ../../db_human/databases db_v20
    cd ..

    # Download MiniConda.
    echo '----- Setting up Python environment --------------------------------------------'
    if [ $NATIVE = false ]; then
        ./install-MiniConda.sh
        cd MiniConda/lib
        ln -s libncursesw.so.5 libtinfow.so.5
        cd ../..
        MiniConda="$PWD/MiniConda/bin/python"
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
            download_list+=$'\nmetaphlan'
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
        rm -f cookies.txt
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
