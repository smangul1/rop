#!/bin/bash

# remove the generated files, if any
echo "Removing files"
rm -f tree.txt annot.txt outtree.txt outimg*.png

# convert it!
echo "Converting to GraPhlAn"
export2graphlan.py -i lefse_input.txt -o lefse_output.txt -t tree.txt -a annot.txt --title "HMP Aerobiosis" --annotations 2,3 --external_annotations 4,5,6 --fname_row 0 --skip_rows 1,2 --ftop 200

# attach annotation to the tree
echo "Running Graphlan"
graphlan_annotate.py --annot annot.txt tree.txt outtree.txt

# generate the beautiful image
graphlan.py --dpi 150 outtree.txt outimg.png --external_legends

