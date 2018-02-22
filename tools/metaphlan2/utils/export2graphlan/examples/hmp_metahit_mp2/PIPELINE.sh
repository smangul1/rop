# !/bin/bash

# remove the generated files, if any
echo "Removing files"
rm -f tree.txt annot.txt outtree.txt outimg*.png

# convert it!
echo "Converting to GraPhlAn"
export2graphlan.py -i merge-very-good.txt -o merge-very-good.txt.out -t tree.txt -a annot.txt --title "MetaHIT vs. HMP (MetaPhlAn2)" --max_clade_size 250 --min_clade_size 40 --annotations 5 --external_annotations 6,7 --abundance_threshold 40.5 --fname_row 0 --ftop 200 --annotation_legend_font_size 11

# attach annotation to the tree
echo "Running Graphlan"
graphlan_annotate.py --annot annot.txt tree.txt outtree.txt

# generate the beautiful image
graphlan.py --dpi 300 --size 7.0 outtree.txt outimg.png --external_legends

