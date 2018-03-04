#!/bin/bash

# remove the generated files, if any
echo "Removing files"
rm -f tree.txt annot.txt outtree.txt outimg*.png

# convert it!
echo "Converting to GraPhlAn"
export2graphlan.py -i merge.txt -o merge.txt.out -t tree.txt -a annot.txt --title "Metabolic pathways" --abundance_threshold 50.0 --external_annotations 3 --background_clades "Metabolism.Metabolism_of_Cofactors_and_Vitamins, Metabolism.Carbohydrate_Metabolism, Metabolism.Amino_Acid_Metabolism, Metabolism.Metabolism_of_Terpenoids_and_Polyketides, Metabolism.Metabolism_of_Other_Amino_Acids, Genetic_Information_Processing.Replication_and_Repair, Environmental_Information_Processing.Membrane_Transport" --background_colors "(150.; 100.; 100.), (55.; 100.; 100.), (280.; 80.; 88.)" --ftop 125

# attach annotation to the tree
echo "Running Graphlan"
graphlan_annotate.py --annot annot.txt tree.txt outtree.txt

# generate the beautiful image
graphlan.py --dpi 300 --size 7.0 outtree.txt outimg.png --external_legends

