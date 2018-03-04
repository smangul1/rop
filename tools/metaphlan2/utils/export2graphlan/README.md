**export2graphlan** is a conversion software tool for producing both annotation and tree file for GraPhlAn. In particular, the annotation file tries to highlight specific sub-trees deriving automatically from input file what nodes are important. The two output file of **export2graphlan** should then be used to run ``graphlan_annotate.py``, in order to attach to the tree the derived annotations, and finally, by executing ``graphlan.py`` the user can get the output image.

# PREREQUISITES #

**export2graphlan** requires the following additional library:

* pandas ver. 0.13.1 ([pandas](http://pandas.pydata.org/index.html))
* BIOM ver. 2.0.1 ([biom-format](http://biom-format.org), only if you have input files in BIOM format)
* SciPy ([scipy](http://www.scipy.org), required by hclust2)

# INSTALLATION #

**export2graphlan** should be obtained using [Mercurial](http://mercurial.selenic.com/) and is available in Bitbucket here: [export2graphlan repository](https://bitbucket.org/CibioCM/export2graphlan).

In a Unix environment you have to type:
```
#!bash

$ hg clone ssh://hg@bitbucket.org/CibioCM/export2graphlan
```
or, alternatively:
```
#!bash

$ hg clone https://hg@bitbucket.org/CibioCM/export2graphlan
```

This will download the **export2graphlan** repository locally in the ``export2graphlan`` subfolder. You then have to put this subfolder into the system path, so that you can use **export2graphlan** from anywhere in your system:
```
#!bash

$ export PATH=`pwd`/export2graphlan/:$PATH
```
Adding the above line into the bash configuration file will make the path addition permanent. For Windows or MacOS systems a similar procedure should be followed.

# USAGE #
```
#!
usage: export2graphlan.py [-h] [-i LEFSE_INPUT] [-o LEFSE_OUTPUT] -t TREE -a
usage: export2graphlan.py [-h] [-i LEFSE_INPUT] [-o LEFSE_OUTPUT] -t TREE -a
                          ANNOTATION [--annotations ANNOTATIONS]
                          [--external_annotations EXTERNAL_ANNOTATIONS]
                          [--background_levels BACKGROUND_LEVELS]
                          [--background_clades BACKGROUND_CLADES]
                          [--background_colors BACKGROUND_COLORS]
                          [--title TITLE] [--title_font_size TITLE_FONT_SIZE]
                          [--def_clade_size DEF_CLADE_SIZE]
                          [--min_clade_size MIN_CLADE_SIZE]
                          [--max_clade_size MAX_CLADE_SIZE]
                          [--def_font_size DEF_FONT_SIZE]
                          [--min_font_size MIN_FONT_SIZE]
                          [--max_font_size MAX_FONT_SIZE]
                          [--annotation_legend_font_size ANNOTATION_LEGEND_FONT_SIZE]
                          [--abundance_threshold ABUNDANCE_THRESHOLD]
                          [--most_abundant MOST_ABUNDANT]
                          [--least_biomarkers LEAST_BIOMARKERS]
                          [--discard_otus] [--internal_levels] [--sep SEP]
                          [--out_table OUT_TABLE] [--fname_row FNAME_ROW]
                          [--sname_row SNAME_ROW]
                          [--metadata_rows METADATA_ROWS]
                          [--skip_rows SKIP_ROWS] [--sperc SPERC]
                          [--fperc FPERC] [--stop STOP] [--ftop FTOP]
                          [--def_na DEF_NA]

export2graphlan.py (ver. 0.17 of 21th August 2014). Convert MetaPhlAn, LEfSe,
and/or HUMAnN output to GraPhlAn input format. Authors: Francesco Asnicar
(francesco.asnicar@gmail.com)

optional arguments:
  -h, --help            show this help message and exit
  --annotations ANNOTATIONS
                        List which levels should be annotated in the tree. Use
                        a comma separate values form, e.g.,
                        --annotation_levels 1,2,3. Default is None
  --external_annotations EXTERNAL_ANNOTATIONS
                        List which levels should use the external legend for
                        the annotation. Use a comma separate values form,
                        e.g., --annotation_levels 1,2,3. Default is None
  --background_levels BACKGROUND_LEVELS
                        List which levels should be highlight with a shaded
                        background. Use a comma separate values form, e.g.,
                        --background_levels 1,2,3
  --background_clades BACKGROUND_CLADES
                        Specify the clades that should be highlight with a
                        shaded background. Use a comma separate values form
                        and surround the string with " if it contains spaces.
                        Example: --background_clades "Bacteria.Actinobacteria,
                        Bacteria.Bacteroidetes.Bacteroidia,
                        Bacteria.Firmicutes.Clostridia.Clostridiales"
  --background_colors BACKGROUND_COLORS
                        Set the color to use for the shaded background. Colors
                        can be either in RGB or HSV (using a semi-colon to
                        separate values, surrounded with ()) format. Use a
                        comma separate values form and surround the string
                        with " if it contains spaces. Example:
                        --background_colors "#29cc36, (150; 100; 100), (280;
                        80; 88)"
  --title TITLE         If specified set the title of the GraPhlAn plot.
                        Surround the string with " if it contains spaces,
                        e.g., --title "Title example"
  --title_font_size TITLE_FONT_SIZE
                        Set the title font size. Default is 15
  --def_clade_size DEF_CLADE_SIZE
                        Set a default size for clades that are not found as
                        biomarkers by LEfSe. Default is 10
  --min_clade_size MIN_CLADE_SIZE
                        Set the minimum value of clades that are biomarkers.
                        Default is 20
  --max_clade_size MAX_CLADE_SIZE
                        Set the maximum value of clades that are biomarkers.
                        Default is 200
  --def_font_size DEF_FONT_SIZE
                        Set a default font size. Default is 10
  --min_font_size MIN_FONT_SIZE
                        Set the minimum font size to use. Default is 8
  --max_font_size MAX_FONT_SIZE
                        Set the maximum font size. Default is 12
  --annotation_legend_font_size ANNOTATION_LEGEND_FONT_SIZE
                        Set the font size for the annotation legend. Default
                        is 10
  --abundance_threshold ABUNDANCE_THRESHOLD
                        Set the minimun abundace value for a clade to be
                        annotated. Default is 20.0
  --most_abundant MOST_ABUNDANT
                        When only lefse_input is provided, you can specify how
                        many clades highlight. Since the biomarkers are
                        missing, they will be chosen from the most abundant
  --least_biomarkers LEAST_BIOMARKERS
                        When only lefse_input is provided, you can specify the
                        minimum number of biomarkers to extract. The taxonomy
                        is parsed, and the level is choosen in order to have
                        at least the specified number of biomarkers
  --discard_otus        If specified the OTU ids will be discarde from the
                        taxonmy. Default behavior keep OTU ids in taxonomy
  --internal_levels     If specified sum-up from leaf to root the abundances
                        values. Default behavior do not sum-up abundances on
                        the internal nodes

input parameters:
  You need to provide at least one of the two arguments

  -i LEFSE_INPUT, --lefse_input LEFSE_INPUT
                        LEfSe input data
  -o LEFSE_OUTPUT, --lefse_output LEFSE_OUTPUT
                        LEfSe output result data

output parameters:
  -t TREE, --tree TREE  Output filename where save the input tree for GraPhlAn
  -a ANNOTATION, --annotation ANNOTATION
                        Output filename where save GraPhlAn annotation

Input data matrix parameters:
  --sep SEP
  --out_table OUT_TABLE
                        Write processed data matrix to file
  --fname_row FNAME_ROW
                        row number containing the names of the features
                        [default 0, specify -1 if no names are present in the
                        matrix
  --sname_row SNAME_ROW
                        column number containing the names of the samples
                        [default 0, specify -1 if no names are present in the
                        matrix
  --metadata_rows METADATA_ROWS
                        Row numbers to use as metadata[default None, meaning
                        no metadata
  --skip_rows SKIP_ROWS
                        Row numbers to skip (0-indexed, comma separated) from
                        the input file[default None, meaning no rows skipped
  --sperc SPERC         Percentile of sample value distribution for sample
                        selection
  --fperc FPERC         Percentile of feature value distribution for sample
                        selection
  --stop STOP           Number of top samples to select (ordering based on
                        percentile specified by --sperc)
  --ftop FTOP           Number of top features to select (ordering based on
                        percentile specified by --fperc)
  --def_na DEF_NA       Set the default value for missing values [default None
                        which means no replacement]
```

*Note*: the last input parameters (``Input data matrix parameters``) refer to the **DataMatrix** class contained in the [hclust2](https://bitbucket.org/nsegata/hclust2/overview) repository.

# EXAMPLES #
The ``examples`` folder contains the following sub-folders: ``hmp_aerobiosis``, ``hmp_metahit_metabolic``, and ``hmp_metahit_mp2``.
Each example should work just by typing in a terminal window (provided that you are inside one of the example folder) the following command:
```
#!bash

$ ./PIPELINE.sh
```

If everything goes well you should find in the same folder of the example six new files: ``annot.txt``, ``outimg.png``, ``outimg_annot.png``, ``outimg_legend.png``, ``outtree.txt``, and ``tree.txt``. Where:

* ``annot.txt``: contains the annotation that will be used by GraPhlAn, produced by the export2graphlan.py script
* ``outimg.png``: is the circular tree produced by GraPhlAn
* ``outimg_annot.png``: contains the annotation legend of the circular tree
* ``outimg_legend.png``: contains the legends of the highlighted biomarkers in the circular tree
* ``outtree.txt``: is the annotated tree produced by graphlan_annotate.py
* ``tree.txt``: is the tree produced by the export2graphlan.py script

# CONTACTS #
Francesco Asnicar ([francescoasnicar@gmail.com](mailto:francescoasnicar@gmail.com))