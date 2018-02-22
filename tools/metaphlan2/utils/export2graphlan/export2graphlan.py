#!/usr/bin/env python


import os
import numpy as np
from argparse import ArgumentParser
from colorsys import hsv_to_rgb
from math import log10
from StringIO import StringIO
from re import compile
from hclust2.hclust2 import DataMatrix


__author__ = 'Francesco Asnicar'
__email__ = 'f.asnicar@unitn.it'
__version__ = '0.20'
__date__ = '29th May 2017'


pre_taxa = compile(".__")


def scale_color((h, s, v), factor=1.0):
    """
    Takes as input a tuple that represents a color in HSV format, and optionally a scale factor.
    Return an RGB string that is the converted HSV color, scaled by the given factor.
    """
    if (h < 0.) or (h > 360.):
        raise Exception('[scale_color()] Hue value out of range (0, 360): ' + str(h))

    if (s < 0.) or (s > 100.):
        raise Exception('[scale_color()] Saturation value out of range (0, 100): ' + str(s))

    if (v < 0.) or (v > 100.):
        raise Exception('[scale_color()] Value value out of range (0, 100): ' + str(v))

    if (factor < 0.) or (factor > 1.):
        raise Exception('[scale_color()] Factor value out of range (0.0, 1.0): ' + str(factor))

    v *= factor
    r, g, b = hsv_to_rgb(h/360., s/100., v/100.)

    return '#{0:02x}{1:02x}{2:02x}'.format(int(round(r*255.)), int(round(g*255.)), int(round(b*255.)))


def read_params():
    """
    Parse the input parameters, performing some validity check.
    Return the parsed arguments.
    """
    parser = ArgumentParser(description="export2graphlan.py (ver. "+__version__+
        " of "+__date__+"). Convert MetaPhlAn, LEfSe, and/or HUMAnN output to GraPhlAn input format. Authors: "+
        __author__+" ("+__email__+")")

    # input parameters group
    group = parser.add_argument_group(title='input parameters',
        description="You need to provide at least one of the two arguments")
    group.add_argument('-i', '--lefse_input',
        type=str,
        required=False,
        help="LEfSe input data. A file that can be given to LEfSe for biomarkers analysis. It can be the result of a "
             "MetaPhlAn or HUMAnN analysis")
    group.add_argument('-o', '--lefse_output',
        type=str,
        required=False,
        help="LEfSe output result data. The result of LEfSe analysis performed on the lefse_input file")

    # output parameters group
    group = parser.add_argument_group(title='output parameters')
    group.add_argument('-t', '--tree',
        type=str,
        required=True,
        help="Output filename where save the input tree for GraPhlAn")
    group.add_argument('-a', '--annotation',
        type=str,
        required=True,
        help="Output filename where save GraPhlAn annotation")

    # annotations
    parser.add_argument('--annotations',
        default=None,
        type=str,
        required=False,
        help="List which levels should be annotated in the tree. Use a comma separate values form, e.g., "
             "--annotation_levels 1,2,3. Default is None")
    parser.add_argument('--external_annotations',
        default=None,
        type=str,
        required=False,
        help="List which levels should use the external legend for the annotation. Use a comma separate values form, "
             "e.g., --annotation_levels 1,2,3. Default is None")
    # shaded background
    parser.add_argument('--background_levels',
        default=None,
        type=str,
        required=False,
        help="List which levels should be highlight with a shaded background. Use a comma separate values form, e.g., "
             "--background_levels 1,2,3. Default is None")
    parser.add_argument('--background_clades',
        default=None,
        type=str,
        required=False,
        help="Specify the clades that should be highlight with a shaded background. Use a comma separate values form "
             "and surround the string with \" if there are spaces. Example: --background_clades \"Bacteria.Actinobacteria, "
             "Bacteria.Bacteroidetes.Bacteroidia, Bacteria.Firmicutes.Clostridia.Clostridiales\". Default is None")
    parser.add_argument('--background_colors',
        default=None,
        type=str,
        required=False,
        help="Set the color to use for the shaded background. Colors can be either in RGB or HSV (using a semi-colon to "
             "separate values, surrounded with ()) format. Use a comma separate values form and surround the string with "
             "\" if it contains spaces. Example: --background_colors \"#29cc36, (150; 100; 100), (280; 80; 88)\". To use "
             "a fixed set of colors associated to a fixed set of clades, you can specify a mapping file in a tab-separated "
             "format, where the first column is the clade (using the same format as for the \"--background_clades\" param) "
             "and the second colum is the color associated. Default is None")
    # title
    parser.add_argument('--title',
        type=str,
        required=False,
        help="If specified set the title of the GraPhlAn plot. Surround the string with \" if it contains spaces, e.g., "
             "--title \"Title example\"")
    # title font size
    parser.add_argument('--title_font_size',
        default=15,
        type=int,
        required=False,
        help="Set the title font size. Default is 15")
    # clade size
    parser.add_argument('--def_clade_size',
        default=10.,
        type=float,
        required=False,
        help="Set a default size for clades that are not found as biomarkers by LEfSe. Default is 10")
    parser.add_argument('--min_clade_size',
        default=20.,
        type=float,
        required=False,
        help="Set the minimum value of clades that are biomarkers. Default is 20")
    parser.add_argument('--max_clade_size',
        default=200.,
        type=float,
        required=False,
        help="Set the maximum value of clades that are biomarkers. Default is 200")
    # font size
    parser.add_argument('--def_font_size',
        default=10,
        type=int,
        required=False,
        help="Set a default font size. Default is 10")
    parser.add_argument('--min_font_size',
        default=8,
        type=int,
        required=False,
        help="Set the minimum font size to use. Default is 8")
    parser.add_argument('--max_font_size',
        default=12,
        type=int,
        required=False,
        help="Set the maximum font size. Default is 12")
    # legend font size
    parser.add_argument('--annotation_legend_font_size',
        default=10,
        type=int,
        required=False,
        help="Set the font size for the annotation legend. Default is 10")
    # abundance threshold
    parser.add_argument('--abundance_threshold',
        default=20.,
        type=float,
        required=False,
        help="Set the minimun abundace value for a clade to be annotated. Default is 20.0")
    # ONLY lefse_input provided
    parser.add_argument('--most_abundant',
        default=10,
        type=int,
        required=False,
        help="When only lefse_input is provided, you can specify how many clades highlight. Since the biomarkers are "
             "missing, they will be chosen from the most abundant. Default is 10")
    parser.add_argument('--least_biomarkers',
        default=3,
        type=int,
        required=False,
        help="When only lefse_input is provided, you can specify the minimum number of biomarkers to extract. The "
             "taxonomy is parsed, and the level is choosen in order to have at least the specified number of biomarkers. "
             "Default is 3")
    # decide to keep the OTU id or to merger at the above taxonomic level
    parser.add_argument('--discard_otus',
        default=True,
        action='store_false',
        help="If specified the OTU ids will be discarde from the taxonmy. Default is True, i.e. keep OTUs IDs in taxonomy")
    # decide to keep the OTU id or to merger at the above taxonomic level
    parser.add_argument('--internal_levels',
        default=False,
        action='store_true',
        help="If specified sum-up from leaf to root the abundances values. Default is False, i.e. do not sum-up abundances "
             "on the internal nodes")
    # path to a mapping file that associates biomarkers to colors
    parser.add_argument('--biomarkers2colors',
    default=None,
    type=str,
    required=False,
    help="Mapping file that associates biomarkers to a specific color... I'll define later the specific format of this file!")

    DataMatrix.input_parameters(parser)
    args = parser.parse_args()

    # check if at least one of the input params is given
    if (not args.lefse_input) and (not args.lefse_output):
        raise Exception("[read_params()] You must provide at least one of the two input parameters: ")

    # check that min_clade_size is less than max_clade_size
    if args.min_clade_size > args.max_clade_size:
        print "[W] min_clade_size cannot be greater than max_clade_size, assigning their default values"
        args.min_clade_size = 20.
        args.max_clade_size = 200.

    # check that min_font_size is less than max_font_size
    if args.min_font_size > args.max_font_size:
        print "[W] min_font_size cannot be greater than max_font_size, assigning their default values"
        args.min_font_size = 8
        args.max_font_size = 12

    return args


def get_file_type(filename):
    """
    Return the extension (if any) of the ``filename`` in lower case.
    """
    return filename[filename.rfind('.')+1:].lower()


def parse_biom(filename, keep_otus=True, internal_levels=False):
    """
    Load a biom table and extract the taxonomy (from metadata), removing the unuseful header.
    Return the input biom in tab-separated format.
    """
    from biom import load_table # avoid to ask for the BIOM library if there is no biom file

    biom_table = load_table(filename)
    strs = biom_table.delimited_self(header_value='TAXA', header_key='taxonomy')
    lst1 = [str(s) for s in strs.split('\n')[1:]] # skip the "# Constructed from biom file" entry
    biom_file = []
    out = [lst1[0]] # save the header
    # pre_taxa = compile(".__")
    classs = compile("\(class\)")

    # consistency check
    i = 0
    while i < (len(lst1)-1):
        if len([s for s in lst1[i].split('\t')]) != len([s for s in lst1[i+1].split('\t')]):
            raise Exception('[parse_biom()] It seems that taxonomic metadata are missing, maybe is the wrong biom file?')

        i += 1

    for l in lst1[1:]:
        otu = None
        lst = [float(s.strip()) for s in l.split('\t')[1:-1]]

        if keep_otus:
            otu = l.split('\t')[0]

        # Clean and move taxa in first place
        taxa = '.'.join([s.strip().replace('[', '').replace('u\'', '').replace(']', '').replace(' ', '').replace('\'', '')
                         for s in l.split('\t')[-1].split(',')])
        taxa = pre_taxa.sub('', taxa) # remove '{k|p|c|o|f|g|s|t}__'
        taxa = classs.sub('', taxa) # remove '(class)'
        taxa = taxa.rstrip('.') # remove trailing dots

        if otu:
            taxa = taxa + '.' + otu

        biom_file.append([taxa] + lst)

    # merge such rows that have the same taxa
    i = 1
    dic = {}

    for l in biom_file[i:]:
        if l[0] not in dic:
            dic[l[0]] = l[1:]

            for k in biom_file[i+1:]:
                if l[0] == k[0]:
                    lst = []
                    lstdic = dic[l[0]]
                    j = 1
                    while j < len(lstdic):
                        lst.append(float(lstdic[j]) + float(k[j]))
                        j += 1

                    dic[l[0]] = lst
        i += 1

    feats = dict(dic)

    if internal_levels:
        feats = add_missing_levels(feats)

    for k in feats:
        out.append('\t'.join([str(s) for s in [k] + feats[k]]))

    return '\n'.join(out)


def add_missing_levels(ff, summ=True):
    """
    Sum-up the internal abundances from leaf to root
    """
    if sum([f.count(".") for f in ff]) < 1:
        return ff

    clades2leaves = {}
    for f in ff:
        fs = f.split(".")

        if len(fs) < 2:
            continue

        for l in range(1, len(fs)+1):
            n = ".".join(fs[:l])

            if n in clades2leaves:
                clades2leaves[n].append(f)
            else:
                clades2leaves[n] = [f]

    ret = {}
    for k in clades2leaves:
        if summ:
            ret[k] = [sum([sum(ff[e]) for e in clades2leaves[k]])]
        else:
            lst = []
            for e in clades2leaves[k]:
                if not lst:
                    for i in ff[e]:
                        lst.append(i)
                else:
                    lst1 = []
                    i = 0
                    while i < len(lst):
                        lst1.append(lst[i] + ff[e][i])
                        i += 1
                    lst = lst1

            ret[k] = lst

    return ret


def get_most_abundant(abundances, xxx):
    """
    Sort by the abundance level all the taxonomy that represent at least two levels.
    Return the first ``xxx`` most abundant.
    """
    abundant = []

    for a in abundances:
        if a.count('|') > 0:
            abundant.append((float(abundances[a]), a.replace('|', '.')))
        elif a.count('.') > 0:
            abundant.append((float(abundances[a]), a))

    abundant.sort(reverse=True)
    return abundant[:xxx]


def get_biomarkes(abundant, xxx):
    """
    Split the taxonomy and then look, level by level, when there are at least ``xxx`` distinct branches.
    Return the set of branches as biomarkers to highlight.
    """
    cc = []
    old_bk = set()
    lvl = 0

    for _, t in abundant:
        cc.append(t.split('.'))

    while lvl < len(max(cc)):
        bk = set()

        for c in cc:
            if lvl < len(c):
                bk |= set([c[lvl]])

        if len(bk) >= xxx:
            break

        if len(old_bk) > len(bk):
            bk = old_bk
            break

        old_bk = bk
        lvl += 1

    return bk

def scale_clade_size(minn, maxx, abu, max_abu):
    """
    Return the value of ``abu`` scaled to ``max_abu`` logarithmically, and then map from ``minn`` to ``maxx``.
    """
    return minn + (maxx-minn) * log10(1. + 9. * (abu/max_abu))


def main():
    """
    """
    # HSV
    colors = [(245., 90., 100.), # blue
              (125., 80., 80.), # green
              (0., 80., 100.), # red
              (195., 100., 100.), # cyan
              (150., 100., 100.), # light green
              (55., 100., 100.), # yellow
              (280., 80., 88.)] # purple
    args = read_params()
    lefse_input = None
    lefse_output = {}
    color = {}
    biomarkers = set()
    taxa = []
    abundances = {}
    max_abundances = None
    max_effect_size = None
    max_log_effect_size = None
    background_list = []
    background_clades = []
    background_colors = {}
    annotations_list = []
    external_annotations_list = []
    lin = False
    lout = False

    # get the levels that should be shaded
    if args.background_levels:
        background_list = [int(i.strip()) for i in args.background_levels.strip().split(',')]

    # get the background_clades
    if args.background_clades:
        if get_file_type(args.background_clades) in ['txt']:
            with open(args.background_clades, 'r') as f:
                background_clades = [str(s.strip()) for s in f]
        else: # it's a string in csv format
            background_clades = [str(s.strip()) for s in args.background_clades.split(',')]

    # read the set of colors to use for the background_clades
    if args.background_colors:
        col = []

        if get_file_type(args.background_colors) in ['txt']: # it's a mapping file
            background_colors = dict([tuple(a.strip().split('\t')) for a in open(args.background_colors, 'r')])
            # with open(args.background_colors, 'r') as f:
            #     col = [str(s.strip()) for s in f]
        else: # it's a string in csv format
            col = [c.strip() for c in args.background_colors.split(',')]

        lst = {}
        i = 0

        if not background_colors:
            for c in background_clades:
                cc = c[:c.find('.')]

                if cc not in lst:
                    background_colors[c] = col[i % len(col)]
                    lst[cc] = col[i % len(col)]
                    i += 1
                else:
                    background_colors[c] = lst[cc]

    # get the levels that will use the internal annotation
    if args.annotations:
        annotations_list = [int(i.strip()) for i in args.annotations.strip().split(',')]

    # get the levels that will use the external legend annotation
    if args.external_annotations:
        external_annotations_list = [int(i.strip()) for i in args.external_annotations.strip().split(',')]

    # check overlapping between internal and external annotations
    if set(annotations_list) & set(external_annotations_list):
        print '[W] Some annotation levels are present in both internal and external params. The shared levels has been removed from the internal list.'
        annotations_list = list(set(annotations_list) - set(external_annotations_list))

    if args.lefse_input:
        # if the lefse_input is in biom format, convert it
        if get_file_type(args.lefse_input) in 'biom':
            try:
                biom = parse_biom(args.lefse_input, args.discard_otus, args.internal_levels)
                lefse_input = DataMatrix(StringIO(biom), args)
            except Exception as e:
                lin = True
                print 'Exception:', e
        else:
            if args.internal_levels:
                aaa = {}
                header = None
                with open(args.lefse_input, 'r') as f:
                    for r in f:
                        if header is None:
                            header = [s.strip() for s in r.split('\t')]
                        else:
                            row = r.split('\t')
                            aaa[row[0].strip().replace('|', '.')] = [float(s.strip()) for s in row[1:]]

                feats = add_missing_levels(aaa, summ=False)
                ss = '\t'.join(header) + '\n'
                ss += '\n'.join(['\t'.join([str(s) for s in [k] + feats[k]]) for k in feats])
                lefse_input = DataMatrix(StringIO(ss), args)
            else:
                lefse_input = DataMatrix(args.lefse_input, args)

        if not lin:
            taxa = [t.replace('|', '.').strip() for t in lefse_input.get_fnames()] # build taxonomy list

            # build all intermediate levels
            inter_lvls = []

            for t in taxa:
                s = ''

                for tt in t.split('.')[:-1]:
                    s = '.'.join([s, tt]) if s else tt

                    if (s not in taxa) and (s not in inter_lvls):
                        inter_lvls.append(s)

            taxa += inter_lvls
            taxa.sort()

            # check for duplicate taxa entries
            if len(taxa) != len(set(taxa)):
                print "There are duplicate taxa entries, please check the input file!"
                exit(1)

            # check if there are abundances to extract
            abundances = dict(lefse_input.get_averages())
            tot_abu = sum([abundances[a] for a in abundances if np.isfinite(abundances[a])])

            if tot_abu > 0:
                max_abundances = max([abundances[x] for x in abundances])
            else:
                abundances = dict()
                lin = False
                print "abundances: empty"
    else: # no lefse_input provided
        lin = True

    if args.lefse_output:
        # if the lefse_output is in biom format... I don't think it's possible!
        if get_file_type(args.lefse_output) in 'biom':
            lout = True
            print "Seriously?? LEfSe output file is not expected to be in biom format!"
        else:
            lst = []

            with open(args.lefse_output, 'r') as out_file:
                for line in out_file:
                    # print
                    # print '>>>'+line+'<<<'
                    # print
                    t, m, bk, es, pv = line.strip().split('\t')
                    lefse_output[t] = (es, bk, m, pv)

                    # get distinct biomarkers
                    if bk:
                        biomarkers |= set([bk])

                    # get all effect size
                    if es:
                        lst.append(float(es))

                max_effect_size = max(lst)

            # no lefse_input file provided!
            if (not taxa) and (not abundances): # build taxonomy list and abundaces map
                for t in lefse_output:
                    _, _, m, _ = lefse_output[t]
                    abundances[t.replace('.', '|')] = float(m)

                max_abundances = max([abundances[x] for x in abundances])

                for t in lefse_output:
                    scaled = scale_clade_size(args.min_clade_size, args.max_clade_size,
                                              abundances[t.replace('.', '|')], max_abundances)

                    if scaled >= args.abundance_threshold:
                        taxa.append(t.replace('|', '.').strip())
    elif not lin: # no lefse_output provided and lefse_input correctly red
        lout = True

        # find the xxx most abundant
        abundant = get_most_abundant(abundances, args.most_abundant)
        # print "abundant:", len(abundant), abundant

        # find the taxonomy level with at least yyy distinct childs from the xxx most abundant
        biomarkers = get_biomarkes(abundant, args.least_biomarkers)
        # print "biomarkers:", len(biomarkers), biomarkers

        # compose lefse_output variable
        for _, t in abundant:
            b = ''

            for bk in biomarkers:
                if bk in t:
                    b = bk

            lefse_output[t] = (2., b, '', '')

        max_effect_size = 2. # It's not gonna work... Maybe now??!?

    # no lefse_output and no lefse_input provided
    if lin and lout:
        print "You must provide at least one input file!"
        exit(1)

    # write the tree
    with open(args.tree, 'w') as tree_file:
        tree_file.write('\n'.join(taxa))

    # for each biomarker assign it to a different color
    if args.biomarkers2colors:
        if os.path.isfile(args.biomarkers2colors): # there exists a mapping file from biomarkers to colors read it
            with open(args.biomarkers2colors) as f:
                for row in f:
                    if not row.startswith('#'):
                        bk = row.strip().split('\t')[0]
                        cl = tuple([float(i.strip()) for i in row.strip().split('\t')[1].split(',')])
                        colors.append(cl)
                        color[bk] = colors.index(cl)
    else: # assign them automagically!
        i = 0

        for bk in biomarkers:
            color[bk] = i % len(colors)
            i += 1

    # print "color:", color

    # find max log abs value of effect size
    if lefse_output:
        lst = []

        for t in lefse_output:
            es, _, _, _ = lefse_output[t]

            if es:
                lst.append(abs(log10(float(es) / max_effect_size)))

        max_log_effect_size = max(lst)

    # write the annotation
    try:
        with open(args.annotation, 'w') as annot_file:
            # set the title
            if args.title:
                annot_file.write('\n'.join(['\t'.join(['title', args.title]),
                                            '\t'.join(['title_font_size', str(args.title_font_size)]), '\n']))

            # write some basic customizations
            annot_file.write('\n'.join(['\t'.join(['clade_separation', '0.5']),
                                        '\t'.join(['branch_bracket_depth', '0.8']),
                                        '\t'.join(['branch_bracket_width', '0.2']),
                                        '\t'.join(['annotation_legend_font_size', str(args.annotation_legend_font_size)]),
                                        '\t'.join(['class_legend_font_size', '10']),
                                        '\t'.join(['class_legend_marker_size', '1.5']), '\n']))

            # write the biomarkers' legend
            for bk in biomarkers:
                biom = pre_taxa.sub('', bk).replace('_', ' ').upper() # remove '{k|p|c|o|f|g|s|t}__'
                # print biom,
                rgb = scale_color(colors[color[bk]])
                # print rgb
                annot_file.write('\n'.join(['\t'.join([biom, 'annotation', biom]),
                                            '\t'.join([biom, 'clade_marker_color', rgb]),
                                            '\t'.join([biom, 'clade_marker_size', '40']), '\n']))

            # write the annotation for the tree
            for taxonomy in taxa:
                level = taxonomy.count('.') + 1 # which level is this taxonomy?
                clean_taxonomy = taxonomy[taxonomy.rfind('.') + 1:] # retrieve the last level in taxonomy
                cleanest_taxonomy = pre_taxa.sub('', clean_taxonomy).replace('_', ' ') # remove '{k|p|c|o|f|g|s|t}__' and substitute '_' with ' '
                scaled = args.def_clade_size

                # scaled the size of the clade by the average abundance
                if (taxonomy in abundances) or (taxonomy.replace('.', '|') in abundances):
                    try:
                        abu = abundances[taxonomy.replace('.', '|')]
                    except:
                        abu = abundances[taxonomy]

                    scaled = scale_clade_size(args.min_clade_size, args.max_clade_size, abu, max_abundances)

                annot_file.write(''.join(['\t'.join([clean_taxonomy, 'clade_marker_size', str(scaled)]), '\n']))

                # put a bakcground annotation to the levels specified by the user
                shaded_background = []

                for l in background_list:
                    if level >= l:
                        lst = [s.strip() for s in taxonomy.strip().split('.')]
                        t = '.'.join(lst[:l])

                        if t not in shaded_background:
                            shaded_background.append(t)

                            font_size = args.min_font_size + ((args.max_font_size - args.min_font_size) / l)

                            annot_file.write('\n'.join(['\t'.join([t, 'annotation_background_color', scale_color(colors[0])]),
                                                        '\t'.join([t, 'annotation', pre_taxa.sub('', t).replace('_', ' ')]), # remove '{k|p|c|o|f|g|s|t}__' and substitute '_' with ' '
                                                        '\t'.join([t, 'annotation_font_size', str(font_size)]), '\n']))

                # put a bakcground annotation to the clades specified by the user
                for c in background_colors:
                    bg_color = background_colors[c]

                    if not bg_color.startswith('#'):
                        bg_color = bg_color.replace('(', '').replace(')', '')
                        h, s, v = bg_color.split(';')
                        bg_color = scale_color((float(h.strip()) , float(s.strip()), float(v.strip())))

                    # check if the taxonomy has more than one level
                    lvls = [str(cc.strip()) for cc in c.split('.')]
                    done_clades = []

                    for l in lvls:
                        if (l in taxonomy) and (l not in done_clades):
                            lvl = taxonomy[:taxonomy.index(l)].count('.') + 1
                            font_size = args.min_font_size + ((args.max_font_size - args.min_font_size) / lvl)

                            annot_file.write('\n'.join(['\t'.join([l, 'annotation_background_color', bg_color]),
                                                        '\t'.join([l, 'annotation', pre_taxa.sub('', l).replace('_', ' ')]), # remove '{k|p|c|o|f|g|s|t}__' and substitute '_' with ' '
                                                        '\t'.join([l, 'annotation_font_size', str(font_size)]), '\n']))

                            done_clades.append(l)

                if lefse_output:
                    if taxonomy in lefse_output:
                        es, bk, _, _ = lefse_output[taxonomy]

                        # if it is a biomarker then color and label it!
                        if bk:
                            fac = log10(1. + 9. * (float(es) / max_effect_size))

                            try:
                                rgbs = scale_color(colors[color[bk]], fac)
                            except Exception as e:
                                print 'Exception:', e
                                print ' '.join(["[W] Assign to", taxonomy, "the default color:", colors[color[bk]]])
                                rgbs = colors[color[bk]]

                            annot_file.write(''.join(['\t'.join([clean_taxonomy, 'clade_marker_color', rgbs]), '\n']))

                            # write the annotation only if the abundance is above a given threshold and it is either internal or external annotation lists
                            if (scaled >= args.abundance_threshold) and \
                               ((level in annotations_list) or (level in external_annotations_list)):
                                font_size = args.min_font_size + ((args.max_font_size - args.min_font_size) / level)
                                annotation = cleanest_taxonomy if level in annotations_list else '*:' + cleanest_taxonomy

                                annot_file.write('\n'.join(['\t'.join([clean_taxonomy, 'annotation_background_color', rgbs]),
                                                            '\t'.join([clean_taxonomy, 'annotation', annotation]),
                                                            '\t'.join([clean_taxonomy, 'annotation_font_size', str(font_size)]), '\n']))
    except Exception as e:
        print 'Exception:', e


if __name__ == '__main__':
    main()

