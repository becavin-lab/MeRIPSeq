#!/usr/bin/env python
from collections import defaultdict
import os
import sys, getopt
from shutil import copyfile
import pandas as pd
import numpy as np
from Bio import SeqIO
import os.path


SOFTWARES = {'meme','dreme','centrimo','meme_tomtom'}
##TYPES_SEARCHES = {'All','Narrow','Const','Max'}


def figure_motif_clustering(path_motif, motif_list):
    '''
    Create SVG file with all motifs included
    :param path_motif:
    :param motif_list:
    :return:
    '''
    print('Create SVG file with all motifs')
    motiflist_filename = path_motif + motif_list
    df_summary = pd.read_csv(motiflist_filename + '.txt', sep='\t')
    svg_text = '<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<!-- Generator: Adobe Illustrator 15.0.0, SVG Export ' \
               'Plug-In . SVG Version: 6.00 Build 0)  -->\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"' \
               'http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n<svg version=\"1.1\" id=\"Calque_1\" ' \
               'xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" x=\"0px\" y=\"0px\"' \
               ' width=\"1600px\" height=\"'+str(len(df_summary.index)*100)+'px\" viewBox=\"0 0 1600 ' \
               + str(len(df_summary.index)*100)+'\" xml:space=\"preserve\">\n'
    i = 0
    for index, row in df_summary.iterrows():
        motif = row['Motif']
        data = row['Data']
        #expdesign = row['Exp_design']
        pvalue = str(row['P-value'])
        #comment = row['Comment']
        print(motif)
        image_file = path_motif + '../motif_figure/' + motif + '_nc.png'
        new_row = '<image overflow=\"visible\" width=\"299\" height=\"176\" xlink:href=\"'+image_file+'\"  ' \
                                                                                                      'transform=\"matrix(0.5184 0 0 0.5184 10.00 '+ str(float(100*i))+')\">\n</image>\n'
        image_file = path_motif + '../motif_figure/' + motif + '_rc.png'
        if os.path.isfile(image_file):
            new_row += '<image overflow=\"visible\" width=\"299\" height=\"176\" xlink:href=\"' + image_file + '\"  ' \
                                                                                                               'transform=\"matrix(0.5184 0 0 0.5184 170.00 '+ str(float(100*i))+')\">\n</image>\n'
        new_row += '<text transform=\"matrix(0.5184 0 0 0.5184 340 '+ str(float(50+100*i))+')\" font-family=' \
                                                                                           '\"\'MyriadPro-Regular\'\" font-size=\"40\">'+data+'</text>\n'
        #new_row += '<text transform=\"matrix(0.5184 0 0 0.5184 900 ' + str(float(50 + 100 * i)) + ')\" font-family=' \
        #                                                                                          '\"\'MyriadPro-Regular\'\" font-size=\"40\">' + expdesign + '</text>\n'
        new_row += '<text transform=\"matrix(0.5184 0 0 0.5184 1300 '+ str(float(50+100*i))+')\" font-family=' \
                                                                                            '\"\'MyriadPro-Regular\'\" font-size=\"40\">'+pvalue+'</text>\n'
        #new_row += '<text transform=\"matrix(0.5184 0 0 0.5184 1000 ' + str(float(50 + 100 * i)) + ')\" font-family=' \
        #            '\"\'MyriadPro-Regular\'\" font-size=\"40\">' + comment + '</text>\n'
        svg_text += new_row
        i += 1
    svg_text += '</svg>\n'

    with open(motiflist_filename + '.svg', 'w') as figure_file:
        figure_file.write(svg_text)
    print('Convert svg to png')
    script = 'convert ' + motiflist_filename + '.svg ' + motiflist_filename + '.png'
    os.system(script)
    os.remove(motiflist_filename + '.txt')


def main(argv):
    '''
    Main function of regroup_figures_motif.py
    :param argv:
    -p --path Path of the working folder
    -e --expdesign name of the exp_design
    -b --bed_name suffix of the bed file
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv, "hp:e:b:",
                                   ["path=", "expdesign=", 'bed_name='])
    except getopt.GetoptError:
        print(
            'Cannot run command - Help: regroup_figures_motif.py -p <path> -e <expdesign> -b <bed_name>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('peak_finalize.py -p <path> -e <expdesign> -b <bed_name>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-e", "--expdesign"):
            exp_design_name = arg
        elif opt in ("-b", "--bed_name"):
            bed_name = arg

    index_files = {'_log2count_list','_log2count_Pos_list', '_log2count_Comp_list', '_log2count_Biocond_list',
                  '_fimoscore_list', '_fimoscore_Pos_list', '_fimoscore_Comp_list', '_fimoscore_Biocond_list'}
    path_motif = path + 'PeakDiffExpression/Motif/' + 'Motif_' + exp_design_name + '/'
    # use clustered list of motif to recreate figure
    for index_file in index_files:
        motif_list_filename = 'Motif_' + exp_design_name + index_file
        figure_motif_clustering(path_motif, motif_list_filename)


if __name__ == "__main__":
    main(sys.argv[1:])
