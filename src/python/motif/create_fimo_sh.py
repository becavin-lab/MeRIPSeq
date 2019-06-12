#!/usr/bin/env python
import sys, getopt


SOFTWARES = {'meme','dreme','centrimo','meme_tomtom'}
##TYPES_SEARCHES = {'All','Narrow','Const','Max'}


def run_fimo(path_analysis, exp_design_name, path_motif, path_motif_list, fimo_file):
    '''
    Creates a fimo.sh file with all the command to run
    :param path_analysis: path where the fasta files, background, and fimo results should be
    :param exp_design_name: name of the peaks : example: CecAm_Raw_1
    :param path_motif: path where the motif will be found
    :param motif_list_filename: path of the motif list
    :return: save file to path_motif + 'Fimo.sh'
    '''
    motif_list = list()
    with open(path_motif_list + '.txt', 'r') as motif_file:
        motif_file.readline()
        for row in motif_file:
            motif_list.append(row.split('\t')[0].strip())
    data_name = exp_design_name


    with open(fimo_file, 'w') as motif_sh:
        for motif in motif_list:
            motif_file = path_motif + 'motif/' + motif + '.meme'
            print(motif_file)
            fasta_file = path_analysis + 'fasta/' + data_name + '.fasta'
            background_file = path_analysis + 'results/' + data_name + '/background'
            folder_fimo = path_motif + 'fimo/' + motif + '_' + data_name
            fimo_sh = 'fimo --parse-genomic-coord --verbosity 1 --thresh 1e-2 --max-stored-scores 1000000 ' \
                      '--oc ' + folder_fimo + ' --bgfile ' + background_file + ' ' + motif_file + ' ' + fasta_file
            print(fimo_sh)
            motif_sh.write(fimo_sh+'\n')


def run_centrimo(path_analysis, exp_design_name, path_motif, motif_list_filename):
    '''
    Creates a centrimo.sh file with all the command to run
    :param path_analysis: path where the fasta files, background, and fimo results should be
    :param exp_design_name: name of the peaks : example: CecAm_Raw_1
    :param path_motif: path where the motif will be found
    :param motif_list_filename: path of the motif list
    :return: save file to path_motif + 'Centrimo.sh'
    '''
    motif_list = list()
    with open(path_motif + motif_list_filename + '.txt', 'r') as motif_file:
        motif_file.readline()
        for row in motif_file:
            motif_list.append(row.split('\t')[0].strip())
    data_name = exp_design_name

    with open(path_motif + 'Centrimo.sh', 'w') as motif_sh:
        for motif in motif_list:
            motif_file = path_analysis + 'motif/' + motif + '.meme'
            #print(motif_file)
            fasta_file = path_analysis + 'fasta/' + data_name + '.fasta'
            background_file = path_analysis + 'results/' + data_name + '/background'
            folder_centrimo = path_analysis + 'centrimo/' + motif + '_' + data_name
            centrimo_sh = 'centrimo -seqlen 150 -verbosity 1 -oc ' + folder_centrimo +\
                          ' --bgfile ' + background_file + ' -score 5.0 '\
                          '-ethresh 10.0 ' + fasta_file + ' '+motif_file
            print(centrimo_sh)
            motif_sh.write(centrimo_sh + '\n')



def main(argv):
    '''
    Main function of create_fimo_sh.py
    :param argv:
    -p --path Path of the working folder
    -e --expdesign Name of the ExpDesign
    -m --motiflist Name of the file with list of motif
    -f --fimo_file file to save for fimosearch
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:e:m:f:", ["path=", "expdesign=", "motiflist=", "fimo_file="])
    except getopt.GetoptError:
        print('Cannot run command - Help: create_fimo_sh.py -p <path> -e <expdesign> -m <motiflist> -f <fimofile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('create_fimo_sh.py -p <path> -e <expdesign> -m <motiflist> -f <fimofile>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-e", "--expdesign"):
            exp_design_name = arg
        elif opt in ("-m", "--motiflist"):
            motif_list = arg
        elif opt in ("-f", "--fimofile"):
            fimo_file = arg

    path_analysis = path + 'PeakDiffExpression/' + exp_design_name + '/Motif/'
    path_motif = path + 'PeakDiffExpression/Motif/'
    path_motif_list = path + 'PeakDiffExpression/Motif/' + motif_list
    # create fimo.sh for searching motif presence
    run_fimo(path_analysis, exp_design_name, path_motif, path_motif_list, fimo_file)


if __name__ == "__main__":
    main(sys.argv[1:])
