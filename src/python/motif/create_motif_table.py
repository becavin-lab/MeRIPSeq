#!/usr/bin/env python
import os
import pandas as pd
import sys, getopt
import numpy as np
from Bio import SeqIO
import os.path

fasta_cutoff = 50

def create_motif_table(path, exp_design_name, path_motif_list):
    '''
    Creates two tables summarizing occurence of motif for each peak
    :param path_analysis: path where the fasta files, background, and fimo results should be
    :param exp_design_name: name of the peaks : example: CecAm_Raw_1
    :param path_motif: path where the motif will be found
    :param motif_list_filename: path of the motif list
    :return: save file to path_analysis + motif_list_filename + '_score.txt'
                and path_analysis + motif_list_filename + '_count.txt'
    '''
    path_analysis = path + 'PeakDiffExpression/' + exp_design_name + '/Motif/'
    path_motif = path + 'PeakDiffExpression/Motif/'
    motif_list = list()
    with open(path_motif_list + '.txt', 'r') as motif_file:
        motif_file.readline()
        for row in motif_file:
            motif_list.append(row.split('\t')[0].strip())

    data_name = exp_design_name

    fasta_to_length = dict()
    print(path_analysis + 'fasta/' + exp_design_name + '.fasta' )
    fasta_seq = SeqIO.to_dict(SeqIO.parse(path_analysis + 'fasta/' + exp_design_name + '.fasta' , "fasta"))
    nb_peaks = len(fasta_seq.keys())
    np_motif = np.zeros((nb_peaks, len(motif_list)))
    df_motif_score = pd.DataFrame(np_motif, index=fasta_seq.keys(), columns=motif_list)
    np_motif = np.zeros((nb_peaks, len(motif_list)))
    df_motif_count = pd.DataFrame(np_motif, index=fasta_seq.keys(), columns=motif_list)
    np_motif = np.ones((nb_peaks, len(motif_list)))
    df_motif_pvalue = pd.DataFrame(np_motif, index=fasta_seq.keys(), columns=motif_list)

    for motif in motif_list:
        print(motif)
        df_motif_peaks = pd.read_csv(path_motif + 'fimo/' + motif + '_' + data_name +
                                     '/fimo.txt', index_col=0, sep='\t')
        grouped = df_motif_peaks[['sequence name','score']].groupby('sequence name')
        score = grouped.sum()
        df_motif_score[motif] = score
        counted = grouped.count()
        df_motif_count[motif] = counted
#           df_motif_pvalue[motif][peak] += pvalue
    df_motif_score.to_csv(path_motif_list + '_score.txt', sep='\t')
    df_motif_count.to_csv(path_motif_list + '_count.txt', sep='\t')

    #df_motif_pvalue.to_csv(path_motif + exp_design_name + '_' + motif_list_filename + '_pvalue.txt', sep='\t')



def get_size_fasta(path, exp_design_name):
    path_analysis = path + 'PeakDiffExpression/' + exp_design_name + '/Motif/'
    with open(path_analysis + 'Fasta_Summary.txt', 'w') as fasta_summary:
        for data_name in os.listdir(path_analysis+'fasta/'):
            if data_name.endswith('.fasta'):
                print(data_name)
                fasta_seq = SeqIO.to_dict(SeqIO.parse(path_analysis+'fasta/'+data_name, "fasta"))
                print(data_name,len(fasta_seq))
                fasta_summary.write(data_name.replace('.fasta','')+'\t'+str(len(fasta_seq))+'\n')


def motif_vs_fasta(path, exp_design_name, path_motif_list):
    '''
    Creates two tables summarizing occurence of motif for each fasta file
    :param path_analysis: path where the fasta files, background, and fimo results should be
    :param exp_design_name: name of the peaks : example: CecAm_Raw_1
    :param path_motif: path where the motif will be found
    :param motif_list_filename: path of the motif list
    :return: save file to path_analysis + motif_list_filename + '_score.txt'
                 and path_analysis + motif_list_filename + '_count.txt'
    '''
    # read fasta summary
    fasta_to_size = dict()
    path_analysis = path + 'PeakDiffExpression/' + exp_design_name + '/Motif/'
    with open(path_analysis + 'Fasta_Summary.txt', 'r') as fasta_summary:
        for row in fasta_summary:
            fasta_name = row.split('\t')[0]
            fasta_size = int(row.split('\t')[1].strip())
            if fasta_size > fasta_cutoff:
                fasta_to_size[fasta_name] = fasta_size
    # read list
    motif_list = list()
    with open(path_motif_list + '.txt', 'r') as motif_file:
        motif_file.readline()
        for row in motif_file:
            motif_list.append(row.split('\t')[0].strip())
    # read motif summarypath_motif_list + '_score.txt'
    score_filename = path_motif_list + '_score.txt'
    df_motif_score = pd.read_csv(score_filename, index_col=0, sep='\t')
    count_filename = path_motif_list + '_count.txt'
    df_motif_count = pd.read_csv(count_filename, index_col=0, sep='\t')

    # create motif vs fasta tables
    nb_fasta = len(fasta_to_size.keys())
    np_fasta = np.zeros((len(motif_list), nb_fasta))
    df_motif_fasta_score = pd.DataFrame(np_fasta, index=motif_list, columns=fasta_to_size.keys())
    np_fasta = np.zeros((len(motif_list), nb_fasta))
    df_motif_fasta_count = pd.DataFrame(np_fasta, index=motif_list, columns=fasta_to_size.keys())

    print(motif_list)
    # go through all fasta file
    for fasta_name, fasta_size in fasta_to_size.items():
        fasta_peaks = SeqIO.to_dict(SeqIO.parse(path_analysis+'fasta/'+fasta_name+'.fasta', 'fasta'))
        #print(fasta_peaks.keys())
        for motif in motif_list:
            motif_score = df_motif_score.loc[fasta_peaks.keys(), motif].sum()
            df_motif_fasta_score[fasta_name][motif] = float(motif_score)/float(fasta_size)
            motif_count = df_motif_count.loc[fasta_peaks.keys(), motif].sum()
            df_motif_fasta_count[fasta_name][motif] = float(motif_count)/float(fasta_size)

    df_motif_fasta_score.to_csv(path_motif_list + '_fasta_score.txt', sep='\t')
    df_motif_fasta_count.to_csv(path_motif_list + '_fasta_count.txt', sep='\t')

    # remove temp files
    os.remove(score_filename)
    os.remove(count_filename)



def main(argv):
    '''
    Main function of create_motif_table.py
    :param argv:
    -p --path Path of the working folder
    -e --expdesign Name of the ExpDesign
    -m --motiflist Name of the file with list of motif
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:e:m:", ["path=", "expdesign=", "motiflist="])
    except getopt.GetoptError:
        print('Cannot run command - Help: create_motif_table.py -p <path> -e <expdesign> -m <motiflist>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('find_peak_max.py -p <path> -e <expdesign> -m <motiflist>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-e", "--expdesign"):
            exp_design_name = arg
        elif opt in ("-m", "--motiflist"):
            motif_list = arg

    #exp_design_name = 'Cecum_all_MaxMaxValues_1'
    path_motif_list = path + 'PeakDiffExpression/Motif/' + motif_list

    # create a table grouping peak marginal score for each motif
    create_motif_table(path, exp_design_name, path_motif_list)

    # get size fasta
    get_size_fasta(path, exp_design_name)

    # Count marginal score and marginal motif presence for every fasta file
    motif_vs_fasta(path, exp_design_name, path_motif_list)


if __name__ == "__main__":
    main(sys.argv[1:])
