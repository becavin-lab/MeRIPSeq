from utils import m6a_utils
from utils import exp_design
from peaks import peak_utils
import csv
import sys
import os
from utils import sequence_python
import HTSeq
import pandas as pd
import math
import numpy
from collections import defaultdict


QVALUE_CUTOFF = 100


def clean_table(col_name, exp_design_name, PATH_PEAKS):
    """
    clean peak table

    INPUT:  PATH_PEAKS + 'MACS2/raw/' + dataset + '.xls'
    OUTPUT: PATH_PEAKS + exp_design_name + '_' + col_name + '_Filter_Biocond.txt'
            PATH_PEAKS + exp_design_name + '_windows_list.txt'

    :param col_name:
    :param exp_design_name:
    :param PATH_PEAKS:
    :return:
    """
    list_data = exp_design.get_data_list('m6aExpDesign_' + exp_design_name)
    for dataset in list_data:
        with open(PATH_PEAKS + col_name + '/raw/' + dataset + '_IPInput_peaks.xls', 'rU') as dataset_file, \
                open(PATH_PEAKS + col_name + '/' + dataset + '_IPInput_clean.txt', 'w') as data_final_file:
            for row in dataset_file:
                if '#' not in row:
                    if len(row) != 1:
                        data_final_file.write(row)


def filter_table(col_name, exp_design_name, PATH_PEAKS):
    """
    Filter peak table to remove peak which have -lof10(qvalue) > 20 for all biocond

    INPUT:  PATH_PEAKS + 'MACS2/raw/' + dataset + '.xls'
    OUTPUT: PATH_PEAKS + exp_design_name + '_' + col_name + '_Filter_Biocond.txt'
            PATH_PEAKS + exp_design_name + '_windows_list.txt'

    :param col_name:
    :param exp_design_name:
    :param PATH_PEAKS:
    :return:
    """
    list_data = exp_design.get_data_list('m6aExpDesign_' + exp_design_name)
    for dataset in list_data:
        with open(PATH_PEAKS + col_name + '/' + dataset + '_IPInput_clean.txt', 'rU') as dataset_file, \
                open(PATH_PEAKS + col_name + '/' + dataset + '_IPInput.txt', 'w') as data_final_file:
            csv_dataset = csv.DictReader(dataset_file, delimiter='\t')
            data_final_file.write('\t'.join(csv_dataset.fieldnames)+'\n')
            print(csv_dataset.fieldnames)
            for row in csv_dataset:
                value = float(row['-log10(qvalue)'])
                if value > QVALUE_CUTOFF:
                    new_row = []
                    for field in csv_dataset.fieldnames:
                        new_row.append(row[field])
                    data_final_file.write('\t'.join(new_row)+'\n')
        os.remove(PATH_PEAKS + col_name + '/' + dataset + '_IPInput_clean.txt')


def regroup_peaks(col_name, exp_design_name, PATH_PEAKS):
    list_data = exp_design.get_data_list('m6aExpDesign_' + exp_design_name)
    peak_position = dict()
    peak_qval = dict()
    peak_fold = dict()
    i = 1
    with open(PATH_PEAKS + col_name + '/' + list_data[0] + '_IPInput.txt', 'rU') as dataset_file:
        csv_dataset = csv.DictReader(dataset_file, delimiter='\t')
        for row in csv_dataset:
            begin = int(row['start'])
            end = int(row['end'])
            chromosome = row['chr']
            qval = float(row['-log10(qvalue)'])
            fold = float(row['fold_enrichment'])
            peak1 = sequence_python.SequencePython('peak' + str(i), begin, end, '+', 'mm10', chromosome)
            peak_id = 'Peak_'+str(i)
            peak_position[peak_id] = peak1
            dict_qval = dict()
            dict_qval[list_data[0]] = qval
            peak_qval[peak_id] = dict_qval
            dict_fold = dict()
            dict_fold[list_data[0]] = fold
            peak_fold[peak_id] = dict_fold
            i += 1
            print(list_data[0] + '\t' + peak_id)

    # search overlap between peak list
    for index in range(1, len(list_data)):
        k=0
        dataset = list_data[index]
        with open(PATH_PEAKS + col_name + '/' + dataset + '_IPInput.txt', 'rU') as dataset_file:
            csv_dataset = csv.DictReader(dataset_file, delimiter='\t')
            for row in csv_dataset:
                k += 1
                begin = int(row['start'])
                end = int(row['end'])
                chromosome = row['chr']
                qval = float(row['-log10(qvalue)'])
                fold = float(row['fold_enrichment'])
                peak2 = sequence_python.SequencePython('peak' + str(i), begin, end, '+', 'mm10', chromosome)
                found = False
                for peak_id, peak1 in peak_position.items():
                    # test if peaks overlap
                    if peak1.isoverlap(peak2):
                        # find unionPeak = the peak overlapping
                        new_peak = peak1.getOverlapPosition(peak2)
                        length_before = peak1.length
                        peak1.begin = new_peak[0]
                        peak1.end = new_peak[1]
                        peak1.beginIntersect = new_peak[2]
                        peak1.endIntersect = new_peak[3]
                        length_after = peak1.length
                        #print(length_before, (new_peak[1] - new_peak[0]), (peak1.end - peak1.begin))
                        peak_qval[peak_id][dataset] = qval
                        peak_fold[peak_id][dataset] = fold
                        found = True

                if not found:
                    peak_id = 'Peak_' + str(i)
                    print(dataset + '\t' + peak_id+'\t'+str(k))
                    peak_position[peak_id] = peak2
                    dict_qval = dict()
                    dict_qval[dataset] = qval
                    peak_qval[peak_id] = dict_qval
                    dict_fold = dict()
                    dict_fold[dataset] = fold
                    peak_fold[peak_id] = dict_fold
                    i += 1

    # save list of allPeaks to a file
    with open(PATH_PEAKS + col_name + '/' + exp_design_name + '_All.txt', 'w') as dataset_file:
        headers = ['Peak_id', 'chr', 'start', 'end', 'length']
        for dataset in list_data:
            headers.append(dataset)
        dataset_file.write('\t'.join(headers)+'\n')
        for peak_id, peak in peak_position.items():
            new_row = [peak_id, peak.chromosome, str(peak.begin), str(peak.end), str((peak.end - peak.begin))]
            if len(peak_fold[peak_id]) > 1:
                for dataset in list_data:
                    if dataset in peak_fold[peak_id]:
                        new_row.append(str(peak_fold[peak_id][dataset]))
                    else:
                        new_row.append(str(1))
                dataset_file.write('\t'.join(new_row)+'\n')


def complete_table(col_name, exp_design_name, PATH_PEAKS):
    with open(PATH_PEAKS + col_name + '/' + exp_design_name + '_All.txt', 'rU') as all_dataset_file, \
            open(PATH_PEAKS + col_name + '/' + exp_design_name + '_All_complete.txt', 'w') as final_all_dataset_file:
        list_data = exp_design.get_data_list('m6aExpDesign_' + exp_design_name)
        csv_all = csv.DictReader(all_dataset_file, delimiter='\t')
        headers = ['Peak_id', 'chr', 'start', 'end', 'length']
        for dataset in list_data:
            headers.append(dataset)
        final_all_dataset_file.write('\t'.join(headers) + '\n')
        ,
        for row in csv_all:
            begin = int(row['start'])
            end = int(row['end'])
            chromosome = row['chr']
            peak1 = sequence_python.SequencePython(row['Peak_id'], begin, end, '+', 'mm10', chromosome)
            for dataset in list_data:
                if row[dataset] == '1':
                    find = False
                    qval = 1
                    with open(PATH_PEAKS + col_name + '/' + dataset + '_IPInput_clean.txt', 'rU') as dataset_file:
                        csv_dataset = csv.DictReader(dataset_file, delimiter='\t')
                        for row_dataset in csv_dataset:
                            begin = int(row_dataset['start'])
                            end = int(row_dataset['end'])
                            chromosome = row_dataset['chr']
                            qval_temp = float(row_dataset['-log10(qvalue)'])
                            peak2 = sequence_python.SequencePython('peak', begin, end, '+', 'mm10', chromosome)
                            if peak1.isoverlap(peak2):
                                qval = qval_temp
                        row[dataset] = qval
            new_row = [peak1.name, peak1.chromosome, str(peak1.begin), str(peak1.end), str(peak1.length)]
            for dataset in list_data:
                new_row.append(str(row[dataset]))
            print(new_row)
            final_all_dataset_file.write('\t'.join(new_row) + '\n')


def annotate_table(exp_design_name, mouse_seq, ref_peaks_list, transcript_list, PATH_PEAKS, col_name):
    """
    Read list of windows_indexes and add :
        - annotation
        - motif presence score
        - relative position on transcript
        - overlapping refpeaks
        -

    INPUT: PATH_PEAKS + exp_design_name + '_windows_indexes.txt'
           PATH_PEAKS + exp_design_name + '_windows_list_Annot.txt'
    OUTPUT: PATH_PEAKS + exp_design_name + '_windows_indexes_Annot.txt'

    :param exp_design_name:
    :param mouse_seq:
    :param ref_peaks_list:
    :param PATH_PEAKS:
    :return:
    """

    with open(PATH_PEAKS + col_name + '/' + exp_design_name + '_All.txt', "rU") as table_all, \
        open(PATH_PEAKS + col_name + '/' + exp_design_name + '_' + col_name + '_Summary_All.txt', "w") as annot_index_file, \
        open(PATH_PEAKS + '/' + exp_design_name + '/' + exp_design_name + '_' + col_name + '_Summary.txt', "w") as annot_final_file:

        csv_table_all = csv.DictReader(table_all, delimiter = '\t')
        df_annot = pd.read_csv(m6a_utils.PATH_ANNOT + 'gencodeVM13/gencode.vM13.annotation.entrez.uniprot_clean.txt',index_col=0, sep='\t')
        headers = ['WindowId']
        list_biocond = exp_design.get_biocond_to_dataset('m6aExpDesign_' + exp_design_name)
        #for biocond in list_biocond:
        #    headers.append(biocond)
        list_data = exp_design.get_data_list('m6aExpDesign_' + exp_design_name)
        for dataset in list_data:
            headers.append(dataset)
        headers.extend(['WindowId','Motif','Relative_pos','Ref_Peaks','Nb_ref_peaks','Classification'])
        headers.extend(['chromo_window','begin_window','end_window','strand_window','type_transcript','index_window'])
        headers.append('Transcript_ID')
        headers.extend(df_annot.columns)
        annot_index_file.write('\t'.join(headers)+'\n')
        annot_final_file.write('\t'.join(headers) + '\n')
        for row in csv_table_all:
            window_id = row['Peak_id']
            peak = HTSeq.GenomicInterval(row['chr'], int(row['start']), int(row['end']), ".")
            transcript_id = ''
            for iv, value in transcript_list[peak].steps():
                if type(value) is HTSeq.GenomicFeature:
                    transcript_id = value.attr['transcript_id']
                    if transcript_id in df_annot.index:
                        transcript_annot = df_annot.loc[transcript_id]
                        chr = df_annot.loc[transcript_id]['chr']
                        begin_tr = int(df_annot.loc[transcript_id]['begin'])
                        end_tr = int(df_annot.loc[transcript_id]['end'])
                        strand_tr = df_annot.loc[transcript_id]['strand']
                        transcript_iv = HTSeq.GenomicInterval(chr, begin_tr, end_tr, strand_tr)
                    else:
                        transcript_id = ''

            # Calculate motif presence score
            sequence_score = 0
            sequence = mouse_seq[peak.chrom][peak.start:peak.end].seq
            for motif in m6a_utils.MOTIF_METH:
                if motif in sequence:
                    sequence_score += m6a_utils.MOTIF_METH[motif]

            # relative position on transcript
            relative_pos = 0.5
            if transcript_id != '':
                diff_start = float(peak.start + peak.length / 2) - transcript_iv.start
                if df_annot['strand'][transcript_id] == '-':
                    diff_start = float(peak.end + peak.length / 2) - transcript_iv.end
                    diff_start = - diff_start

                relative_pos = diff_start / transcript_iv.length

            # lengthTranscript = math.fabs(float(df_annot['begin'][values[0]]) - float(df_annot['end'][values[0]]))
            # diffStart = float(begin + length / 2) - df_annot['begin'][values[0]]
            # if strand == '-':
            #     diffStart = float(begin + length / 2) - df_annot['end'][values[0]]
            #     diffStart = - diffStart
            # relative_pos = diffStart / lengthTranscript

            # Look at overlapping refpeaks
            ref_peaks = ''
            for iv, value in ref_peaks_list[peak].steps():
                if len(value):
                    ref_peaks = ref_peaks + str(value) + ';;'
                    # print(ref_peaks)
            nb_ref_peaks = len(ref_peaks.split(';;')) - 1

            # apply classification
            classification = 0
            if (relative_pos < 0.3) or (relative_pos > 0.7):
                if sequence_score > 1:
                    # if nb_ref_peaks > 0:
                    classification = 1

            new_row = [window_id]
            #for biocond in list_biocond:
            #    new_row.append(row[biocond])
            for dataset in list_data:
                new_row.append(row[dataset])
            new_row.extend([window_id, str(sequence_score), str(relative_pos), ref_peaks, str(nb_ref_peaks), classification,
                   peak.chrom, str(peak.start), str(peak.end), str(peak.length)])
            if transcript_id != '':
                new_row.extend(['protein_coding', '1', transcript_id])
                for i in range(0, len(df_annot.columns)):
                    header = df_annot.columns[i]
                    if header == 'UniprotIDs':
                        uniprot = df_annot[header][transcript_id]
                        #print(df_annot[header][values[0]].isnull())
                        #print(type(uniprot))
                        if not uniprot == 'None' and not type(uniprot) == numpy.float:
                            uniprot = df_annot[header][transcript_id].split(';')[0]
                            new_row.append(uniprot)
                        else:
                            new_row.append('none')
                    else:
                        new_row.append(df_annot[header][transcript_id])
                annot_final_file.write('\t'.join([str(i) for i in new_row]) + '\n')
            annot_index_file.write('\t'.join([str(i) for i in new_row]) + '\n')


def read_ref_peaks_GTF():
    """
    General function which read ref peaks file and return them
    :return: ref_peaks_list[feature.iv] = feature.get_gff_line().strip().replace('\t',';')
    """
    gtf_peaks_file = HTSeq.GFF_Reader(m6a_utils.PATH + 'ReferencePeaks/All_Ref_Peaks.gtf')
    ref_peaks_list = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    # load GTF file with all transcripts
    for feature in gtf_peaks_file:
        ref_peaks_list[feature.iv] = feature.get_gff_line().strip().replace('\t',';')
    print('ref peaks GTF')
    return ref_peaks_list


def read_transcript_GTF():
    """
    General function which read ref peaks file and return them
    :return: ref_peaks_list[feature.iv] = feature.get_gff_line().strip().replace('\t',';')
    """
    gtf_peaks_file = HTSeq.GFF_Reader(m6a_utils.PATH_ANNOT + 'gencodeVM13/gencode.vM13.annotation.transcript.gtf')
    transcript_list = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    # load GTF file with all transcripts
    for feature in gtf_peaks_file:
        if feature.attr['gene_type'] == 'protein_coding':
            transcript_list[feature.iv] = feature
    print('protein coding transcript GTF')
    return transcript_list


def filter_motif_presence(col_name, exp_design_name, protein_coding, PATH_PEAKS):
    """
    Keep only motif with a motif present
    INPUT:  PATH_PEAKS + exp_design_name + '_' + col_name + '_Summary.txt'

    OUTPUT: PATH_PEAKS + exp_design_name + '_' + col_name + '_Summary_Motif.txt'

    :param col_name:
    :param exp_design_name:
    :param protein_coding:
    :param PATH_PEAKS:
    :return:
    """
    df_biocond = pd.read_csv(PATH_PEAKS + '/' + exp_design_name + '/' + exp_design_name + '_' + col_name + '_Summary.txt', index_col=0, sep='\t')
    list_peak = []
    print('Nb peaks: ',len(df_biocond))
    for row_name, row in df_biocond.iterrows():
        if row['Motif'] > 1:
            list_peak.append(row_name)

    df_filter = df_biocond.loc[list_peak]
    print('Nb peaks with motif: ', len(df_filter))
    df_filter.to_csv(PATH_PEAKS + '/' + exp_design_name + '/' + exp_design_name + '_' + col_name + '_Summary_Motif.txt', sep='\t')


def filter_per_gene(col_name, exp_design_name, protein_coding, PATH_PEAKS):
    """
    Keep only one peak with a motif per gene : keep the one with the highest POI value on all biocond

    INPUT:  PATH_PEAKS + exp_design_name + '_windows_indexes_motif.txt'
            PATH_PEAKS + exp_design_name + '_' + col_name + '_Summary_Motif.txt'

    OUTPUT: PATH_PEAKS + exp_design_name + '_peak_per_gene.txt' (for stats)
            PATH_PEAKS + exp_design_name + '_' + col_name + '_Summary_Gene.txt'

    :param col_name:
    :param exp_design_name:
    :param protein_coding:
    :param PATH_PEAKS:
    :return:
    """
    gene_group_dict = defaultdict(list)
    with open(PATH_PEAKS + '/' + exp_design_name + '/' + exp_design_name + '_' + col_name + '_Summary_Motif.txt', "rU") as table_file:
        gene_group_dict = defaultdict(list)
        csv_table = csv.DictReader(table_file, delimiter='\t')
        for row in csv_table:
            window_id = row['WindowId']
            gene_name = row['transcript_name']
            gene_group_dict[gene_name].append(window_id)

    with open(PATH_PEAKS + '/' + exp_design_name + '/' + exp_design_name + '_' + col_name + '_peak_per_gene.txt', 'w') as biocond_summary_file:
        for gene, list_peak in gene_group_dict.items():
            #print(len(list_peak))
            biocond_summary_file.write(str(len(list_peak))+'\n')


    df_biocond = pd.read_csv(PATH_PEAKS + '/' + exp_design_name + '/' + exp_design_name + '_' + col_name +'_Summary_Motif.txt', index_col=0, sep='\t')
    list_biocond = exp_design.get_biocond_to_dataset('m6aExpDesign_' + exp_design_name)
    list_peak = []
    for gene_name in gene_group_dict:
        list_window = gene_group_dict[gene_name]
        max_window = list_window[0]
        values = []
        for biocond in list_biocond:
            values.append(df_biocond[biocond][max_window])
        max_value = numpy.median(numpy.array(values))

        # find max value of all the window
        for i in range (1,len(list_window)):
            window = list_window[1]
            values = []
            for biocond in list_biocond:
                values.append(df_biocond[biocond][max_window])
            median_value = numpy.median(numpy.array(values))
            if median_value > max_value:
                max_value = median_value
                max_window = window

        list_peak.append(max_window)

    df_biocond = df_biocond.loc[list_peak]
    print('nb peaks: ', len(df_biocond))
    df_biocond.to_csv(PATH_PEAKS + '/' + exp_design_name + '/' + exp_design_name + '_' + col_name + '_Summary_Gene.txt', sep='\t')



PATH_PEAKS = m6a_utils.PATH_PEAKS + '/'

type_value = 'MACS2'
exp_design_name = 'Caecum'

# # remove peak in which p_value > 0.05
#clean_table(type_value, exp_design_name, PATH_PEAKS)
filter_table(type_value, exp_design_name, PATH_PEAKS)
# regroup all fold enrichment value
regroup_peaks(type_value, exp_design_name, PATH_PEAKS)
# add missing fold change value if existing
#complete_table(type_value, exp_design_name, PATH_PEAKS)

# calculate ConvMouse and ermFree median fold_enrichment (with Excel for the moment)

# annotate windows_indexes + add relat_pos, motif, ref_peak, annotation
#mouse_seq = m6a_utils.read_mouse_seq()
# #mouse_seq = ''
#ref_peaks_list = read_ref_peaks_GTF()
# #ref_peaks_list = ''
#transcript_list = read_transcript_GTF()
# # add annotation to the table
#annotate_table(exp_design_name, mouse_seq, ref_peaks_list, transcript_list, PATH_PEAKS, type_value)
# # keep only peak with a motif
#filter_motif_presence(type_value, exp_design_name, True, PATH_PEAKS)
# # keep only one peak per gene
# filter_per_gene(type_value, exp_design_name, True, PATH_PEAKS)

print('end')