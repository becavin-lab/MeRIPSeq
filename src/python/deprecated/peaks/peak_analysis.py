import pytz
from utils import m6a_utils
from utils import exp_design
import subprocess
import csv
from collections import defaultdict
import collections
import numpy
import HTSeq
import pandas as pd
import math
import os
from utils import sequence_python

QVALUE_CUTOFF = 2
OCCURENCE_CUTOFF = 2

GENE_TYPE = ['protein_coding','processed_pseudogene','lincRNA','TEC','unprocessed_pseudogene']


def filter_MACS_peaks(exp_design_name, PATH_PEAKS):
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
    list_data = exp_design.get_biocond_list(exp_design_name)
    for dataset in list_data:
        with open(PATH_PEAKS + 'MACS2/raw/' + dataset + '_IPInput_peaks.xls', 'rU') as dataset_file, \
                open(PATH_PEAKS + 'MACS2/' + dataset + '_IPInput_clean.txt', 'w') as data_final_file:
            for row in dataset_file:
                if '#' not in row: # remove file informartion
                    if 'chr' in row: # keep only reference chromosome
                        if len(row) != 1:
                            data_final_file.write(row)
    with open(PATH_PEAKS + 'MACS2/All_' + exp_design_name + '_raw.bed', 'w') as final_bed:
        for dataset in list_data:
            with open(PATH_PEAKS + 'MACS2/' + dataset + '_IPInput_clean.txt', 'rU') as dataset_file, \
                    open(PATH_PEAKS + 'MACS2/' + dataset + '.bed', 'w') as bed_file:
                csv_dataset = csv.DictReader(dataset_file, delimiter='\t')
                fields = ['chr','start','end','-log10(qvalue)']
                print(dataset)
                for row in csv_dataset:
                    value = float(row['-log10(qvalue)'])
                    if value > QVALUE_CUTOFF:
                        new_row = []
                        for field in fields:
                            new_row.append(row[field])
                        bed_file.write('\t'.join(new_row)+'\n')
                        final_bed.write('\t'.join(new_row)+'\n')
            os.remove(PATH_PEAKS + 'MACS2/' + dataset + '_IPInput_clean.txt')

    # merge all peaks together
    sort_command = 'sort -k1,1 -k2,2n '+PATH_PEAKS+'MACS2/All_'+exp_design_name+'_raw.bed > '+PATH_PEAKS+'MACS2/All_'+exp_design_name+'_sort.bed'
    print sort_command
    merge_command = 'bedtools merge -i '+PATH_PEAKS+'MACS2/All_'+exp_design_name+'_sort.bed -c 4 -o distinct > '+PATH_PEAKS+'MACS2/All_'+exp_design_name+'.bed'
    print(merge_command)
    wc_command = 'wc -l '+PATH_PEAKS+'MACS2/*.bed'
    print(wc_command)

    # filter by occurence
    with open(PATH_PEAKS+'MACS2/All_'+exp_design_name+'.bed', 'rU') as peaks_list, \
            open(PATH_PEAKS + exp_design_name + '/Detection/' + exp_design_name + '_MACS2_occurence.txt', 'w') as occurence_file, \
            open(PATH_PEAKS + exp_design_name + '/Detection/' + exp_design_name + '_MACS2_windows_merge.bed', 'w') as final_bed:
        occurence_file.write('Peak_id\tlogQ\n')
        index = 1
        for row in peaks_list:
            log_q = row.split('\t')[3].strip().split(',')
            occurence_file.write('Peak_'+str(index)+'\t'+str(len(log_q))+'\n')
            index += 1
            if len(log_q) > QVALUE_CUTOFF:
                new_row = [row.split('\t')[0],row.split('\t')[1],row.split('\t')[2],'MACS2']
                final_bed.write('\t'.join(new_row)+'\n')
        print('MACS2 filtered')


def get_peak_position(exp_design_name, PATH_PEAKS, technique):
    dict_peaks = dict()
    with open(PATH_PEAKS + exp_design_name + '_'+technique+'_occurence.txt', 'rU') as occurence_file:
        occurence_file.readline()
        for row in occurence_file:
            occurence = int(row.split('\t')[1].strip())
            if occurence > OCCURENCE_CUTOFF:
                dict_peaks[row.split('\t')[0]] = 0
        print('windows_peaks: ', len(dict_peaks))

    # annotate peak_window
    dict_peak_window = defaultdict(list)
    if technique == 'POI':
        PATH_WINDOWS = m6a_utils.PATH + 'Genome/gencodeVM13/gencode.vM13.exon.slidingwindow.gtf'
    else:
        PATH_WINDOWS = m6a_utils.PATH + 'Genome/gencodeVM13/gencode.vM13.gene.slidingwindow.gtf'

    gtf_peaks_file = HTSeq.GFF_Reader(PATH_WINDOWS)
    k = 0
    for feature in gtf_peaks_file:
        window_id = feature.attr['ID']
        if window_id in dict_peaks:
            if window_id not in dict_peak_window:
                k += 1
                if k%1000 == 0:
                    print(k,len(dict_peaks),feature.iv)
                begin_window = feature.iv.start
                end_window = feature.iv.end
                chromo_window = feature.iv.chrom
                strand_window = feature.iv.strand
                dict_peak_window[window_id].append(chromo_window)
                dict_peak_window[window_id].append(str(begin_window))
                dict_peak_window[window_id].append(str(end_window))
                type_transcript = window_id.split('_window_')[1].rsplit('_', 1)[0]
                dict_peak_window[window_id].append(type_transcript)
                transcript_name = window_id.split('_window_')[0].rsplit('-', 1)[0]
                dict_peak_window[window_id].append(transcript_name)
                dict_peak_window[window_id].append(strand_window)

    with open(PATH_PEAKS +  exp_design_name + '_'+technique+'_windows.bed', 'w') as annot_file:
        for key, value in dict_peak_window.items():
            #print(key,value)
            if key == 'None':
                print('none',value)
            annot_file.write('\t'.join(value) + '\n')


def merge_peaks(exp_design_name, PATH_PEAKS, technique):
    '''
    Need to install bedtools
    '''
    # sort
    bedfile = PATH_PEAKS +  exp_design_name + '_'+technique+'_windows.bed'
    sort_bedfile = PATH_PEAKS +  exp_design_name + '_'+technique+'_windows_sort.bed'
    merge_bedfile = PATH_PEAKS + exp_design_name + '_' + technique + '_windows_merge.bed'

    sort_command = 'sort -k1,1 -k2,2n '+bedfile+' > '+sort_bedfile
    print sort_command
    #return_code = subprocess.call(sort_command)

    # merge
    merge_command = 'bedtools merge -s -c 4,5 -o distinct -i '+sort_bedfile+' > '+merge_bedfile
    #return_code = subprocess.call(sort_command)
    print merge_command
    print('wc -l ' + bedfile)
    print('wc -l ' + sort_bedfile)
    print('wc -l ' + merge_bedfile)


def annotate_peaks(exp_design_name, mouse_seq, gene_list, PATH_PEAKS, technique):
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
    if technique == '':
        bed_filename = PATH_PEAKS + '/' + exp_design_name+'_All_Techniques_Motif.bed'
        fasta_filename = PATH_PEAKS + '/' + exp_design_name + '.fasta'
        annot_filename = PATH_PEAKS + '/' + exp_design_name + '_Annot.txt'
        #bed_filename = PATH_PEAKS + '/' + exp_design_name + '_All_Techniques_Motif.bed'
        #fasta_filename = PATH_PEAKS + '/' + exp_design_name + '.fasta'
        #annot_filename = PATH_PEAKS + '/' + exp_design_name + '_Annot.txt'
    else:
        bed_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_windows_merge.bed'
        fasta_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '.fasta'
        annot_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_Annot.txt'


    with open(bed_filename, "rU") as bed_file, \
            open(fasta_filename, "w") as fasta_file, \
            open(annot_filename, "w") as annot_file:

        df_annot = pd.read_csv(m6a_utils.PATH_ANNOT + 'gencodeVM13/gencode.vM13.gene.annotation.txt',index_col=0, sep='\t')
        headers = ['WindowId','chromo_peak','begin_peak','end_peak','length_peak',
                   'Motif', 'Peak_Presence', 'Type_overlap', 'Ref_peak', 'Relat_pos','Gene_ID']
        headers.extend(df_annot.columns)
        annot_file.write('\t'.join(headers)+'\n')
        index_coding = 0
        index_nc = 0
        for row in bed_file:
            peak_technique = row.split('\t')[3].strip()
            peak = HTSeq.GenomicInterval(row.split('\t')[0], int(row.split('\t')[1]), int(row.split('\t')[2]), ".")
            gene_ids = []
            gene_ivs = []
            gene_types = []
            for gene_id in gene_list:
                gene = gene_list[gene_id]
                if peak.overlaps(gene.iv):
                    if technique == '':
                        if gene.attr['gene_type'] in GENE_TYPE:
                            # if gene.attr['gene_type'] == 'protein_coding' && gene.attr['gene_type'] in GENE_TYPE:
                            gene_ids.append(gene_id)
                            gene_ivs.append(gene.iv)
                            gene_types.append(gene.attr['gene_type'])
                    else:
                        gene_ids.append(gene_id)
                        gene_ivs.append(gene.iv)
                        gene_types.append(gene.attr['gene_type'])

            #print(len(transcript_ids),peak)
            if len(gene_ids) == 0:
                index_nc += 1
                window_id = 'Peak_NC_' + str(index_nc)
                if index_nc % 1000 == 0:
                    print(window_id)
                print('no ',peak)
            else:
                for i in range(len(gene_ids)):
                    gene_id = gene_ids[i]
                    gene_iv = gene_ivs[i]
                    gene_type = gene_types[i]
                    if gene_type == 'protein_coding':
                        index_coding += 1
                        window_id = 'Peak_' + str(index_coding)
                        if index_coding % 1000 == 0:
                            print(window_id)
                    else:
                        index_nc += 1
                        window_id = 'Peak_NC_' + str(index_nc)
                        if index_nc % 1000 == 0:
                            print(window_id)

                    # Calculate motif presence score
                    sequence_score = 0
                    sequence = mouse_seq[peak.chrom][peak.start:peak.end].seq
                    for motif in m6a_utils.MOTIF_METH:
                        if motif in sequence:
                            sequence_score += m6a_utils.MOTIF_METH[motif]

                    # motif search sequence
                    sequence_for_motif = mouse_seq[peak.chrom][(peak.start-100):(peak.end+100)].seq
                    fasta_file.write('>'+window_id+'\n'+str(sequence_for_motif)+'\n')

                    # relative position on transcript
                    relative_pos_start = 0.5
                    relative_pos_end = 0.5
                    if gene_id != '':
                        relative_pos_start = float(peak.start - gene_iv.start)
                        if gene_iv.strand == '-':
                            relative_pos_start = float(gene_iv.end - peak.end)
                        relative_pos_start = relative_pos_start / gene_iv.length

                        relative_pos_end = float(peak.end - gene_iv.start)
                        if gene_iv.strand == '-':
                            relative_pos_end = float(gene_iv.end - peak.start)
                        relative_pos_end = relative_pos_end / gene_iv.length
                    relat_pos = float(relative_pos_end - relative_pos_start)/2 + relative_pos_start

                    new_row = [window_id, peak.chrom, str(peak.start), str(peak.end), str(peak.length),
                               str(sequence_score), peak_technique, '_Over_', '_Ref_', str(relat_pos), gene_id]

                    if gene_id != '':
                        for i in range(0, len(df_annot.columns)):
                            header = df_annot.columns[i]
                            new_row.append(df_annot[header][gene_id])
                    else:
                        print('not found '+window_id)

                    annot_file.write('\t'.join([str(i) for i in new_row]) + '\n')
                    #annot_index_file.write('\t'.join([str(i) for i in new_row]) + '\n')


def peak_relative_pos(PATH_PEAKS, exp_design_name, technique):
    """
    Study the poisition in the transcript of peaks : UTR, Start_codon, CDS, exon, intron

    INPUT:  m6a_utils.PATH + 'Genome/gencodeVM13/gencode.vM13.transcript.description.txt'
            PATH_PEAKS + exp_design_name + '_' + col_name + '_Summary_Gene.txt'

    OUTPUT: PATH_PEAKS + exp_design_name + '_' + col_name + '_Relat_pos.txt'

    :param col_name:
    :param exp_design_name:
    :param PATH_PEAKS:
    :return:
    """
    if technique == '':
        relat_filename = PATH_PEAKS + '/' + exp_design_name + '_Relat_pos.txt'
        annot_filename = PATH_PEAKS + '/' + exp_design_name + '_Peaks_Ref.txt'
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_Peaks.txt'
    else:
        relat_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_Relat_pos.txt'
        annot_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_Annot.txt'
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_Peaks.txt'

    peak_to_overlap = dict()
    with open(m6a_utils.PATH + 'Genome/gencodeVM13/gencode.vM13.gene.description.txt', 'r') as transcript_description, \
            open(relat_filename, 'w') as relat_pos_file:
        transcript_to_details = defaultdict(list)
        headers = ['WindowId','Transcript_ID','Type','Pos']
        relat_pos_file.write('\t'.join(headers) + '\n')
        for row in transcript_description:
            transcript_id = row.split('\t')[0]
            transcript_to_details[transcript_id].append(row.strip().split('\t'))

        df_biocond = pd.read_csv(annot_filename, index_col=0, sep='\t')
        not_match = 0
        nb_start = 0
        nb_stop = 0
        nb_5utr = 0
        nb_3utr = 0
        nb_cds = 0
        nb_exon = 0
        nb_intron = 0
        for row_name, row in df_biocond.iterrows():
            if not '_NC_' in row_name:
                #print(row_name,row['Transcript_ID'],isinstance(row['Transcript_ID'],basestring))
                if (isinstance(row['Gene_ID'], basestring)) and (row['Gene_ID'] != ''):
                    transcript_id = row['Gene_ID']
                    relat_pos_start = int(row['begin_peak'])
                    relat_pos_end = int(row['end_peak'])
                    strand = row['strand']
                    chromo = row['chr']
                    if transcript_id in transcript_to_details:
                        all_elements = transcript_to_details[transcript_id]
                        match_element = []
                        start_codon_present = False
                        utr_present = False
                        cds_present = False
                        exon_present = False

                        for element in all_elements:
                            begin_relat = int(element[2])
                            end_relat = int(element[3])
                            element_gene = sequence_python.SequencePython('element', begin_relat, end_relat, strand, 'mm10', chromo)
                            peak = sequence_python.SequencePython('peak', relat_pos_start, relat_pos_end, strand, 'mm10', chromo)
                            if element_gene.isoverlap(peak):
                                #print(begin_relat,end_relat,relat_pos_start,relat_pos_end)

                                type_element = element[1]
                                relat_relat_pos = float(relat_pos_start-begin_relat / (end_relat-begin_relat))
                                element.append(str(relat_relat_pos))
                                match_element.append(element)
                                if type_element == '5UTR' or type_element == '3UTR':
                                    utr_present = True
                                #elif type_element == 'start_codon' or type_element == 'stop_codon':
                                #    start_codon_present = True
                                elif type_element == 'CDS':
                                    cds_present = True
                                elif type_element == 'exon':
                                    exon_present = True

                        # pos peak
                        pos_peak = 'intron'
                        relat_pos = relat_pos_start
                        intersection = [row_name, transcript_id, 'intron', str(relat_pos)]
                        if start_codon_present:
                            for element in match_element:
                                type_element = element[1]
                                if type_element == 'start_codon':
                                    intersection = [row_name, transcript_id, 'start_codon', str(relat_pos)]
                                    nb_start += 1
                                if type_element == 'stop_codon':
                                    intersection = [row_name, transcript_id, 'stop_codon', str(relat_pos),str(relat_pos_end)]
                                    nb_stop += 1
                        elif utr_present:
                            for element in match_element:
                                type_element = element[1]
                                if type_element == '5UTR':
                                    intersection = [row_name, transcript_id, '5UTR', str(relat_pos)]
                                    nb_5utr += 1
                                if type_element == '3UTR':
                                    intersection = [row_name, transcript_id, '3UTR', str(relat_pos)]
                                    nb_3utr += 1
                        elif cds_present:
                            for element in match_element:
                                type_element = element[1]
                                if type_element == 'CDS':
                                    intersection = [row_name, transcript_id, 'CDS', str(relat_pos)]
                                    nb_cds += 1
                        elif exon_present:
                            for element in match_element:
                                type_element = element[1]
                                if type_element == 'exon':
                                    intersection = [row_name, transcript_id, 'exon', str(relat_pos)]
                                    nb_exon += 1
                        else:
                            nb_intron += 1
                        #print(transcript_id,intersection)
                        #df_biocond['Type_overlap'][row_name] = intersection[2]
                        relat_pos_file.write('\t'.join(intersection)+'\n')
                    else:
                        #print('No: ', row_name)
                        not_match += 1
                        intersection = [row_name, '', 'no match', str(0)]
                        relat_pos_file.write('\t'.join(intersection) + '\n')
                else:
                    #print('No: ', row_name)
                    not_match += 1
                    intersection = [row_name, '', 'no match', str(0)]
                    relat_pos_file.write('\t'.join(intersection) + '\n')

                peak_to_overlap[intersection[0]] = intersection[2]

    # with open(peaks_filename, 'w') as peaks_file, \
    #     open(annot_filename, 'rU') as annot_file:
    #     peaks_file.write(annot_file.readline())
    #     for row in annot_file:
    #         peak_id = row.split('\t')[0]
    #         overlap = peak_to_overlap[peak_id]
    #         peaks_file.write(row.replace('_Over_',overlap))

    with open(relat_filename, 'w') as sum_overlap_file:
        new_row = 'Technique '+technique+'\n5\'UTR\t' + str(nb_5utr) + '\n3\'UTR\t' + str(nb_3utr) + \
            '\nCDS\t' + str(nb_cds) + '\nIntron \t' + str(nb_intron) + '\nNo Match \t' + str(not_match)
        print(str(new_row))
        sum_overlap_file.write(str(new_row))


def search_ref_peak(PATH_PEAKS, exp_design_name, technique):
    temp_PATH = m6a_utils.PATH + 'temp/'
    if technique == '':
        ref_filename = m6a_utils.PATH + 'ReferencePeaks/pc_ep_mouse_mm10.bed'
        clip_filename = m6a_utils.PATH + 'ReferencePeaks/sb_m6a_mouse_mm10.bed'
        trew_filename = m6a_utils.PATH + 'ReferencePeaks/trew_mouse_mm10.bed'

        final_filename = PATH_PEAKS + '/' + exp_design_name + '_Peaks_Ref.txt'
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_Peaks.txt'
    else:
        relat_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_Relat_pos.txt'
        annot_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_Annot.txt'
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_Peaks.txt'
    temp_PATH = m6a_utils.PATH+'temp/'
    # prepare bed
    chromo_to_peak = defaultdict(list)
    peak_to_pos = dict()
    peak_to_ref = defaultdict(list)
    peak_to_clip = defaultdict(list)
    peak_to_trew = defaultdict(list)
    with open(peaks_filename, 'rU') as peak_file:
        csv_peaks = csv.DictReader(peak_file, delimiter = '\t')
        for row in csv_peaks:
            ref_gene = sequence_python.SequencePython('peak', int(row['begin_peak']), int(row['end_peak']), '+', 'mm10',\
                                                   row['chromo_peak'])
            peak_to_pos[row['WindowId']] = ref_gene
            chromo_to_peak[row['chromo_peak']].append(row['WindowId'])
    print(len(chromo_to_peak))
    print('Go through Ref')

    with open(ref_filename, 'rU') as ref_file:
        ref_file.readline()
        k = 0
        for element in ref_file:
            #print(element)
            begin = int(element.split('\t')[1])
            end = int(element.split('\t')[2])
            chromo = element.split('\t')[0]
            info = element.split('\t')[4]
            #if chromo == 'chr2':
            #    break
            ref_peak = sequence_python.SequencePython('element', begin, end, '+', 'mm10', chromo)
            k+=1
            if k%10000==0:
                print(k)
            for peak_id in chromo_to_peak.get(chromo):
                peak = peak_to_pos[peak_id]
                if (peak.begin < ref_peak.begin):
                    overlapLength = peak.end - ref_peak.begin;
                else:
                    overlapLength = ref_peak.end - peak.begin
                if (overlapLength > 0):
                    #print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
                    peak_to_ref[peak_id].append(info)

    with open(temp_PATH+exp_design_name+'_Peak_to_ref.txt', 'w') as ref_peak_file:
        for peak_id, value in peak_to_ref.items():
            ref_peak_file.write(peak_id + '\t' + str(len(value))+'\t'+','.join(value)+'\n')

    with open(clip_filename, 'rU') as ref_file:
        ref_file.readline()
        k = 0
        for element in ref_file:
            #print(element)
            begin = int(element.split('\t')[1])
            end = int(element.split('\t')[2])
            chromo = element.split('\t')[0]
            #if chromo == 'chr2':
            #    break
            ref_peak = sequence_python.SequencePython('element', begin, end, '+', 'mm10', chromo)
            k+=1
            if k%10000==0:
                print(k)
            for peak_id in chromo_to_peak.get(chromo):
                peak = peak_to_pos[peak_id]
                if (peak.begin < ref_peak.begin):
                    overlapLength = peak.end - ref_peak.begin;
                else:
                    overlapLength = ref_peak.end - peak.begin
                if (overlapLength > 0):
                    #print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
                    peak_to_clip[peak_id].append(1)

    with open(temp_PATH+exp_design_name+'_Peak_to_clip.txt', 'w') as ref_peak_file:
        for peak_id, value in peak_to_clip.items():
            ref_peak_file.write(peak_id + '\t' + str(len(value))+'\n')

    print('trew')
    with open(trew_filename, 'rU') as ref_file:
        ref_file.readline()
        k = 0
        for element in ref_file:
            #print(element)
            begin = int(element.split('\t')[1])
            end = int(element.split('\t')[2])
            chromo = element.split('\t')[0]
            info = element.split('\t')[4]
            ref_peak = sequence_python.SequencePython('element', begin, end, '+', 'mm10', chromo)
            k+=1
            #if chromo == 'chr2':
            #    break
            if k%10000==0:
                print(k)
            for peak_id in chromo_to_peak.get(chromo):
                peak = peak_to_pos[peak_id]
                if (peak.begin < ref_peak.begin):
                    overlapLength = peak.end - ref_peak.begin;
                else:
                    overlapLength = ref_peak.end - peak.begin
                if (overlapLength > 0):
                    #print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
                    peak_to_trew[peak_id].append(info)

    with open(temp_PATH+exp_design_name+'_Peak_to_trew.txt', 'w') as ref_peak_file:
        for peak_id, value in peak_to_trew.items():
            ref_peak_file.write(peak_id + '\t' + str(len(value))+'\t'+','.join(value)+'\n')


def add_ref_peak(PATH_PEAKS, exp_design_name, technique):
    temp_PATH = m6a_utils.PATH + 'temp/'
    if technique == '':
        ref_filename = temp_PATH + '/' + exp_design_name+'_Peak_to_ref.txt'
        clip_filename = temp_PATH + '/' + exp_design_name+'_Peak_to_clip.txt'
        trew_filename = temp_PATH + '/' + exp_design_name+'_Peak_to_trew.txt'

        final_filename = PATH_PEAKS + '/' + exp_design_name + '_Peaks_Ref.txt'
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_Peaks.txt'
    else:
        relat_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_Relat_pos.txt'
        annot_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_Annot.txt'
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_Peaks.txt'

    # prepare dict
    peak_to_ref = dict()
    peak_to_clip = dict()
    peak_to_trew = dict()
    with open(ref_filename, 'rU') as ref_file:
        for row in ref_file:
            peak_to_ref[row.split('\t')[0]] = row.split('\t')[1].strip()
    with open(clip_filename, 'rU') as ref_file:
        for row in ref_file:
            peak_to_clip[row.split('\t')[0]] = row.split('\t')[1].strip()
    with open(trew_filename, 'rU') as ref_file:
        for row in ref_file:
            peak_to_trew[row.split('\t')[0]] = row.split('\t')[1].strip()
    print('Ref',len(peak_to_ref))
    print('Clip',len(peak_to_clip))
    print('TREW',len(peak_to_trew))

    peak_number = 0
    with open(peaks_filename, 'rU') as peak_file, \
          open(final_filename, 'w') as final_file:
        header = peak_file.readline().replace('Ref_peak','Ref_MeRIP\tCLIP\tTREW')
        final_file.write(header)
        for row in peak_file:
            peak_number += 1
            peak_id = row.split('\t')[0]
            ref=''
            trew=''
            clip=''
            if peak_id in peak_to_ref:
                ref = peak_to_ref[peak_id]
            if peak_id in peak_to_clip:
                clip = peak_to_clip[peak_id]
            if peak_id in peak_to_trew:
                trew = peak_to_trew[peak_id]

            new_row = row.replace('_Ref_',ref+'\t'+clip+'\t'+trew)
            final_file.write(new_row)
        print('Peak',peak_number,float(100 * len(peak_to_ref) / peak_number))


def find_utr(PATH_PEAKS, exp_design_name, technique):
    peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_Peaks_Ref.txt'
    utr_filename = PATH_PEAKS + '/' + exp_design_name + '_Peaks_UTR.txt'
    # Load UTR 5'
    utr5_gtf_filename = m6a_utils.PATH + 'Genome/gencodeVM13/gencode.vM13.annotation.5utr.gtf'
    utr5_list = dict()
    gene_to_utr5 = defaultdict(list)
    for feature in HTSeq.GFF_Reader(utr5_gtf_filename):
        # print(feature.attr['transcript_id'])
        utr5_list[feature.attr['utr_id']] = feature
        gene_to_utr5[feature.attr['gene_id']].append(feature.attr['utr_id'])
    # Load UTR 3'
    utr3_gtf_filename = m6a_utils.PATH + 'Genome/gencodeVM13/gencode.vM13.annotation.3utr.gtf'
    utr3_list = dict()
    gene_to_utr3 = defaultdict(list)
    for feature in HTSeq.GFF_Reader(utr3_gtf_filename):
        # print(feature.attr['transcript_id'])
        utr3_list[feature.attr['utr_id']] = feature
        gene_to_utr3[feature.attr['gene_id']].append(feature.attr['utr_id'])

    csv_peaks = csv.DictReader(open(peaks_filename, 'rU'), delimiter='\t')
    with open(utr_filename, 'w') as utr_file:
        utr_file.write('UTR_name\ttype\tpeak_id\n')
        for peak_row in csv_peaks:
            peak = sequence_python.SequencePython('peak', int(peak_row['begin_peak']), int(peak_row['end_peak']), '+', 'mm10', \
                                                  peak_row['chromo_peak'])
            gene_id = peak_row['Gene_ID']
            if gene_id != '':
                for utr5_name in gene_to_utr5[gene_id]:
                    utr_feature = utr5_list[utr5_name]
                    utr = sequence_python.SequencePython('utr', utr_feature.iv.start, utr_feature.iv.end,'+', 'mm10',utr_feature.iv.chrom)
                    if peak.isoverlap(utr):
                        print(utr5_name)
                        utr_file.write(utr5_name+'\tUTR5\t'+peak_row['WindowId']+'\n')
                for utr3_name in gene_to_utr3[gene_id]:
                    utr_feature = utr3_list[utr3_name]
                    utr = sequence_python.SequencePython('utr', utr_feature.iv.start, utr_feature.iv.end,'+', 'mm10',utr_feature.iv.chrom)
                    if peak.isoverlap(utr):
                        print(utr3_name)
                        utr_file.write(utr3_name+'\tUTR3\t'+peak_row['WindowId']+'\n')


def create_venn_diagram(PATH_PEAKS, exp_design_name, list_technique):
    '''
    Create files to do a venn diagram if technique overlaps
    :param PATH_PEAKS:
    :param exp_design_name:
    :param list_technique:
    :return:
    '''
    for technique in list_technique:
        with open(PATH_PEAKS+'/' + exp_design_name+'_All_Techniques_Motif.bed', 'rU') as annot_file, \
                open(PATH_PEAKS + '/' + exp_design_name + '_'+technique+'_Motif_List.txt', 'w') as list_peak:
            index = 1
            for row in annot_file:
                peak_id = 'Peak_' + str(index)
                peak_technique = row.split('\t')[3]
                print(peak_technique)
                if technique in peak_technique:
                    list_peak.write(peak_id+'\n')
                index += 1


def filter_motif_presence(exp_design_name, PATH_PEAKS, technique):
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
    annot_filename = PATH_PEAKS + exp_design_name + '_' + technique + '_Peaks.txt'
    bed_filename = PATH_PEAKS + exp_design_name + '_' + technique + '_Motif.bed'
    bed_nc_filename = PATH_PEAKS + exp_design_name + '_' + technique + '_Motif_NC.bed'
    motif_filename = PATH_PEAKS + exp_design_name + '_' + technique + '_Peaks_Motif.txt'
    motif_prot_filename = PATH_PEAKS + exp_design_name + '_' + technique + '_Peaks_Prot.txt'


    df_biocond = pd.read_csv(annot_filename, index_col=0, sep='\t')
    list_peak_motif = []
    list_peak_motif_prot_coding = []
    print('Nb peaks: ',len(df_biocond))
    with open(bed_nc_filename, 'w') as bed_file_nc, \
            open(bed_filename, 'w') as bed_file:
        for row_name, row in df_biocond.iterrows():
            if row['Motif'] > 1:
                bed_row = row['chromo_peak'] + '\t' + str(row['begin_peak']) + '\t' + str(
                    row['end_peak']) + '\t' + technique + '\t1\t+\n'
                if row['type'] in GENE_TYPE:
                    bed_file_nc.write(bed_row)
                    bed_file.write(bed_row)
                    list_peak_motif.append(row_name)
                if row['type'] == 'protein_coding':
                    list_peak_motif.append(row_name)
                    bed_file.write(bed_row)
                    list_peak_motif_prot_coding.append(row_name)

    df_filter = df_biocond.loc[list_peak_motif]
    print('Nb peaks with motif: ', len(df_filter))
    df_filter.to_csv(motif_filename, sep='\t')
    df_filter = df_biocond.loc[list_peak_motif_prot_coding]
    print('Nb peaks with motif and prot_coding: ', len(df_filter))
    df_filter.to_csv(motif_prot_filename, sep='\t')


def create_all_venn_diagram(PATH_PEAKS, list_technique):
    for technique in list_technique:
        with open(PATH_PEAKS + 'Caecum_Liver_Ref.bed', 'rU') as all_peak, \
                open(PATH_PEAKS + 'All_Peaks_'+technique+'_List.txt', 'w') as list_peak:
            index=0
            for row in all_peak:
                peak_id = 'Peak_'+str(index)
                index+=1
                if technique in row.split('\t')[3]:
                    list_peak.write(peak_id+'\n')


def read_gene_GTF():
    """
    General function which read ref peaks file and return them
    :return: ref_peaks_list[feature.iv] = feature.get_gff_line().strip().replace('\t',';')
    """
    gtf_peaks_file = HTSeq.GFF_Reader(m6a_utils.PATH_ANNOT + 'gencodeVM13/gencode.vM13.annotation.gene.gtf')
    gene_list = dict()
    # load GTF file with all transcripts
    for feature in gtf_peaks_file:
            #print(feature.attr['transcript_id'])
            gene_list[feature.attr['gene_id']] = feature
    print('Gene GTF - gencodeVM13/gencode.vM13.annotation.gene.gtf')
    print(len(gene_list))
    return gene_list



def create_peak_gtf_bed(PATH_PEAKS, exp_design_name):
    """
    Read all PATH_PEAKS+'/'+exp_design_name+'_'+technique+'_'+Final.txt
    Combine peaks
    and save to GFF
    :param list_technique:
    :return:
    """
    with open(PATH_PEAKS +'/'+exp_design_name+'_Final.gtf', 'w') as gtf_file, \
            open(PATH_PEAKS + '/' + exp_design_name + '_Final.bed', 'w') as bed_file, \
            open(PATH_PEAKS + '/' + exp_design_name + '_Peaks.txt', 'rU') as peak_file:
        csv_peak_file = csv.DictReader(peak_file, delimiter = '\t')
        for row in csv_peak_file:
            interval = HTSeq.GenomicInterval(row['chromo_peak'],int(row['begin_peak']),int(row['end_peak']), '.')
            peak_id = row['WindowId']
            feature = HTSeq.GenomicFeature(peak_id, 'exon', interval)
            print(feature.get_gff_line().strip() + '; gene_id \"'+peak_id+'\"')
            gtf_file.write(feature.get_gff_line().strip() + '; gene_id \"'+peak_id+'\"'+'\n')
            print(row['Peak_Presence'])
            bed_file.write(feature.iv.chrom+'\t'+str(feature.iv.start)+'\t'+str(feature.iv.end)+'\t'+row['Peak_Presence']
                           +'1\t+\t'+exp_design_name+'__'+peak_id+'\n')


def create_peak_fasta(PATH_PEAKS, exp_design_name):
    """
    Read all PATH_PEAKS+'/'+exp_design_name+'_'+technique+'_'+Final.txt
    Combine peaks
    and save to GFF
    :param list_technique:
    :return:
    """
    with open(PATH_PEAKS +'/'+exp_design_name+'.fasta', 'w') as fasta_file, \
            open(PATH_PEAKS + '/' + exp_design_name + '_Peaks.txt', 'rU') as peak_file:
        csv_peak_file = csv.DictReader(peak_file, delimiter = '\t')
        mouse_seq = m6a_utils.read_mouse_seq()
        for row in csv_peak_file:
            sequence_for_motif = mouse_seq[row['chromo_peak']][(int(row['begin_peak'])-100):(int(row['end_peak'])+100)].seq
            fasta_file.write('>' + row['WindowId'] + '\n' + str(sequence_for_motif) + '\n')


def merge_tissue_peaks():
    caecum_filename = m6a_utils.PATH_PEAKS + 'Caecum/Caecum_Peaks_Ref.txt'
    liver_filename = m6a_utils.PATH_PEAKS + 'Liver/Liver_Peaks_Ref.txt'
    ref_filename = m6a_utils.PATH_PEAKS + '../ReferencePeaks/Ref_Peaks.txt'
    final_filename = m6a_utils.PATH_PEAKS + 'Liver_Caecum_Peaks.txt'
    with open(final_filename, 'w') as final_file:
        df_caecum = pd.read_csv(caecum_filename, index_col=0, sep='\t')
        df_liver = pd.read_csv(liver_filename, index_col=0, sep='\t')
        df_ref = pd.read_csv(ref_filename, index_col=0, sep='\t')
        liver_peaks = dict()
        ref_peaks = dict()
        cecum_peaks = dict()
        peak_cecum_number = 0
        peak_ref_number = 0
        peak_liver_number = 0


        # Load peaks
        chrom_to_peaks_liver = defaultdict(list)
        chrom_to_peaks_ref = defaultdict(list)
        chrom_to_peaks_cecum = defaultdict(list)
        for peak_id_liver, row_liver in df_liver.iterrows():
            if not '_NC_' in peak_id_liver:
                peak_liver_number += 1
                liver_peak = sequence_python.SequencePython('Peak', row_liver['begin_peak'], row_liver['end_peak'], '+',
                                                            'mm10',row_liver['chr'])
                liver_peaks[peak_id_liver] = liver_peak
                chrom_to_peaks_liver[row_liver['chr']].append(peak_id_liver)
        for peak_id_ref, row_ref in df_ref.iterrows():
            if not '_NC_' in peak_id_ref:
                peak_ref_number += 1
                ref_peak = sequence_python.SequencePython('Peak', row_ref['begin_peak'], row_ref['end_peak'], '+',
                                                            'mm10',row_ref['chr'])
                ref_peaks[peak_id_ref] = ref_peak
                chrom_to_peaks_ref[row_ref['chr']].append(peak_id_ref)
        for peak_id_cecum, row_cecum in df_caecum.iterrows():
            if not '_NC_' in peak_id_cecum:
                peak_cecum_number += 1
                cecum_peak = sequence_python.SequencePython('Peak', row_cecum['begin_peak'], row_cecum['end_peak'], '+',
                                                            'mm10',row_cecum['chr'])
                cecum_peaks[peak_id_cecum] = cecum_peak
                chrom_to_peaks_cecum[row_cecum['chr']].append(peak_id_cecum)

        # create lists

        # liver_only = set()
        # ref_only() = set()
        # cecum_only() = set()
        #
        # liver_cecum = set()
        # liver_ref  = set()
        # cecum_ref = set()
        #
        # liver_cecum_ref = set()

        # cecum only
        nb_peak = 0
        index = 0

        chrom_to_2 = chrom_to_peaks_liver
        chrom_to_3 = chrom_to_peaks_ref
        dict_peaks_2 = liver_peaks
        dict_peaks_1 = cecum_peaks
        dict_peaks_3 = ref_peaks

        for peak_id, peak_pos in dict_peaks_1.items():
            index += 1
            if index%1000==0:
                print(index,peak_id)
            found_2 = False
            for peak_id_2 in chrom_to_2[peak_pos.chromosome]:
                peak_pos_2 = dict_peaks_2[peak_id_2]
                if (peak_pos_2.begin < peak_pos.begin):
                    overlapLength = peak_pos_2.end - peak_pos.begin;
                else:
                    overlapLength = peak_pos.end - peak_pos_2.begin
                if (overlapLength > 0):
                    # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end
                    found_2 = True
                    break
            found_3 = False
            for peak_id_3 in chrom_to_3[peak_pos.chromosome]:
                peak_pos_3 = dict_peaks_3[peak_id_3]
                if (peak_pos_3.begin < peak_pos.begin):
                    overlapLength = peak_pos_3.end - peak_pos.begin;
                else:
                    overlapLength = peak_pos.end - peak_pos_3.begin
                if (overlapLength > 0):
                    # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end
                    found_3 = True
                    break
            if (found_2) and (not found_3):
                nb_peak += 1

        print(nb_peak)
        print(len(dict_peaks_1))
        print(len(dict_peaks_2))
        print(len(dict_peaks_3))

        # overlap_liver = set()
        # overlap_cecum_liver = 0
        # overlap_liver_cecum = 0
        # overlap_ref_cecum = 0
        # overlap_ref_cecum_liver = 0
        # overlap_ref_liver = 0
        # final_file.write('Peak_Id\t'+'\t'.join(df_caecum.columns)+'\tLiver_peak\tNb_Liver_Peak\tgeneral_id\n')
        # index_general_id = 1
        # for peak_id_caecum, row in df_caecum.iterrows():
        #     #print(peak_id_caecum)
        #     if not '_NC_' in peak_id_caecum:
        #         caecum_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10', row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_liver[row['chr']]:
        #             liver_peak = liver_peaks[peak_id_liver]
        #             if (liver_peak.begin < caecum_peak.begin):
        #                 overlapLength = liver_peak.end - caecum_peak.begin;
        #             else:
        #                 overlapLength = caecum_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_liver)
        #                 overlap_liver.add(peak_id_liver)
        #
        #         print_row = ''
        #         # for colname in df_caecum.columns:
        #         #     print_row += str(row[colname])+'\t'
        #         # general_id = 'Peak_' + str(index_general_id)
        #         # new_row = 'Caecum_'+peak_id_caecum+'\t'+print_row+';'.join(peak_overlap)+'\t'+str(len(peak_overlap))+'\t'+general_id+'\n'
        #         if len(peak_overlap) != 0:
        #             overlap_liver_cecum += 1
        #             #print('no overlap')
        #         #final_file.write(new_row)
        #         index_general_id += 1
        #
        # for peak_id_caecum, row in df_liver.iterrows():
        #     #print(peak_id_caecum)
        #     if not '_NC_' in peak_id_caecum:
        #         caecum_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10', row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_cecum[row['chr']]:
        #             liver_peak = cecum_peaks[peak_id_liver]
        #             if (liver_peak.begin < caecum_peak.begin):
        #                 overlapLength = liver_peak.end - caecum_peak.begin;
        #             else:
        #                 overlapLength = caecum_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_liver)
        #                 overlap_liver.add(peak_id_liver)
        #
        #         print_row = ''
        #         # for colname in df_caecum.columns:
        #         #     print_row += str(row[colname])+'\t'
        #         # general_id = 'Peak_' + str(index_general_id)
        #         # new_row = 'Caecum_'+peak_id_caecum+'\t'+print_row+';'.join(peak_overlap)+'\t'+str(len(peak_overlap))+'\t'+general_id+'\n'
        #         if len(peak_overlap) != 0:
        #             overlap_cecum_liver += 1
        #             #print('no overlap')
        #         #final_file.write(new_row)
        #         index_general_id += 1
        #
        #
        # for peak_id_caecum, row in df_caecum.iterrows():
        #     #print(peak_id_caecum)
        #     if not '_NC_' in peak_id_caecum:
        #         caecum_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10', row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_ref[row['chr']]:
        #             liver_peak = ref_peaks[peak_id_liver]
        #             if (liver_peak.begin < caecum_peak.begin):
        #                 overlapLength = liver_peak.end - caecum_peak.begin;
        #             else:
        #                 overlapLength = caecum_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_liver)
        #                 overlap_liver.add(peak_id_liver)
        #
        #         if len(peak_overlap) != 0:
        #             overlap_ref_cecum += 1
        #
        #
        # for peak_id_caecum, row in df_liver.iterrows():
        #     #print(peak_id_caecum)
        #     if not '_NC_' in peak_id_caecum:
        #         caecum_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10', row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_ref[row['chr']]:
        #             liver_peak = ref_peaks[peak_id_liver]
        #             if (liver_peak.begin < caecum_peak.begin):
        #                 overlapLength = liver_peak.end - caecum_peak.begin;
        #             else:
        #                 overlapLength = caecum_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_liver)
        #                 overlap_liver.add(peak_id_liver)
        #
        #         if len(peak_overlap) != 0:
        #             overlap_ref_liver += 1
        #
        # for peak_id_caecum, row in df_liver.iterrows():
        #     #print(peak_id_caecum)
        #     if not '_NC_' in peak_id_caecum:
        #         caecum_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10', row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_ref[row['chr']]:
        #             liver_peak = ref_peaks[peak_id_liver]
        #             if (liver_peak.begin < caecum_peak.begin):
        #                 overlapLength = liver_peak.end - caecum_peak.begin;
        #             else:
        #                 overlapLength = caecum_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 for peak_id_liver in chrom_to_peaks_cecum[row['chr']]:
        #                     liver_peak = cecum_peaks[peak_id_liver]
        #                     if (liver_peak.begin < caecum_peak.begin):
        #                         overlapLength = liver_peak.end - caecum_peak.begin;
        #                     else:
        #                         overlapLength = caecum_peak.end - liver_peak.begin
        #                     if (overlapLength > 0):
        #                         # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                         peak_overlap.append(peak_id_liver)
        #                         overlap_liver.add(peak_id_liver)
        #
        #         if len(peak_overlap) != 0:
        #             overlap_ref_cecum_liver += 1


        # for peak_id_ref, row in df_ref.iterrows():
        #     print('Ref',peak_id_ref)
        #     if not '_NC_' in peak_id_ref:
        #         peak_ref_number += 1
        #         ref_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10',
        #                                                      row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_liver[row['chr']]:
        #             liver_peak = liver_peaks[peak_id_liver]
        #             if (liver_peak.begin < ref_peak.begin):
        #                 overlapLength = liver_peak.end - ref_peak.begin;
        #             else:
        #                 overlapLength = ref_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_liver)
        #         if len(peak_overlap) != 0:
        #             overlap_ref_liver += 1
        #             # print('no overlap')
        #
        #         peak_overlap = list()
        #         for peak_id_cecum in chrom_to_peaks_cecum[row['chr']]:
        #             cecum_peak = cecum_peaks[peak_id_cecum]
        #             if (cecum_peak.begin < ref_peak.begin):
        #                 overlapLength = cecum_peak.end - ref_peak.begin;
        #             else:
        #                 overlapLength = ref_peak.end - cecum_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_cecum)
        #         if len(peak_overlap) != 0:
        #             overlap_ref_cecum += 1
        #
        #         peak_overlap = list()
        #         for peak_id_cecum in chrom_to_peaks_cecum[row['chr']]:
        #             cecum_peak = cecum_peaks[peak_id_cecum]
        #             if (cecum_peak.begin < ref_peak.begin):
        #                 overlapLength = cecum_peak.end - ref_peak.begin;
        #             else:
        #                 overlapLength = ref_peak.end - cecum_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 for peak_id_liver in chrom_to_peaks_liver[row['chr']]:
        #                     liver_peak = liver_peaks[peak_id_liver]
        #                     if (liver_peak.begin < ref_peak.begin):
        #                         overlapLength = liver_peak.end - ref_peak.begin;
        #                     else:
        #                         overlapLength = ref_peak.end - liver_peak.begin
        #                     if (overlapLength > 0):
        #                         # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                         peak_overlap.append(peak_id_liver)
        #         if len(peak_overlap) != 0:
        #             overlap_ref_cecum_liver += 1
        # for peak_id_liver, row_liver in df_liver.iterrows():
        #     peak_liver_number += 1
        #     if peak_id_liver != overlap_liver:
        #         print_row = ''
        #         for colname in df_liver.columns:
        #             print_row += str(row_liver[colname]) + '\t'
        #         general_id = 'Peak_' + str(index_general_id)
        #         new_row = 'Liver_' + peak_id_liver + '\t' + print_row + '\t\t\t'+general_id+'\n'
        #         final_file.write(new_row)
        #         index_general_id += 1
        #
        # print(general_id,index_general_id)
        # print('Overlap Caecum + Liver', overlap_liver_cecum, peak_cecum_number, float(100*overlap_liver_cecum/peak_cecum_number))
        # print('Overlap Liver + Caecum', overlap_cecum_liver, peak_liver_number, float(100*len(overlap_liver)/peak_liver_number))
        # print('Overlap Ref + Liver', overlap_ref_liver, peak_liver_number, float(100 * overlap_ref_liver / peak_liver_number))
        # print('Overlap Ref + Cecum', overlap_ref_cecum, peak_cecum_number, float(100 * overlap_ref_cecum / peak_cecum_number))
        # print('Overlap Ref (Liver - Cecum)', overlap_ref_cecum_liver, overlap_ref_liver, overlap_ref_cecum, peak_ref_number, float(100 * overlap_ref_liver / peak_ref_number), float(100 * overlap_ref_cecum / peak_ref_number))


def get_merge_id_to_real_id(PATH_PEAKS):
    with open(PATH_PEAKS + '/Liver_to_generalID.txt', "w") as id_ref_file_liver, \
            open(PATH_PEAKS + '/Caecum_to_generalID.txt', "w") as id_ref_file_caecum, \
            open(PATH_PEAKS + '/Caecum_Liver_Final.bed', 'rU') as bed_file:
        id_ref_file_caecum.write('WindowId\tGeneral_id\n')
        id_ref_file_liver.write('WindowId\tGeneral_id\n')
        index = 1
        number_common = 1
        number_peak = 1
        for row in bed_file:
            peak_id = 'Peak_'+str(index)
            technique_peak_id = row.strip().split('\t')[3].split(',')
            print('tech',technique_peak_id)
            for peak in technique_peak_id:
                if 'Liver' in peak:
                    print(peak_id, 'Liver', peak)
                    peak = peak.replace('Liver__','')
                    id_ref_file_liver.write(peak+'\t'+peak_id+'\n')
                    number_peak+=1
                    if ('Caecum' in row):
                        number_common+=1
                else:
                    print(peak_id, 'Caecum', peak)
                    peak = peak.replace('Caecum__', '')
                    id_ref_file_caecum.write(peak + '\t' + peak_id + '\n')
                    number_peak += 1
                    if ('Liver' in row):
                        number_common+=1
            index+=1
        print('common peaks',number_common)
        print('number peaks',number_peak)



def filter_gene_counts(exp_design_name):
    '''
    Remove non protein coding genes in the HTseq count files
    :param exp_design_name:
    :return:
    '''
    list_coding_gene = dict()
    list_dataset = list()
    gene_list = read_gene_GTF()
    for gene_id, gene in gene_list.items():
        gene_type = gene.attr['gene_type']
        if gene_type == 'protein_coding':
            list_coding_gene[gene_id] = 0
        if gene_type in GENE_TYPE:
            list_coding_gene[gene_id] = 0
            #print(gene_id,gene.attr['gene_type'])

    print('len',len(list_coding_gene))
    with open(m6a_utils.PATH + 'ExpDesign/target_' + exp_design_name + '_Epi.txt', 'rU') as list_file:
        list_file.readline()
        for row in list_file:
            list_dataset.append(row.split('\t')[0].strip())
    list_yo = list()
    for dataset in list_dataset:
        # with open(m6a_utils.PATH + 'DiffExpression/HTSeq/Htseq_' + dataset + '.txt', 'rU') as list_file, \
        #         open(m6a_utils.PATH + 'DiffExpression/HTSeq/Htseq_' + dataset + '_all.txt', 'w') as list_file_all:
        #     for row in list_file:
        #         print(row)
        #         list_file_all.write(row)
        print(dataset)
        with open(m6a_utils.PATH + 'DiffExpression/HTSeq/Htseq_' + dataset + '_all.txt', 'rU') as list_file, \
            open(m6a_utils.PATH + 'DiffExpression/HTSeq/Htseq_' + dataset + '.txt', 'w') as list_file_all:
            for row in list_file:
                gene = row.split('\t')[0]
                #row = row.replace('ENSMU','\nENSMU')
                #row = row.replace('__', '\n__')
                if gene in list_coding_gene:
                    list_file_all.write(row)
                if '__' in gene:
                    list_file_all.write(row)

    print(len(list_coding_gene),len(list_dataset))


def modify_peak_counts(exp_design_name):
    '''
    Remove non protein coding genes in the HTseq count files
    :param exp_design_name:
    :return:
    '''
    list_dataset = list()
    with open(m6a_utils.PATH + 'ExpDesign/target_' + exp_design_name + '_Epi.txt', 'rU') as list_file:
        list_file.readline()
        for row in list_file:
            list_dataset.append(row.split('\t')[0].strip())
    list_yo = list()
    for dataset in list_dataset:
        # with open(m6a_utils.PATH + 'PeakDetection/'+exp_design_name+'/Peak_Count/Htseq_' + dataset + '_peak.txt', 'rU') as list_file, \
        #          open(m6a_utils.PATH + 'PeakDetection/'+exp_design_name+'/Peak_Count/Htseq_' + dataset + '_peak_all.txt', 'w') as list_file_all:
        #      for row in list_file:
        #          #print(row)
        #          list_file_all.write(row)
        with open(m6a_utils.PATH + 'PeakDetection/'+exp_design_name+'/Peak_Count/Htseq_' + dataset + '_peak_all.txt', 'rU') as list_file, \
                open(m6a_utils.PATH + 'PeakDetection/' + exp_design_name + '/Peak_Count/Htseq_' + dataset + '_peak_NC.txt','rU') as list_nc_file, \
                open(m6a_utils.PATH + 'PeakDetection/' + exp_design_name + '/Peak_Count/Htseq_' + dataset + '_peak.txt', 'w') as list_file_all:
            for row in list_file:
                gene = row.split('\t')[0]
                if '__' not in gene:
                    list_file_all.write(row)
            for row in list_nc_file:
                gene = row.split('\t')[0]
                if '__' not in gene:
                    list_file_all.write(row)


def parse_reads_distrib():
    with open(m6a_utils.PATH + 'rna_distrib.txt', 'w') as final_distrib_file:
        final_distrib_file.write('Data\tType\ttotal_reads\tTSS_up_10kb\t5UTR\tExon\tIntron\t3UTR\tTES_down_10kb\tIntergenic\n')
        for filename in os.listdir(m6a_utils.PATH+'RNADistrib/'):
            with open(m6a_utils.PATH+'RNADistrib/'+filename, 'rU') as distrib_file:
                data = filename.replace('_distrib.txt','')
                type = 'Input'
                if 'IP' in data:
                    type = 'IP'
                for row in distrib_file:
                    if 'Total Tags' in row:
                        total_tags = row.replace('Total Tags','').strip()
                    if 'Total Assigned Tags' in row:
                        total_assign = row.replace('Total Assigned Tags', '').strip()
                    if 'CDS_Exons' in row:
                        exon = get_number_tag(row)
                    if '5\'UTR_Exons' in row:
                        utr_5 = get_number_tag(row)
                    if '3\'UTR_Exons' in row:
                        utr_3 = get_number_tag(row)
                    if 'Introns' in row:
                        intron = get_number_tag(row)
                    if 'TSS_up_10kb' in row:
                        tss = get_number_tag(row)
                    if 'TES_down_10kb' in row:
                        tes = get_number_tag(row)
                unassigned = int(total_tags) - int(total_assign)
                new_row = [data,type,total_tags,tss,utr_5,exon,intron,utr_3,tes,str(unassigned)]
                final_distrib_file.write('\t'.join(new_row)+'\n')


def get_number_tag(row):
    temp_split = row.split(' ')
    temp2_split = list()
    for temp in temp_split:
        if temp != '':
            temp2_split.append(temp)
    return temp2_split[2]



    #print(len(list_coding_gene),len(list_dataset))



################################################################
################################################################
#   RUN ALL SCRIPTS HERE


#parse_reads_distrib()

#exp_design_name = 'Liver'
exp_design_name = 'Caecum'

print(exp_design_name)

# take care of MACS peaks
#PATH_PEAKS = m6a_utils.PATH_PEAKS + '/'
#filter_MACS_peaks(exp_design_name, PATH_PEAKS)

PATH_PEAKS = m6a_utils.PATH_PEAKS + exp_design_name + '/Detection/'
#os.mkdir(PATH_PEAKS)
# need : PATH_PEAKS/exp_design_name + '_' + col_name + '_Biocond.txt'

# Regroup windows for RPMF, POI and Fisher
#list_technique = ['RPMF', 'Fisher', 'POI']
#for technique in list_technique:
#    get_peak_position(exp_design_name, PATH_PEAKS, technique)
#for technique in list_technique:
#    merge_peaks(exp_design_name, PATH_PEAKS, technique)

# annotate peaks
# list_technique = ['RPMF','Fisher', 'POI','MACS2']
# mouse_seq = m6a_utils.read_mouse_seq()
# #mouse_seq = ''
# gene_list = read_gene_GTF()
# for technique in list_technique:
#     print(technique)
#     annotate_peaks(exp_design_name, mouse_seq, gene_list, PATH_PEAKS, technique)
# for technique in list_technique:
#     print(technique)
#     peak_relative_pos(PATH_PEAKS, exp_design_name, technique)
# for technique in list_technique:
#     print(technique)
#     filter_motif_presence(exp_design_name, PATH_PEAKS, technique)


#### Withouth motif cutoff
# list_technique = ['RPMF', 'Fisher', 'POI', 'MACS2']
# cat_command = 'cat '
# for technique in list_technique:
#     awk_command = 'awk \'{OFS = "\\t"; print $1,$2,$3,"+","'+technique+'"}\' '+PATH_PEAKS+'/'+ \
#                   exp_design_name+'_'+technique+'_windows_merge.bed > ' + PATH_PEAKS+'/'+ \
#                   exp_design_name+'_'+technique+'.bed'
#     print(awk_command)
#
# for technique in list_technique:
#     cat_command += PATH_PEAKS+'/' + exp_design_name+'_'+technique+'.bed '
# cat_command += ' > '+PATH_PEAKS+'/' + exp_design_name+'_All_Techniques_raw.bed'
# print(cat_command)
# sort_command = 'sort -k1,1 -k2,2n ' +PATH_PEAKS+'/' + exp_design_name+'_All_Techniques_raw.bed > ' + \
#                PATH_PEAKS + '/' + exp_design_name + '_All_Techniques_sort.bed'
# print sort_command
# merge_command = 'bedtools merge -i '+PATH_PEAKS+'/' + exp_design_name+\
#                 '_All_Techniques_sort.bed -c 4,5 -o distinct > '+PATH_PEAKS+'/' + exp_design_name+'_All_Techniques.bed'
# print(merge_command)
#
#create_venn_diagram(PATH_PEAKS, exp_design_name, list_technique)

#####  With motif cutoff
list_technique = ['RPMF', 'Fisher', 'POI', 'MACS2']
# cat_command = 'cat '
# for technique in list_technique:
#     awk_command = 'awk \'{OFS = "\\t"; print $1,$2,$3,"+","'+technique+'"}\' '+PATH_PEAKS+'/'+ \
#                   exp_design_name+'_'+technique+'_Motif.bed > ' + PATH_PEAKS+'/'+ \
#                   exp_design_name+'_'+technique+'.bed'
#     print(awk_command)
#
# for technique in list_technique:
#     cat_command += PATH_PEAKS+'/' + exp_design_name+'_'+technique+'_Motif.bed '
# cat_command += ' > '+PATH_PEAKS+'/' + exp_design_name+'_All_Techniques_Motif_raw.bed'
# print(cat_command)
# sort_command = 'sort -k1,1 -k2,2n ' +PATH_PEAKS+'/' + exp_design_name+'_All_Techniques_Motif_raw.bed > ' + \
#                PATH_PEAKS + '/' + exp_design_name + '_All_Techniques_Motif_sort.bed'
# print sort_command
# merge_command = 'bedtools merge -i '+PATH_PEAKS+'/' + exp_design_name+\
#                 '_All_Techniques_Motif_sort.bed -c 4,5 -o distinct > '+PATH_PEAKS+'/' + exp_design_name+'_All_Techniques_Motif.bed'
# print(merge_command)
#
#create_venn_diagram(PATH_PEAKS, exp_design_name, list_technique)

#####  For non_coding peaks
# list_technique = ['RPMF', 'Fisher', 'POI', 'MACS2']
# cat_command = 'cat '
# for technique in list_technique:
#     awk_command = 'awk \'{OFS = "\\t"; print $1,$2,$3,"+","'+technique+'"}\' '+PATH_PEAKS+'/'+ \
#                   exp_design_name+'_'+technique+'_Motif_NC.bed > ' + PATH_PEAKS+'/'+ \
#                   exp_design_name+'_'+technique+'_NC.bed'
#     print(awk_command)
#
# for technique in list_technique:
#     cat_command += PATH_PEAKS+'/' + exp_design_name+'_'+technique+'_Motif_NC.bed '
# cat_command += ' > '+PATH_PEAKS+'/' + exp_design_name+'_All_Techniques_Motif_raw_NC.bed'
# print(cat_command)
# sort_command = 'sort -k1,1 -k2,2n ' +PATH_PEAKS+'/' + exp_design_name+'_All_Techniques_Motif_raw_NC.bed > ' + \
#                PATH_PEAKS + '/' + exp_design_name + '_All_Techniques_Motif_sort_NC.bed'
# print sort_command
# merge_command = 'bedtools merge -i '+PATH_PEAKS+'/' + exp_design_name+\
#                 '_All_Techniques_Motif_sort_NC.bed -c 4,5 -o distinct > '+PATH_PEAKS+'/' + exp_design_name+'_All_Techniques_Motif_NC.bed'
# print(merge_command)

#create_venn_diagram(PATH_PEAKS, exp_design_name, list_technique)


PATH_PEAKS = m6a_utils.PATH_PEAKS + exp_design_name + '/'


# annotate peaks
#mouse_seq = m6a_utils.read_mouse_seq()
#mouse_seq = ''
#gene_list = read_gene_GTF()
#annotate_peaks(exp_design_name, mouse_seq, gene_list, PATH_PEAKS, '')
#peak_relative_pos(PATH_PEAKS, exp_design_name, '')

#search_ref_peak(PATH_PEAKS, exp_design_name, '')
#add_ref_peak(PATH_PEAKS, exp_design_name, '')
#find_utr(PATH_PEAKS, exp_design_name, '')
#create_peak_gtf_bed(PATH_PEAKS, exp_design_name)
#create_peak_fasta(PATH_PEAKS, exp_design_name)

#PATH_PEAKS = m6a_utils.PATH_PEAKS + exp_design_name

#filter_gene_counts(exp_design_name)
#modify_peak_counts(exp_design_name)

#get common list of peak between Caecum, Liver
# catCommand = 'cat '+m6a_utils.PATH_PEAKS+'Caecum//Caecum_Final.bed '+m6a_utils.PATH_PEAKS+'Liver/Liver_Final.bed > '+\
#              m6a_utils.PATH_PEAKS+'/Caecum_Liver_Final_raw.bed'
# print(catCommand)
sort_command = 'sort -k1,1 -k2,2n ' +m6a_utils.PATH_PEAKS+'/Caecum_Liver_Final_raw.bed > ' + \
                m6a_utils.PATH_PEAKS + '/Caecum_Liver_Final_sort.bed'
print sort_command
merge_command = 'bedtools merge -i '+m6a_utils.PATH_PEAKS+'/Caecum_Liver_Final_sort.bed -c 6 -o distinct > '+ \
                 m6a_utils.PATH_PEAKS + '/Caecum_Liver_Final.bed'
print(merge_command)

#awk_command = 'awk \'{OFS = "\\t"; print $1,$2,$3,"Ref","1","+","Ref"}\' '+m6a_utils.PATH + \
#              'ReferencePeaks/pc_ep_mouse_mm10_merge.bed > '+m6a_utils.PATH + 'ReferencePeaks/All_Peaks_List.bed'
#merge_tissue_peaks()
#get_merge_id_to_real_id(m6a_utils.PATH_PEAKS)




#### MOTIF search

#all_type = ['Meth','Meth+Gene+','Meth+Gene-','Meth-Gene-','Meth-Gene+','Meth_Epi']
#os.mkdir(m6a_utils.PATH_PEAKS + exp_design_name + '/Motif/')
# PATH_MOTIF = m6a_utils.PATH_PEAKS + exp_design_name + '/Motif/All_Peaks/'
# #os.mkdir(PATH_MOTIF)
# PATH_fasta = m6a_utils.PATH_PEAKS + exp_design_name + '/' + exp_design_name + '.fasta'
# meme_command = 'meme-chip -oc '+PATH_MOTIF+' -index-name meme-chip.html -time 300 -order 1 -db ' \
#             '/opt/meme/db/EUKARYOTE/jolma2013.meme -db /opt/meme/db/JASPAR/JASPAR_CORE_2014_vertebrates.meme -db ' \
#             '/opt/meme/db/MOUSE/uniprobe_mouse.meme -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 ' \
#             '-dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 ' + PATH_fasta
# print(meme_command)


# all_type = ['Meth','MethGene','MethNoGene','Meth_Epi']
# for type in all_type:
#     PATH_MOTIF = m6a_utils.PATH_PEAKS + exp_design_name + '/Motif/' + type +'/'
#     #os.mkdir(PATH_MOTIF)
#     PATH_fasta = m6a_utils.PATH_PEAKS + exp_design_name + '/DiffSummary/Sequence/' + exp_design_name + '_'+type+'.fasta'
#     meme_command = 'meme-chip -oc '+PATH_MOTIF+' -index-name meme-chip.html -time 300 -order 1 -db ' \
#             '/opt/meme/db/EUKARYOTE/jolma2013.meme -db /opt/meme/db/JASPAR/JASPAR_CORE_2014_vertebrates.meme -db ' \
#             '/opt/meme/db/MOUSE/uniprobe_mouse.meme -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 ' \
#             '-dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 ' + PATH_fasta
#     #meme-chip -oc BestPeaks/ -index-name meme-chip.html -time 300 -order 1 -db /opt/meme/db/EUKARYOTE/jolma2013.meme -db /opt/meme/db/JASPAR/JASPAR_CORE_2014_vertebrates.meme -db /opt/meme/db/MOUSE/uniprobe_mouse.meme -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 3 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 All_IPvsInput_Best.fasta
#     print(meme_command)

#os.mkdir(m6a_utils.PATH_PEAKS + exp_design_name + '/Motif/')
# PATH_MOTIF = m6a_utils.PATH_PEAKS + exp_design_name + '/Motif/5UTR/'
# #os.mkdir(PATH_MOTIF)
# PATH_fasta = m6a_utils.PATH + 'Genome/gencodeVM13/gencode.vM13.annotation.5utr.fasta'
# meme_command = 'meme-chip -oc '+PATH_MOTIF+' -index-name meme-chip.html -time 300 -order 1 -db ' \
#             '/opt/meme/db/EUKARYOTE/jolma2013.meme -db /opt/meme/db/JASPAR/JASPAR_CORE_2014_vertebrates.meme -db ' \
#             '/opt/meme/db/MOUSE/uniprobe_mouse.meme -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 ' \
#             '-dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 ' + PATH_fasta
# print(meme_command)
#
# #os.mkdir(m6a_utils.PATH_PEAKS + exp_design_name + '/Motif/')
# PATH_MOTIF = m6a_utils.PATH_PEAKS + exp_design_name + '/Motif/3UTR/'
# #os.mkdir(PATH_MOTIF)
# PATH_fasta = m6a_utils.PATH + 'Genome/gencodeVM13/gencode.vM13.annotation.3utr.fasta'
# meme_command = 'meme-chip -oc '+PATH_MOTIF+' -index-name meme-chip.html -time 300 -order 1 -db ' \
#             '/opt/meme/db/EUKARYOTE/jolma2013.meme -db /opt/meme/db/JASPAR/JASPAR_CORE_2014_vertebrates.meme -db ' \
#             '/opt/meme/db/MOUSE/uniprobe_mouse.meme -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 ' \
#             '-dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 ' + PATH_fasta
# print(meme_command)
