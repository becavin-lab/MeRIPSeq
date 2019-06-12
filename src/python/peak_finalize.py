#!/usr/bin/env python
import csv
import sys, getopt
import os
import HTSeq
from collections import defaultdict

###  Setup all cut-off
# Fisher cut-off
PVALUE_CUTOFF = 0.001
# MACS2 cut-off
QVALUE_CUTOFF = 2
# POI cut-off
POM_CUTOFF = 4
POI_CUTOFF = 2
# RPMF cut-off
RPMF_CUTOFF = 10
# window occurence
OCCURENCE_CUTOFF = 2


def get_data_list(PATH, exp_design_name):
    file_name = PATH + '/ExpDesign/' + exp_design_name + '_exp_design.txt'
    bioconds_set = set()
    bioconds = list()
    with open(file_name, "rU") as exp_design_file:
        exp_design = csv.DictReader(exp_design_file, delimiter='\t')
        listdata_to_biocond = dict()
        for row in exp_design:
            bioconds_set.add(row.get('DataName'))
        for row in bioconds_set:
            bioconds.append(row)
        #print('bioconds',bioconds)
        return bioconds


def peaks_occurence(path, exp_design_name, peak_technique):
    print(exp_design_name)
    bioconds = get_data_list(path, exp_design_name)
    print(bioconds)
    print(len(bioconds))
    window_to_occurence = dict()
    biocond_to_occurence = dict()
    for dataset in bioconds:
        biocond_to_occurence[dataset] = 0

    PATH_PEAKS = path + '/PeakDetection/temp/'
    with open(PATH_PEAKS + peak_technique + '_' + bioconds[0] + '.txt', 'rU') as dataset_file:
        print(bioconds[0])
        dict_dataset = csv.DictReader(dataset_file, delimiter='\t')
        for row in dict_dataset:
            window_id = row['WindowId']
            window_to_occurence[window_id] = 0
    
    print(len(biocond_to_occurence))
    #print('nb windows ',len(window_to_occurence))

    headers = ['Windows_id']
    #dataset = listdata_to_biocond.keys()[0]
    for dataset in bioconds:
        print(dataset)
        with open(PATH_PEAKS + peak_technique + '_' + dataset + '.txt', 'rU') as dataset_file:
            dict_dataset = csv.DictReader(dataset_file, delimiter='\t')
            for row in dict_dataset:
                window_id = row['WindowId']
                if peak_technique == 'RPMF':
                    if row['RPMF'] != '':
                        rpmf = float(row['RPMF'])
                    else:
                        rpmf = 0
                    if rpmf >= RPMF_CUTOFF:
                        window_to_occurence[window_id] += 1
                        biocond_to_occurence[dataset] += 1
                elif peak_technique == 'POI':
                    pom = row['POM']
                    poi = row['POI']
                    if not (pom == '0' or pom == 'Nan'):
                        if not (poi == '0' or poi == 'Nan'):
                            pom = float(pom)
                            if pom > POM_CUTOFF:
                                poi = float(poi)
                                if poi > POI_CUTOFF:
                                    window_to_occurence[window_id] += 1
                                    biocond_to_occurence[dataset] += 1
                elif peak_technique == 'Fisher':
                    fisher = float(row['Correct_pvalues'])
                    if fisher < PVALUE_CUTOFF:
                        window_to_occurence[window_id] += 1
                        biocond_to_occurence[dataset] += 1

    print('Write file with ' + str(len(window_to_occurence)) + ' windows')
    PATH_PEAKS = path + '/PeakDetection/' + peak_technique + '/'
    with open(PATH_PEAKS + exp_design_name + '_' + peak_technique + '_occurence.txt', 'w') as peak_occurence_windows_file:
        if peak_technique == 'RPMF':
            header = 'RPMF>'+str(RPMF_CUTOFF)
        elif peak_technique == 'POI':
            header = 'POM>' + str(POM_CUTOFF) + 'POI>' + str(POI_CUTOFF)
        elif peak_technique == 'Fisher':
            header = 'pvalue<' + str(PVALUE_CUTOFF)
        headers = ['Window_id',header]
        peak_occurence_windows_file.write('\t'.join(headers)+'\n')
        for key, value in window_to_occurence.items():
            if value != 0:
                row = key+'\t'+str(value)+'\n'
                peak_occurence_windows_file.write(row)
    
    print('Write file with ' + str(len(biocond_to_occurence)) + ' windows')
    with open(PATH_PEAKS + exp_design_name + '_' + peak_technique + '_BioCond_occurence.txt', 'w') as peak_occurence_windows_file:
        if peak_technique == 'RPMF':
            header = 'RPMF>'+str(RPMF_CUTOFF)
        elif peak_technique == 'POI':
            header = 'POM>' + str(POM_CUTOFF) + 'POI>' + str(POI_CUTOFF)
        elif peak_technique == 'Fisher':
            header = 'pvalue<' + str(PVALUE_CUTOFF)
        headers = ['Window_id',header]
        peak_occurence_windows_file.write('\t'.join(headers)+'\n')
        for key, value in biocond_to_occurence.items():
            if value != 0:
                row = key+'\t'+str(value)+'\n'
                peak_occurence_windows_file.write(row)


def get_peak_position(path, exp_design_name, peak_technique, annotation_file):
    dict_peaks = dict()
    PATH_PEAKS = path + '/PeakDetection/' + peak_technique + '/'
    with open(PATH_PEAKS + exp_design_name + '_' + peak_technique + '_occurence.txt', 'rU') as occurence_file:
        occurence_file.readline()
        for row in occurence_file:
            occurence = int(row.split('\t')[1].strip())
            if occurence > OCCURENCE_CUTOFF:
                dict_peaks[row.split('\t')[0]] = 0
        #print('windows_peaks: ', len(dict_peaks))

    # annotate peak_window
    dict_peak_window = defaultdict(list)
    if peak_technique == 'POI':
        PATH_WINDOWS = path + '/Genome/' + annotation_file + '.exon.slidingwindows.gtf'
    else:
        PATH_WINDOWS = path + '/Genome/' + annotation_file + '.gene.slidingwindows.gtf'

    gtf_peaks_file = HTSeq.GFF_Reader(PATH_WINDOWS)
    k = 0
    for feature in gtf_peaks_file:
        window_id = feature.attr['ID']
        if window_id in dict_peaks:
            if window_id not in dict_peak_window:
                k += 1
                #if k%5000 == 0:
                #    print(k,len(dict_peaks),feature.iv)
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

    with open(PATH_PEAKS + exp_design_name + '_' + peak_technique + '_windows.bed', 'w') as annot_file:
        for key, value in dict_peak_window.items():
            #print(key,value)
            if key == 'None':
                print('none',value)
            annot_file.write('\t'.join(value) + '\n')


def merge_peaks(path, exp_design_name, peak_technique, bed_name):
    '''
    Need to install bedtools
    '''
    PATH_PEAKS = path + '/PeakDetection/' + peak_technique + '/'
    # sort
    bedfile = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_windows.bed'
    sort_bedfile = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_windows_sort.bed'
    merge_bedfile = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_windows_merge.bed'
    final_bedfile = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_' + bed_name + '.bed'

    sort_command = 'sort -k1,1 -k2,2n '+bedfile+' > '+sort_bedfile
    #print sort_command
    os.system(sort_command)
    #return_code = subprocess.call(sort_command)

    # merge
    merge_command = 'bedtools merge -c 4,5 -o distinct -i '+sort_bedfile+' > '+merge_bedfile
    #return_code = subprocess.call(sort_command)
    #print merge_command
    os.system(merge_command)
    # awk
    awk_command = 'awk \'{OFS = "\\t"; print $1,$2,$3,$4,"' + peak_technique + '","+"}\' '+merge_bedfile + ' > ' + final_bedfile
    os.system(awk_command)
    os.system('wc -l ' + bedfile)
    os.system('wc -l ' + sort_bedfile)
    os.system('wc -l ' + final_bedfile)
    os.remove(bedfile)
    os.remove(merge_bedfile)
    os.remove(sort_bedfile)


def filter_MACS_peaks(path, exp_design_name, bed_name):
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

    PATH_PEAKS = path + '/PeakDetection/temp/'
    list_data = get_data_list(path, exp_design_name)
    for dataset in list_data:
        with open(PATH_PEAKS + dataset + '_peaks.xls', 'rU') as dataset_file, \
                open(PATH_PEAKS + 'MACS2_' + dataset + '_clean.txt', 'w') as data_final_file:
            for row in dataset_file:
                if '#' not in row: # remove file informartion
                    if 'chr' in row: # keep only reference chromosome
                        if len(row) != 1:
                            data_final_file.write(row)
    peak_technique = "MACS2"
    
    PATH_PEAKS = path + '/PeakDetection/' + peak_technique + '/'
    bedfile_name = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_windows.bed'
    sort_bedfile_name = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_windows_sort.bed'
    merge_bedfile_name = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_windows_merge.bed'
    final_bedfile_name = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_' + bed_name + '.bed'

    biocond_to_occurence = dict()
    for dataset in list_data:
        biocond_to_occurence[dataset] = 0

    with open(bedfile_name, 'w') as final_bed:
        for dataset in list_data:
            with open(path + '/PeakDetection/temp/MACS2_' + dataset + '_clean.txt', 'rU') as dataset_file:
                    #open(path + '/PeakDetection/temp/MACS2_' + dataset + '.bed', 'w') as bed_file:
                csv_dataset = csv.DictReader(dataset_file, delimiter='\t')
                fields = ['chr','start','end','-log10(qvalue)']
                print(dataset)
                for row in csv_dataset:
                    value = float(row['-log10(qvalue)'])
                    if value > QVALUE_CUTOFF:
                        biocond_to_occurence[dataset] += 1
                        new_row = []
                        for field in fields:
                            new_row.append(row[field])
                        #bed_file.write('\t'.join(new_row)+'\n')
                        final_bed.write('\t'.join(new_row)+'\n')
            os.remove(path + '/PeakDetection/temp/MACS2_' + dataset + '_clean.txt')

    # merge all peaks together
    sort_command = 'sort -k1,1 -k2,2n ' + bedfile_name + ' > ' + sort_bedfile_name
    os.system(sort_command)
    merge_command = 'bedtools merge -i ' + sort_bedfile_name + ' -c 4 -o distinct > ' + merge_bedfile_name
    os.system(merge_command)

    print('Write file with ' + str(len(biocond_to_occurence)) + ' bioconds')
    with open(PATH_PEAKS + exp_design_name + '_' + peak_technique + '_BioCond_occurence.txt', 'w') as peak_occurence_windows_file:
        header = 'MACS2'
        headers = ['Window_id',header]
        peak_occurence_windows_file.write('\t'.join(headers)+'\n')
        for key, value in biocond_to_occurence.items():
            if value != 0:
                row = key+'\t'+str(value)+'\n'
                peak_occurence_windows_file.write(row)


    # filter by occurence
    with open(merge_bedfile_name, 'rU') as peaks_list, \
            open(PATH_PEAKS + exp_design_name + '_MACS2_occurence.txt', 'w') as occurence_file, \
            open(final_bedfile_name, 'w') as final_bed:
        occurence_file.write('Peak_id\tlogQ\n')
        index = 1
        for row in peaks_list:
            log_q = row.split('\t')[3].strip().split(',')
            occurence_file.write('Peak_'+str(index)+'\t'+str(len(log_q))+'\n')
            index += 1
            if len(log_q) > OCCURENCE_CUTOFF:
                new_row = [row.split('\t')[0],row.split('\t')[1],row.split('\t')[2],'.','MACS2','+']
                final_bed.write('\t'.join(new_row)+'\n')
    os.system('wc -l ' + bedfile_name)
    os.system('wc -l ' + sort_bedfile_name)
    os.system('wc -l ' + final_bedfile_name)
    os.remove(bedfile_name)
    os.remove(merge_bedfile_name)
    os.remove(sort_bedfile_name)


def main(argv):
    '''
    Main function of peak_calling_fisher.py
    :param argv:
    -p --path Path of the working folder
    -a --annotation Annotation file in GTF format
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:e:t:a:b:", ["path=", "expdesign=","peak_technique=","annotation=",'bed_name='])
    except getopt.GetoptError:
        print('Cannot run command - Help: peak_finalize.py -p <path> -e <expdesign> -t <peak_technique> -a <annotation> -b <bed_name>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('peak_finalize.py -p <path> -e <expdesign> -t <peak_technique> -a <annotation> -b <bed_name>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-e", "--expdesign"):
            exp_design_name = arg
        elif opt in ("-t", "--peaktechnique"):
            peak_technique = arg
        elif opt in ("-a", "--annotation"):
            annotation_file = arg
        elif opt in ("-b", "--bed_name"):
            bed_name = arg

    if peak_technique == 'MACS2':
        filter_MACS_peaks(path, exp_design_name, bed_name)
    else:
        peaks_occurence(path, exp_design_name, peak_technique)
        get_peak_position(path, exp_design_name, peak_technique, annotation_file)
        merge_peaks(path, exp_design_name, peak_technique, bed_name)


if __name__ == "__main__":
    main(sys.argv[1:])
