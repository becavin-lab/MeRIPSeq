#!/usr/bin/env python
import HTSeq
import csv
import sys, getopt
import numpy
import os
from fisher import pvalue
import pandas as pd
import re
from statsmodels.stats import multitest


def calculate_coverage(PATH_PEAKS, PATH_RAW_COV, data):
    with open(PATH_PEAKS + 'Fisher_' + data + '.txt',"w") as listPeaks:
        header = 'WindowId\tWindowcov'
        listPeaks.write(header+'\n')
        fileName = PATH_RAW_COV + 'HTSeq_Count_' + data + '_windows.txt'
        print 'read '+fileName
        with open(fileName) as windowCovFile:
            for row in windowCovFile:
                window_id = list(filter(None, re.split('\t| *', row)))[0] #row.split('\t')[0]
                if len(list(filter(None, re.split('\t| *', row)))) == 0: #len(row.split('\t')) == 0:
                    window_cov = 0
                else:
                    window_cov = float(list(filter(None, re.split('\t| *', row)))[1]) #float(row.split('\t')[1])
                listPeaks.write(window_id+'\t'+str(window_cov)+'\n')
    print("Cov calculated")


def get_library_size(PATH_LIBRARY_SIZE, ip_data, input_data):
    with open(PATH_LIBRARY_SIZE, "r") as file_library_size:
        library_size_ip = ''
        library_size_input = ''
        for row in file_library_size:
            bam_name = re.split('\t| *', row)[0] #row.split('\t')[0]
            if bam_name == ip_data:
                library_size_ip = float(list(filter(None, re.split('\t| *', row)))[1].strip()) #float(row.split('\t')[1].strip())
                print('IP Factor: ', bam_name, library_size_ip)
            if bam_name == input_data:
                library_size_input = float(list(filter(None, re.split('\t| *', row)))[1].strip())  #float(row.split('\t')[1].strip())
                print('Input Factor: ', bam_name, library_size_input)
        return [library_size_ip, library_size_input]


def calculate_fisher(PATH_PEAKS, biocond, ip_data, input_data, library_size_ip, library_size_input, window_cutoff):
    with open(PATH_PEAKS + 'Fisher_' + ip_data + '.txt', "r") as ipfile, \
            open(PATH_PEAKS + 'Fisher_' + input_data + '.txt', "r") as inputFile, \
            open(PATH_PEAKS + 'Fisher_' + biocond + '.txt', "w") as bed_file:
        # Read transcript result
        header = ['WindowId','Windowcov','Windowcov_Input','Ratio_windowcov','pvalue']
        bed_file.write('\t'.join(header)+ '\n')
        ipfile.readline()   
        inputFile.readline()
        window_name_to_row = dict()
        for rowInput in inputFile:
            #print(rowInput)
            window_id_input = re.split('\t| *', rowInput)[0]  #rowInput.split('\t')[0]
            window_name_to_row[window_id_input] = rowInput

        # WindowId 0
        # Windowcov 1
        # RPM 2
        index = 1
        for row_ip in ipfile:
            window_id_ip = re.split('\t| *', row_ip)[0]  #row_ip.split('\t')[0]
            row_input = window_name_to_row[window_id_ip]
            row_input = row_input.replace(window_id_ip,'').strip()
            new_row = row_ip.strip() + '\t' + row_input.strip()

            # Calc ratio window
            window_cov = float(list(filter(None, re.split('\t| *', row_ip)))[1].strip())  #float(row_ip.split('\t')[1].strip())
            window_input_cov = float(re.split('\t| *', row_input)[0].strip())  #float(row_input.split('\t')[0].strip())
            if window_input_cov == 0:
                new_row += '\t'+str(window_cov)
            else:
                ratio_windows = window_cov / window_input_cov
                new_row += '\t'+str(ratio_windows)

            # Calc fisher-test
            #print('library size input')
            #print(library_size_input)
            if window_cov > window_cutoff:
                p = pvalue(int(window_cov), library_size_ip, int(window_input_cov), library_size_input)
                new_row += '\t'+str(p.right_tail)
            else:
                new_row += '\t1'
            
            # write file
            bed_file.write(new_row + '\n')
            if index % 1000000 == 0:
                print('Windows '+str(index), '/25000000')
            index += 1
    print("Fisher test calculated")


def correct_p_values(PATH_PEAKS, biocond):
    print('Correct p-values')
    pvalues = pd.read_csv(PATH_PEAKS + 'Fisher_' + biocond + '.txt', sep = '\t')
    correct_pvalues = multitest.multipletests(pvalues['pvalue'], method='fdr_bh')
    pvalues['Correct_pvalues'] = correct_pvalues[1]
    pvalues.to_csv(PATH_PEAKS+'Fisher_' + biocond + '.txt', header=True, sep = '\t')
    print('P-value correction calculated')


def clean_table(PATH_PEAKS, biocond):
    with open(PATH_PEAKS+'Fisher_' + biocond + '.txt', 'r') as input, \
        open(PATH_PEAKS+'Fisher_' + biocond + '_clean.txt', 'w') as output:
        for row in input:
            newrow = row.split('\t',1) #re.split('\t| *', row)  #row.split('\t',1)
            output.write(newrow[1])
    with open(PATH_PEAKS + 'Fisher_' + biocond + '.txt', 'w') as ouput_new, \
        open(PATH_PEAKS + 'Fisher_' + biocond + '_clean.txt', 'r') as input_new:
        for row in input_new:
            ouput_new.write(row)
    os.remove(PATH_PEAKS + 'Fisher_' + biocond + '_clean.txt')


def main(argv):
    '''
    Main function of peak_calling_fisher.py
    :param argv:
    -p --path Path of the working folder
    -a --annotation Annotation file in GTF format
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:b:i:c:w:", ["path=", "biocond=","ip_data=", "input_data=", "window_cutoff="])
    except getopt.GetoptError:
        print('Cannot run command - Help: peak_calling_fisher.py -p <path> -b <biocond> -i <ip_data> -c <input_data> -w <window_cutoff>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('peak_calling_fisher.py -p <path> -b <biocond> -i <ip_data> -c <input_data> -w <window_cutoff>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-b", "--biocond"):
            biocond = arg
        elif opt in ("-i", "--ipdata"):
            ip_data = arg
        elif opt in ("-c", "--inputdata"):
            input_data = arg
        elif opt in ("-w", "--windcutoff"):
            window_cutoff = int(arg)

    PATH_RAW_COV = path + '/PeakDetection/RawCov/'
    PATH_PEAKS = path + '/PeakDetection/temp/'
    PATH_LIBRARY_SIZE = path + '/Seqdepth/STAR_nbReads.txt'

    calculate_coverage(PATH_PEAKS, PATH_RAW_COV, ip_data)
    calculate_coverage(PATH_PEAKS, PATH_RAW_COV, input_data)
    library_size = get_library_size(PATH_LIBRARY_SIZE, ip_data, input_data)
    calculate_fisher(PATH_PEAKS, biocond, ip_data, input_data, library_size[0], library_size[1], window_cutoff)
    correct_p_values(PATH_PEAKS, biocond)
    clean_table(PATH_PEAKS, biocond)
    os.remove(PATH_PEAKS + 'Fisher_' + ip_data + '.txt')
    os.remove(PATH_PEAKS + 'Fisher_' + input_data + '.txt')


if __name__ == "__main__":
    main(sys.argv[1:])
