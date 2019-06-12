#!/usr/bin/env python
import HTSeq
import csv
import sys, getopt
import numpy
import os
import re


def calculate_rpm(PATH_PEAKS, PATH_RAW_COV, PATH_LIBRARY_SIZE, data):
    # Read ibrary size
    with open(PATH_LIBRARY_SIZE, "r") as file_library_size:
        library_size = 0.0
        for row in file_library_size:

            bam_name = re.split('\t| *', row)[0] #row.split('\t')[0]
            if bam_name == data:
                library_size = float(list(filter(None, re.split('\t| *', row)))[1].strip())/1000000 #float(row.split('\t')[1].strip()) / 1000000
                print('RPM FACTOR :', bam_name, library_size)
    
    with open(PATH_PEAKS + 'RPMF_' + data+'.txt', "w") as listPeaks:
        header = 'WindowId\tWindowcov\tRPM'
        listPeaks.write(header+'\n')
        fileName = PATH_RAW_COV + 'HTSeq_Count_' + data + '_windows.txt'
        print 'read '+ fileName
        with open(fileName) as windowCovFile:
            for row in windowCovFile:
                window_id = row.split('\t')[0]
                if len(re.split('\t| *', row)) == 1:
                    window_cov = 0
                else:
                    window_cov = float(re.split('\t| *', row)[1]) #float(row.split('\t')[1])
                # RPM calculation
                norm_cov = (window_cov / library_size)
                listPeaks.write(window_id+'\t'+str(window_cov)+'\t'+str(norm_cov)+'\n')
                #else:
                #    listPeaks.write(window_id + '\t' + str(window_cov) + '\tNan\n')
    print("RPM calculated")


def calculate_rpm_fold(PATH_PEAKS, biocond, ip_data, input_data, window_cutoff):
    with open(PATH_PEAKS + 'RPMF_' + ip_data + '.txt', "r") as ipfile, \
            open(PATH_PEAKS + 'RPMF_' + input_data + '.txt', "r") as inputFile, \
            open(PATH_PEAKS + 'RPMF_' + biocond + '.txt', "w") as bed_file:
        # Read transcript result
        header = ['WindowId','Windowcov','RPM','Windowcov_Input','RPM_Input','Ratio_windowcov','RPMF']
        bed_file.write('\t'.join(header)+ '\n')
        ipfile.readline()   
        inputFile.readline()
        window_name_to_row = dict()
        for rowInput in inputFile:
            #print(rowInput)
            window_id_input =  list(filter(None, re.split('\t| *', rowInput)))[0] #rowInput.split('\t')[0]
            window_name_to_row[window_id_input] = rowInput

        # WindowId 0
        # Windowcov 1
        # RPM 2
        for row_ip in ipfile:
            window_id_ip = list(filter(None, re.split('\t| *', row_ip)))[0] #row_ip.split('\t')[0]
            row_input = window_name_to_row[window_id_ip]
            row_input = row_input.replace(window_id_ip,'').strip()
            new_row = row_ip.strip() + '\t' + row_input.strip()

            # Calc ratio window
            #window_cov = float(row_ip.split('\t')[1].strip())
            window_cov = float(list(filter(None, re.split('\t| *', row_ip)))[1].strip()) 

            #window_input_cov = float(row_input.split('\t')[0].strip())
            window_input_cov = float(list(filter(None, re.split('\t| *', row_input)))[0].strip()) 

            if window_input_cov == 0:
                new_row += '\t'+str(window_cov)
            else:
                ratio_windows = window_cov / window_input_cov
                new_row += '\t'+str(ratio_windows)
            
            
            # Calc ratio RPM 
            if window_cov > window_cutoff:
                rpm_ip = re.split('\t| *', row_ip)[-1].strip() #row_ip.split('\t')[-1].strip()
                rpm_input = re.split('\t| *', row_input)[-1].strip() #row_input.split('\t')[-1].strip()
                rpm_ip = float(rpm_ip)
                rpm_input = float(rpm_input)
                if rpm_input > 0:
                    rpm = float(rpm_ip) / float(rpm_input)
                else:
                    rpm = float(rpm_ip) / 0.1
                new_row += '\t'+str(rpm)
            else:
                new_row += '\t0'
            
            # write file
            bed_file.write(new_row + '\n')
    print("RPMF Score calculated")


def main(argv):
    '''
    Main function of peak_calling_rpmf.py
    :param argv:
    -p --path Path of the working folder
    -a --annotation Annotation file in GTF format
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:b:i:c:w:", ["path=", "biocond=","ip_data=", "input_data=", "window_cutoff="])
    except getopt.GetoptError:
        print('Cannot run command - Help: peak_calling_rpmf.py -p <path> -b <biocond> -i <ip_data> -c <input_data> -w <window_cutoff>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('peak_calling_rpmf.py -p <path> -b <biocond> -i <ip_data> -c <input_data> -w <window_cutoff>')
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

    calculate_rpm(PATH_PEAKS, PATH_RAW_COV, PATH_LIBRARY_SIZE, ip_data)
    calculate_rpm(PATH_PEAKS, PATH_RAW_COV, PATH_LIBRARY_SIZE, input_data)
    calculate_rpm_fold(PATH_PEAKS, biocond, ip_data, input_data, window_cutoff)
    os.remove(PATH_PEAKS + 'RPMF_' + ip_data + '.txt')
    os.remove(PATH_PEAKS + 'RPMF_' + input_data + '.txt')


if __name__ == "__main__":
    main(sys.argv[1:])
