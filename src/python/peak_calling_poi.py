#!/usr/bin/env python
import HTSeq
import csv
import sys, getopt
import numpy
import os
import re


def calculate_pom_score(PATH_PEAKS, PATH_RAW_COV, data, window_cutoff):
    """
    Calculation of the Peak Over Median Score (POM score)
    load : 'HTSeq_' + dataName + '_transcripts.txt'
    load all : 'HTSeq_' + dataName + '_windows_'+chromo+'.txt'
    save in : 'HTSeq_' + dataName+'.txt'
    :param dataName:
    :return:
    """
    # Read transcript result
    transcriptCovDict = {}
    with open(PATH_RAW_COV + 'HTSeq_Median_' + data + '_genes.txt', "r") as transcriptsfile:
        for row in transcriptsfile:
            transcriptCovDict[row.split('\t')[0]] = row.split('\t')[1].strip()
    with open(PATH_PEAKS + 'POI_' + data +'.txt',"w") as listPeaks:
        header = 'WindowId\tWindowcov\tTranscriptCov\tPOM'
        listPeaks.write(header+'\n')
        fileName = PATH_RAW_COV + 'HTSeq_Median_' + data + '_windows.txt'
        print 'read '+fileName
        with open(fileName) as windowCovFile:
            for row in windowCovFile:
                window_id = row.split('\t')[0]
                window_cov = float(row.split('\t')[1])
                transcript_id = window_id.split('_')[0]
                if transcriptCovDict.has_key(transcript_id):
                    transcript_cov = float(transcriptCovDict.get(transcript_id))
                    if window_cov > window_cutoff:
                        if transcript_cov > 0:
                            # POM Score calculation
                            pom = window_cov / transcript_cov
                            listPeaks.write(window_id+'\t'+str(window_cov)+'\t'+str(transcript_cov)+'\t'+str(pom)+'\n')
                        elif transcript_cov == 0:
                            # POM Score calculation
                            pom = window_cov / 0.25
                            listPeaks.write(window_id+'\t'+str(window_cov)+'\t'+str(transcript_cov)+'\t'+str(pom)+'\n')
                        else:
                            listPeaks.write(window_id + '\t' + str(window_cov) + '\t' + str(transcript_cov) + '\t' +
                                        'Nan' + '\n')
                    else:
                        listPeaks.write(window_id + '\t' + str(window_cov) + '\t' + str(transcript_cov) + '\t' +
                                        'Nan' + '\n')
                #print rowBegin,transcriptCov,fold
                else:
                    listPeaks.write(window_id + '\t' + str(window_cov) + '\tNan\tNan\n')
    print("POM Score calculated")


def calculate_poi_score(PATH_PEAKS, biocond, ip_data, input_data):
    """
    Calculation of the Peak Over Input Score (POI score)
    load : 'HTSeq_' + data_name + '_IP.txt'
    load : 'HTSeq_' + data_name + '_Input.txt'
    save in : 'HTSeq_' + data_name + '_IPInput.txt'
    :param dataName:
    :return:
    """
    with open(PATH_PEAKS + 'POI_' + ip_data +'.txt', "r") as ipfile, open(PATH_PEAKS + 'POI_' + input_data + '.txt', "r") as inputFile, open(PATH_PEAKS + 'POI_' + biocond + '.txt', "w") as wigfile:
        # Read transcript result
        header = ['WindowId','Windowcov','TranscriptCov','POM','Windowcov_Input','TranscriptCov_Input','POM_Input','POI']
        wigfile.write('\t'.join(header)+ '\n')
        ipfile.readline()   
        inputFile.readline()
        window_name_to_row = dict()
        for rowInput in inputFile:
            #print(rowInput)
            window_id_input = list(filter(None, re.split('\t| *', rowInput)))[0] #rowInput.split('\t')[0]
            window_name_to_row[window_id_input] = rowInput

        # WindowId 0
        # Windowcov 1
        # TranscriptCov 2
        # POM 3
        for row_ip in ipfile:
            window_id_ip = re.split('\t| *', row_ip)[0]  #row_ip.split('\t')[0]
            row_input = window_name_to_row[window_id_ip]
            row_input = row_input.replace(window_id_ip,'').strip()
            new_row = row_ip.strip() + '\t' + row_input.strip()

            # Calc ratio window
            # window_cov = float(row_ip.split('\t')[1].strip())
            # window_input_cov = float(row_input.split('\t')[0].strip())
            # if window_input_cov == 0:
            #     new_row += '\t'+str(window_cov)
            # else:
            #     ratio_windows = window_cov / window_input_cov
            #     new_row += '\t'+str(ratio_windows)

            # # Calc ratio transcript
            # transcript_cov = float(row_ip.split('\t')[-2].strip())
            # transcript_input_cov = float(row_input.split('\t')[-2].strip())
            # if transcript_input_cov == 0:
            #     new_row += '\t'+str(transcript_cov)
            # else:
            #     ratio_transcript = transcript_cov / transcript_input_cov
            #     new_row += '\t'+str(ratio_transcript)

            # Calc ratio 
            pom_ip = row_ip.split('\t')[-1].strip()
            pom_input = row_input.split('\t')[-1].strip()
            if (not pom_ip == 'Nan'):
                if (pom_input == 'Nan'):
                    pom_input = '0'
                pom_input = float(pom_input)
                if pom_input > 0:
                    poi = float(pom_ip) / float(pom_input)
                else:
                    poi = float(pom_ip) / 0.25
                new_row += '\t'+str(poi)
            else:
                new_row += '\t0'

            # write file
            wigfile.write(new_row + '\n')

    print("POI Score calculated")


def main(argv):
    '''
    Main function of peak_calling_poi.py
    :param argv:
    -p --path Path of the working folder
    -a --annotation Annotation file in GTF format
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:b:i:c:w:", ["path=", "biocond=","ip_data=", "input_data=", "window_cutoff="])
    except getopt.GetoptError:
        print('Cannot run command - Help: peak_calling_poi.py -p <path> -b <biocond> -i <ip_data> -c <input_data> -w <window_cutoff>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('peak_calling_poi.py -p <path> -b <biocond> -i <ip_data> -c <input_data> -w <window_cutoff>')
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

    calculate_pom_score(PATH_PEAKS, PATH_RAW_COV, ip_data, window_cutoff)
    calculate_pom_score(PATH_PEAKS, PATH_RAW_COV, input_data, window_cutoff)
    calculate_poi_score(PATH_PEAKS, biocond, ip_data, input_data)
    os.remove(PATH_PEAKS + 'POI_' + ip_data + '.txt')
    os.remove(PATH_PEAKS + 'POI_' + input_data + '.txt')


if __name__ == "__main__":
    main(sys.argv[1:])
