#!/usr/bin/env python
import csv
import sys, getopt
import os


def create_mean_wig(path, exp_design_name, biocond, data='IP'):
    exp_design = csv.DictReader(open(path + '/ExpDesign/' + exp_design_name + '_exp_design.txt', 'r'), delimiter='\t')

    datasets = ''
    for row in exp_design:
        if biocond == 'all':
            datasets += path + '/Mapping/' + row[data+'_dataset'] + '.bw '
        else:
            if row['BioCond'] == biocond:
                datasets += path + '/Mapping/' + row[data+'_dataset'] + '.bw '

    if biocond == 'all':
        wig_filename = path + '/Mapping/' + exp_design_name + '_'+data+'_Mean.wig'
    else:
        wig_filename = path + '/Mapping/' + biocond + '_'+data+'_Mean.wig'
    commandline = 'wiggletools write ' + wig_filename + ' mean ' + datasets
    print(commandline)
    os.system(commandline)



def clean_wig(path, exp_design_name, data='IP'):
    wig_filename = path + '/Mapping/' + exp_design_name + '_'+data+'_Mean.wig'
    clean_wig_filename = path + '/Mapping/' + exp_design_name + '_'+data+'_Mean_clean.wig'
    k = 0
    with open(wig_filename, 'rU') as wig_file, \
            open(clean_wig_filename, 'w') as new_wig_file:
        for row in wig_file:
            k += 1
            if k % 1000000 == 0:
                print(k)

            if len(row.split('\t')) == 4:
                chromo = row.split('\t')[0]
                start = int(row.split('\t')[1])
                end = int(row.split('\t')[2])
                coverage = int(float(row.split('\t')[3].strip()) * 1000)
                for i in range(start, end+1):
                    new_wig_file.write(chromo+'\t'+str(i)+'\t'+str(coverage)+'\n')
            elif 'fixedStep' in row:
                chromo = row.split(' ')[1].replace('chrom=', '')
                start = int(row.split(' ')[2].replace('start=', ''))
            else:
                coverage = int(float(row.strip()) * 1000)
                new_wig_file.write(chromo + '\t' + str(start) + '\t' + str(coverage) + '\n')
                start += 1



def main(argv):
    '''
    Main function of peak_calling_fisher.py
    :param argv:
    -p --path Path of the working folder
    -a --annotation Annotation file in GTF format
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv, "hp:e:b:", ["path=", "expdesign=", "biocond="])
    except getopt.GetoptError:
        print('Cannot run command - Help: create_mean_wig.py -p <path> -e <expdesign> -b <biocond>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('create_mean_wig.py -p <path> -e <expdesign> -b <biocond>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-e", "--expdesign"):
            exp_design_name = arg
        elif opt in ("-b", "--biocond"):
            biocond = arg

    create_mean_wig(path, exp_design_name, biocond)
    if biocond == 'all':
        clean_wig(path, exp_design_name)
    else:
        clean_wig(path, biocond)


if __name__ == "__main__":
    main(sys.argv[1:])

