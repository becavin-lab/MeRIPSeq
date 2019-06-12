#!/usr/bin/env python

import HTSeq
import collections
import sys, getopt
import numpy

#
#   Calculation of the Peak Over Median Score (POM score)
#
def calculate_htseq(annotation_file, input_file, output_file):
    '''
    Calculation of the Peak Over Median Score (POM score)
    '''
    # load wig file
    cvg = HTSeq.GenomicArray("auto", stranded=True, typecode='i')
    with open(input_file, 'rU') as wig_file:
        for row in wig_file:
            #print(row)
            chromo = row.split('\t')[0]
            start = int(row.split('\t')[1])
            nb_reads = int(float(row.split('\t')[2].strip()))
            interval = HTSeq.GenomicInterval(chromo, start, start + 1, "+")
            cvg[interval] = nb_reads

    # Count reads for all windows
    gtf_windows_file = HTSeq.GFF_Reader(annotation_file)
    windows = dict()
    # load GTF file with all transcripts
    k=1
    for feature in gtf_windows_file:
    #for feature in itertools.islice(gtf_windows_file, 100000):
        #print(feature)
        #if feature.iv.chrom == 'chr1':
        #if feature.name in windows.keys():
        #    print(feature.name)
        windows[feature.name] = feature.iv
        k += 1
        if k % 1000000 == 0:
            print('window ' + str(k)+'/3900000')

    print('Read window GTF '+str(len(windows))+' loaded')

    with open(output_file, "w") as htseqfile:
        for name, windowiv in windows.items():
            #Calculate median coverage for each window
            if len(name) != 0:
                #print(name,windowiv.length)
                values = []
                lengs = []
                for iv, value in cvg[windowiv].steps():
                    #print(iv.length)
                    #lengs.append(iv.length)
                    for i in range(iv.length):
                        values.append(value)
                medianvalue = numpy.median(numpy.array(values))
                #if medianvalue > 0:
                #print(windowiv, medianvalue, values)
                #print(windowiv, medianvalue, lengs)
                row = str(name)+'\t'+str(medianvalue)+'\n'
                htseqfile.write(row)


def main(argv):
    '''
    Main function of htseq_window.py

    :param argv:
    -p --path Path of the working folder
    -a --annotation Annotation file in GTF format
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv, "ha:i:o:", ["annotation=", "input=", "output="])
    except getopt.GetoptError:
        print('Cannot run command - Help: htseq_window.py -a <annotation> -i <input> -o <ouput>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('htseq_window.py -a <annotation> -i <input> -o <ouput>')
            sys.exit()
        elif opt in ("-a", "--annotation"):
            annotation_file = arg
        elif opt in ("-i", "--input"):
            input_file = arg
        elif opt in ("-o", "--output"):
            output_file = arg

    calculate_htseq(annotation_file, input_file, output_file)


if __name__ == "__main__":
    main(sys.argv[1:])
