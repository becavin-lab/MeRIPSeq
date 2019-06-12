#!/usr/bin/env python
import HTSeq
import csv
from collections import defaultdict
import sys, getopt

WINDOW_SIZE = 100

def create_sliding_gene_window_GTF(path_annot, annotation_file):
    '''
    Prepare window file in GTF for Peak detection
    Every overlapping window of 100bp on every gene is calculated
    :return:
    '''
    print('Create ' + annotation_file + '.gene.slidingwindows.gtf file for Peak detection')
    print('Every overlapping window of 100bp on every gene is calculated')
    overlap = WINDOW_SIZE/2
    gtf_file = HTSeq.GFF_Reader(path_annot + '/' + annotation_file + '.annotation.gene.gtf', end_included=True)
    windows = HTSeq.GenomicArrayOfSets( "auto", stranded=True)
    transcriptID = 1
    with open(path_annot + '/' + annotation_file + '.gene.slidingwindows.gtf', "w") as slidingGTF:
        for feature in gtf_file:
            interval=feature.iv
            transcriptID += 1
            if transcriptID%10000 == 0:
                print('Gene: '+str(transcriptID)+'/ 100000')
            windowID = 1
            if interval.strand == '+':
                begin = interval.start_d
                end = begin + WINDOW_SIZE
                while end < interval.end_d:
                    window = HTSeq.GenomicInterval(interval.chrom, begin, end+1, interval.strand)
                    featureWindow = HTSeq.GenomicFeature(feature.attr['gene_id']+'_window_' +
                                                       feature.attr['gene_type']+'_'+str(windowID),
                                                         'window', window)
                    windowID += 1
                    #print(featureWindow.get_gff_line())
                    begin += overlap
                    end = begin + WINDOW_SIZE
                    slidingGTF.write(featureWindow.get_gff_line())

                end = interval.end_d
                window = HTSeq.GenomicInterval(interval.chrom, begin, end, interval.strand)
                featureWindow = HTSeq.GenomicFeature(feature.attr['gene_id'] + '_window_' +
                                                    feature.attr['gene_type'] + '_' + str(windowID),
                                                     'window', window)
                slidingGTF.write(featureWindow.get_gff_line())
            else:
                #print(str(interval.end_d) + '  ' + str(interval.start_d))
                begin = interval.start_d
                end = begin - WINDOW_SIZE
                while end > interval.end_d:
                    window = HTSeq.GenomicInterval(interval.chrom, end, begin+1, interval.strand)
                    featureWindow = HTSeq.GenomicFeature(feature.attr['gene_id'] +'_window_'+
                                                         feature.attr['gene_type']+'_'+str(windowID),
                                                         'window', window)
                    windowID += 1
                    #print(featureWindow.get_gff_line())
                    begin -= overlap
                    end = begin - WINDOW_SIZE
                    slidingGTF.write(featureWindow.get_gff_line())

                end = interval.end_d
                window = HTSeq.GenomicInterval(interval.chrom, end+1, begin, interval.strand)
                featureWindow = HTSeq.GenomicFeature(feature.attr['gene_id'] + '_window_' +
                                                         feature.attr['gene_type'] + '_' + str(windowID),
                                                         'window', window)
                slidingGTF.write(featureWindow.get_gff_line())


def create_sliding_exon_window_GTF(path_annot, annotation_file):
    '''
        Prepare window file in GTF for Peak detection
        Every overlapping window of 100bp on every exon is calculated
        :return:
    '''
    print('Create ' + annotation_file + '.exon.slidingwindows.gtf file for Peak detection')
    print('Every overlapping window of 100bp on every exon is calculated')
    overlap = WINDOW_SIZE/2
    gtf_file = HTSeq.GFF_Reader(path_annot + '/' + annotation_file + '.annotation.exon.gtf', end_included=True)
    windows = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    transcriptID = 1
    with open(path_annot + '/' + annotation_file + '.exon.slidingwindows.gtf', "w") as slidingGTF:
        for feature in gtf_file:
            if feature.type == 'exon':
                interval = feature.iv
                transcriptID += 1
                if transcriptID%100000 == 0:
                    print('Gene: '+str(transcriptID)+'/ 100000')
                windowID = 1
                if interval.strand == '+':
                    begin = interval.start_d
                    end = begin + WINDOW_SIZE
                    while end < interval.end_d:
                        window = HTSeq.GenomicInterval(interval.chrom, begin, end+1, interval.strand)
                        featureWindow = HTSeq.GenomicFeature(feature.attr['transcript_name']+'_'+feature.attr['exon_number']
                                                           +'_window_'+
                                                           feature.attr['gene_type']+'_'+str(windowID),
                                                             'window', window)
                        windowID += 1
                        #print(featureWindow.get_gff_line())
                        begin += overlap
                        end = begin + WINDOW_SIZE
                        slidingGTF.write(featureWindow.get_gff_line())

                    end = interval.end_d
                    window = HTSeq.GenomicInterval(interval.chrom, begin, end, interval.strand)
                    featureWindow = HTSeq.GenomicFeature(feature.attr['transcript_name'] +'_'+feature.attr['exon_number']
                                                         + '_window_' +
                                                        feature.attr['gene_type'] + '_' + str(windowID),
                                                         'window', window)
                    if window.length > 1:
                        slidingGTF.write(featureWindow.get_gff_line())
                else:
                    #print(str(interval.end_d) + '  ' + str(interval.start_d))
                    begin = interval.start_d
                    end = begin - WINDOW_SIZE
                    while end > interval.end_d:
                        window = HTSeq.GenomicInterval(interval.chrom, end, begin+1, interval.strand)
                        featureWindow = HTSeq.GenomicFeature(feature.attr['transcript_name'] +'_'+feature.attr['exon_number']
                                                             +'_window_'+
                                                             feature.attr['gene_type']+'_'+str(windowID),
                                                             'window', window)
                        windowID += 1
                        #print(featureWindow.get_gff_line())
                        begin -= overlap
                        end = begin - WINDOW_SIZE
                        slidingGTF.write(featureWindow.get_gff_line())

                    end = interval.end_d
                    window = HTSeq.GenomicInterval(interval.chrom, end+1, begin, interval.strand)
                    featureWindow = HTSeq.GenomicFeature(feature.attr['transcript_name']+'_'+feature.attr['exon_number']+
                                                         '_window_' +
                                                             feature.attr['gene_type'] + '_' + str(windowID),
                                                             'window', window)
                    if window.length > 1:
                        slidingGTF.write(featureWindow.get_gff_line())



def main(argv):
    '''
    Main function of parse_annotation.py
    4 step in the workflow
    Step 1 - Parse GTF to extract gene and exon of protein_coding genes
    Step 2 -
    :param argv:
    -p --path Path of the working folder
    -a --annotation Annotation file in GTF format
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:a:", ["path=", "annotation="])
    except getopt.GetoptError:
        print('Cannot run command - Help: parse_annotation.py -p <path> -a <annotation_file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('parse_annotation.py -p <path> -a <annotation_file>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path_annot = arg
        elif opt in ("-a", "--annotation"):
            annotation_file = arg

    create_sliding_gene_window_GTF(path_annot, annotation_file)
    create_sliding_exon_window_GTF(path_annot, annotation_file)


if __name__ == "__main__":
    main(sys.argv[1:])

