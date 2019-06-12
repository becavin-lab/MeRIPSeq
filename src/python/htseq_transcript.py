#!/usr/bin/env python
import HTSeq
import sys, getopt
import numpy
from collections import defaultdict


def calculate_htseq(annotation_file, input_file, output_file):
    """
    Calculation of the Median transcript value
    :param dataName:
    :return:
    """
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

    # Count reads for all exon
    gtf_transcript_file = HTSeq.GFF_Reader(annotation_file)

    # Count reads for all transcripts
    exons = dict()
    gene_to_exon = defaultdict(list)
    index = 1
    # load GTF file with all transcripts
    for feature in gtf_transcript_file:
        if feature.type == 'exon':
            gene_id = feature.attr['transcript_name']
            exon_id = gene_id + str(feature.attr['exon_number'])
            #print(exon_id,feature.iv)
            exons[exon_id] = feature.iv
            gene_to_exon[gene_id].append(exon_id)
            index += 1

    print('Read exon GTF '+str(len(exons))+' exons loaded')

    with open(output_file, "w") as htseqfile:
        for gene in gene_to_exon.keys():
            values = []
            lentho = []
            for exon in gene_to_exon[gene]:
                exoniv = exons[exon]
                lentho.append(exoniv)
                for iv, value in cvg[exoniv].steps():
                    for i in range(iv.length):
                        values.append(value)

            medianvalue = numpy.median(numpy.array(values))
            #if medianvalue > 0:
            #    print(medianvalue, values)
            #    print(medianvalue, len(values), lentho)
            row = str(gene)+'\t'+str(medianvalue)+'\n'
            #print(row)
            htseqfile.write(row)
    print('Median transcript expression calculated')


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
        opts, args = getopt.getopt(argv,"ha:i:o:", ["annotation=", "input=", "output="])
    except getopt.GetoptError:
        print('Cannot run command - Help: htseq_transcript.py -a <annotation> -i <input> -o <ouput>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('htseq_transcript.py -a <annotation> -i <input> -o <ouput>')
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