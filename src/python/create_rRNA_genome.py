#!/usr/bin/env python
import HTSeq
import csv
from collections import defaultdict
import sys, getopt
from Bio import SeqIO


def get_rna_bed(path_annot, annotation_file):
    print('Select rRNA in bed file')
    gtf_file = HTSeq.GFF_Reader(path_annot + '/' + annotation_file + '.annotation.gtf', end_included=True)
    with open(path_annot + '/' + annotation_file + '.annotation.rRNA.bed', "w") as bed_file:
        for feature in gtf_file:
            if feature.attr['gene_type'] == 'rRNA':
                #print(feature.get_gff_line())
                bed_file.write(feature.iv.chrom + '\t' + str(feature.iv.start) + '\t' + str(feature.iv.end) +
                               '\t1\t1\t'+feature.iv.strand+'\n')


def get_rna_fasta(path_annot, annotation_file, genome_file):
    print('Save rRNA in fasta file')
    mouse_seq = read_mouse_seq(path_annot, genome_file)
    with open(path_annot + '/' + annotation_file + '.annotation.rRNA_merge.bed', "r") as bed_file, \
            open(path_annot + '/' + annotation_file + '.rRNA.fasta', "w") as fasta_file:
        fasta_file.write('>mm10_GRCm38_rRNA\n')
        for row in bed_file:
            chromo = row.split('\t')[0]
            begin = int(row.split('\t')[1])
            end = int(row.split('\t')[2])
            strand = row.split('\t')[3]
            sequence = mouse_seq[chromo][begin:end].seq
            if strand == '+':
                fasta_file.write(str(sequence))
            else:
                compl = sequence.reverse_complement()
                fasta_file.write(str(compl))
                #print(str(compl))


def read_mouse_seq(path, genome_file):
    with open(path + genome_file, 'r') as mousefile:
        mouse_seq = SeqIO.to_dict(SeqIO.parse(mousefile, "fasta"))
        print('Mouse genome read')
        return mouse_seq


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
        opts, args = getopt.getopt(argv,"hp:a:g:o:", ["path=", "annotation=","genome=","option="])
    except getopt.GetoptError:
        print('Cannot run command - Help: create_rRNA_genome.py -p <path> -a <annotation_file> -g <genome> -o <option>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('create_rRNA_genome.py -p <path> -a <annotation_file> -g <genome> -o <option>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path_annot = arg
        elif opt in ("-a", "--annotation"):
            annotation_file = arg
        elif opt in ("-g", "--genome"):
            genome_file = arg
        elif opt in ("-o", "--option"):
            option = arg
    
    if option=='bed':
        get_rna_bed(path_annot, annotation_file)
    else:
        get_rna_fasta(path_annot, annotation_file, genome_file)


if __name__ == "__main__":
    main(sys.argv[1:])

