#!/usr/bin/env python

import HTSeq
import subprocess
from Bio import SeqIO

PATH = '/Users/cbecavin/Documents/m6aAkker/'
PATH_MAPPING = 'Volumes/USBHUB/m6aAkker/Mapping/'
PATH_PEAKS = PATH + 'PeakDetection/'
PATH_EXPRESSION = PATH + 'DiffExpression/Star/'
PATH_PEAKS_POM = PATH + 'PeakDetection/POMScore/'
PATH_DIFF = PATH + 'DiffChipSeq/'
PATH_ExpDesign = PATH + 'ExpDesign/'
PATH_ANNOT = PATH + 'Genome/'
PATH_IGV = '/Users/cbecavin/Documents/IGV/'

MOTIF_METH = {'GAACA': 2, 'GGACA': 3, 'GAACT': 5,'GGACT': 8}
WINDOW_SIZE = 100

MIN_COVERAGE = 9
MIN_FOLD = 4


DESIGN_ALL = "m6aExpDesign_All_datasets"
DESIGN_LIVER = "m6aExpDesign_Liver_all"
DESIGN_CMGF = 'm6aExpDesign_CM_GF'
DESIGN_SEQ78 = 'm6aExpDesign_CM_GF_Seq78'
DESIGN_CAECUM = 'm6aExpDesign_Caecum'


def shellCommand(command):
    p = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return(stdout, stderr, p.returncode)


def read_mouse_seq():
    with open(PATH_ANNOT + 'GRCm38.primary_assembly.genome.fa', 'r') as mousefile:
        mouse_seq = SeqIO.to_dict(SeqIO.parse(mousefile, "fasta"))
        print('Mouse genome read')
        return mouse_seq

def read_genome(filename):
    with open(filename, 'r') as mousefile:
        mouse_seq = SeqIO.to_dict(SeqIO.parse(mousefile, "fasta"))
        print('Genome read')
        return mouse_seq


def read_transcript_gtf():
    gtf_transcript_file = HTSeq.GFF_Reader(PATH_ANNOT + "gencodeVM13/gencode.vM13.annotation.transcript.gtf")
    transcipt_id_to_attr = dict()
    for feature in gtf_transcript_file:
        info_transcript = [feature.iv.start, feature.iv.end, feature.iv.strand, feature.iv.chrom,
                           feature.attr['gene_name'], feature.attr['transcript_type'], feature.attr['gene_type'],
                           feature.attr['transcript_id']]
        transcipt_id_to_attr[feature.attr['transcript_name']] = info_transcript
    return transcipt_id_to_attr


def get_protein_coding_gtf():
    gtf_file = HTSeq.GFF_Reader(PATH_ANNOT + "gencodeVM13/gencode.vM13.annotation.transcript.gtf", end_included=True)
    windows = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    with open(PATH_ANNOT + "gencode.vM9.annotation.protein_coding.gtf", "w") as slidingGTF:
        for feature in gtf_file:
            if feature.attr['gene_type'] == "protein_coding":
                slidingGTF.write(feature.get_gff_line())