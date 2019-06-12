import HTSeq
import sys
import math
from utils import m6a_utils
import csv
from collections import defaultdict
import glob, os
from sets import Set
from Bio import SeqIO
from Bio import ExPASy
import re
from Bio import Entrez
from os import listdir
from os.path import isfile, join



PATH = '/Users/cbecavin/Documents/m6aAkker/'

PATH_ANNOT = PATH + 'Genome/'
PATH_DIFFCHIP = PATH + 'DiffChipSeq/'
PATH_DIFFBINDCOUNT = PATH + 'DiffChipSeq/DiffCount'
PATH_IGV = PATH + 'igv/'

GENE_TYPE = ['protein_coding','processed_pseudogene','lincRNA','TEC','unprocessed_pseudogene']


# Go through GTF file from GenCode and replace every Gene_id by their uniprot ID
# def parseGTFqdqz_annotation():
#     transcr_to_entrez = dict()
#     entrez_set = set()
#     for row in open(PATH + 'Genome/gencodeVM13/gencode.vM13.metadata.EntrezGene', 'rU'):
#         transcr = row.split('\t')[0]
#         entrez = row.split('\t')[1].strip()
#         entrez_set.add(entrez)
#         transcr_to_entrez[transcr] = entrez
#         # get all UniprotIDs for each EntrezID
#         file_uniprot_ID.readline()
#         entrezid_to_uniprotid = dict()
#         for row in file_uniprot_ID:
#             entrezid = row.split('\t')[0]
#             if not entrezid_to_uniprotid.has_key(entrezid):
#                 entrezid_to_uniprotid[entrezid] = [row.split('\t')[1]]
#             else:
#                 uniprot_info = entrezid_to_uniprotid[entrezid]
#                 uniprot_info.append(row.split('\t')[1])
#                 entrezid_to_uniprotid[entrezid] = uniprot_info
#         print('Add ', len(entrezid_to_uniprotid), ' uniprotID')
#     print('Transcript with entrez_id: '+str(len(transcr_to_entrez)) + ' nb entrez_id: '+str(len(entrez_set)))
#     transcr_to_swissprot = dict()
#     set_swissport = set()
#     for row in open(PATH + 'Genome/gencodeVM13/gencode.vM13.metadata.SwissProt', 'rU'):
#         transcr = row.split('\t')[0]
#         swissprot = row.split('\t')[1].strip()
#         set_swissport.add(swissprot)
#         transcr_to_swissprot[transcr] = swissprot
#     print('Transcript with swissprot: '+str(len(transcr_to_swissprot)) + ' nb entrez_id: '+str(len(set_swissport)))
#
#     headers = ['Transcript_ID','begin','end','chr','strand','type','gene_begin','gene_end','exon_number',
#                'transcript_name','gene_name','gene_id','EntrezID','SwissProtID']
#     with open(PATH + 'Genome/gencodeVM13/gencode.vM13.annotation.gtf', 'rU') as annotation_gtf, \
#             open(PATH_ANNOT + "gencodeVM13/DAVID.entrezID.to.uniprotID.txt", "rU") as file_uniprot_ID, \
#             open(PATH + 'Genome/gencodeVM13/gencode.vM13.annotation.txt', 'w') as annotation_txt:
#         gtf_peaks_file = HTSeq.GFF_Reader(annotation_gtf)
#         retained_gene = 0
#         number_gene = 0
#         annotation_txt.write('\t'.join(headers)+'\n')
#         for feature in gtf_peaks_file:
#             if(feature.type=='gene'):
#                 pos_gene = [str(feature.iv.start),str(feature.iv.end)]
#             if ('transcript_id' in feature.attr) and (feature.type == 'transcript'):
#                 number_gene += 1
#                 transcript_id = feature.attr['transcript_id']
#                 row = [transcript_id,str(feature.iv.start),str(feature.iv.end),str(feature.iv.chrom),
#                        str(feature.iv.strand),feature.type,pos_gene[0],pos_gene[1]]
#                 if 'exon_number' in feature.attr:
#                     row.append(feature.attr['exon_number'])
#                 else:
#                     row.append('')
#
#                 if 'transcript_name' in feature.attr:
#                     row.append(feature.attr['transcript_name'])
#                 else:
#                     row.append('')
#
#                 if 'gene_name' in feature.attr:
#                     row.append(feature.attr['gene_name'])
#                     row.append(feature.attr['gene_id'])
#                 else:
#                     row.append('')
#                     row.append('')
#
#                 if transcript_id in transcr_to_entrez.keys():
#                     row.append(transcr_to_entrez[transcript_id])
#                 else:
#                     row.append('')
#                 if transcript_id in transcr_to_swissprot.keys():
#                     row.append(transcr_to_swissprot[transcript_id])
#                 else:
#                     row.append('')
#
#                 # print(feature.attr['gene_id'],feature.attr['gene_name'])
#                 # print(feature.attr['gene_type'])
#                 annotation_txt.write('\t'.join(row)+'\n')
#
#                 if number_gene % 10000 == 0:
#                     print(feature.get_gff_line().strip())
#         print(number_gene)

#
#  keep only transcripts in the gtf
#
def create_gene_and_exon_gtf():
    '''
    Keep only gene in the GTF for count and median expression calculation
    :return:
    '''
    gtf_file = HTSeq.GFF_Reader(PATH_ANNOT + "/gencodeVM13/gencode.vM13.annotation.gtf", end_included=True)
    with open(PATH_ANNOT + "/gencodeVM13/gencode.vM13.annotation.gene.gtf", "w") as cleanGTF:
        for feature in gtf_file:
            if feature.type == 'gene':
                #print(feature.get_gff_line())
                cleanGTF.write(feature.get_gff_line())

    gtf_file = HTSeq.GFF_Reader(PATH_ANNOT + "/gencodeVM13/gencode.vM13.annotation.gtf", end_included=True)
    with open(PATH_ANNOT + "/gencodeVM13/gencode.vM13.annotation.exon.gtf", "w") as cleanGTF:
        for feature in gtf_file:
            if feature.type == 'gene':
                #print(feature.get_gff_line())
                cleanGTF.write(feature.get_gff_line())
            if feature.type == 'exon':
                #print(feature.get_gff_line())
                cleanGTF.write(feature.get_gff_line())
    gtf_file = HTSeq.GFF_Reader(PATH_ANNOT + "/gencodeVM13/gencode.vM13.annotation.gtf", end_included=True)
    with open(PATH_ANNOT + "/gencodeVM13/gencode.vM13.annotation.transcript.gtf", "w") as cleanGTF:
        for feature in gtf_file:
            if feature.type == 'transcript':
                # print(feature.get_gff_line())
                cleanGTF.write(feature.get_gff_line())
            if feature.type == 'exon':
                # print(feature.get_gff_line())
                cleanGTF.write(feature.get_gff_line())


def create_utr_gtf():
    with open(PATH+'Genome/gencodeVM13/gencode.vM13.annotation.gtf','r') as annotation_gtf, \
            open(PATH + 'Genome/gencodeVM13/gencode.vM13.annotation.3utr.gtf', 'w') as utr_3_gtf, \
            open(PATH + 'Genome/gencodeVM13/gencode.vM13.annotation.5utr.gtf', 'w') as utr_5_gtf, \
            open(PATH + 'Genome/gencodeVM13/utr5_description.txt', 'w') as utr_descr_5, \
            open(PATH + 'Genome/gencodeVM13/utr3_description.txt', 'w') as utr_descr_3:
        gtf_peaks_file = HTSeq.GFF_Reader(annotation_gtf)
        transcript_to_details = defaultdict(list)
        utr_descr_5.write('utr_id\tGene_ID\n')
        utr_descr_3.write('utr_id\tGene_ID\n')
        k=0
        for feature in gtf_peaks_file:
            if 'gene_id' in feature.attr.keys():
                if feature.attr['gene_type'] in GENE_TYPE:
                    transcr_id = feature.attr['gene_id']
                    #print(feature.get_gff_line())
                    transcript_to_details[transcr_id].append(feature)
                    k+=1
                    if k%10000 == 0:
                        print(k)
        no_utr = 0
        for transcr_id in transcript_to_details:
            #print(transcr_id+'\t'+str(transcript_to_details[transcr_id]))
            list_feature = transcript_to_details[transcr_id]
            list_utr = []

            for feature in list_feature:
                if feature.type == 'gene':
                    transcript = feature
                elif feature.type == 'UTR':
                    list_utr.append(feature)

            # get gene
            if len(list_utr)==0:
                print(transcr_id,'no UTR')
                #no_utr += 1
                #for feature in list_feature:
                #    print(transcript.name, transcript.attr['gene_type'], feature.type)
            else:
                utr_5_index = 1
                utr_3_index = 1
                for feature in list_feature:
                    type = 'exon'
                    diffStart = float(feature.iv.start - transcript.iv.start)
                    diffEnd = float(feature.iv.end - transcript.iv.start)
                    if transcript.iv.strand == '-':
                        diffStart = float(transcript.iv.end - feature.iv.end)
                        diffEnd = float(transcript.iv.end - feature.iv.start)

                    relative_pos_begin = diffStart / transcript.iv.length
                    relative_pos_end = diffEnd / transcript.iv.length
                    type = feature.type
                    if feature.type == 'UTR':
                        if relative_pos_begin < 0.5:
                            type = '5UTR'
                            utr_id = transcr_id + "_5UTR_"+str(utr_5_index)
                            utr_5_index += 1
                            feature.attr['utr_id'] = utr_id
                            if (feature.iv.length > 10) and (feature.iv.length < 30000):
                                utr_5_gtf.write(feature.get_gff_line())
                                utr_descr_5.write(utr_id+'\t'+transcr_id+'\n')
                        else:
                            type = '3UTR'
                            utr_id = transcr_id + "_3UTR_"+str(utr_3_index)
                            utr_3_index += 1
                            feature.attr['utr_id'] = utr_id
                            if (feature.iv.length > 10) and (feature.iv.length < 30000):
                                utr_3_gtf.write(feature.get_gff_line())
                                utr_descr_3.write(utr_id + '\t' + transcr_id + '\n')

                        print(type, transcr_id,feature.iv.strand,feature.type, relative_pos_begin,relative_pos_end)

        print('no UTR total: ',no_utr,len(transcript_to_details))


def create_utr_fasta():
    with open(PATH+'Genome/gencodeVM13/gencode.vM13.annotation.gtf','r') as annotation_gtf, \
            open(PATH + 'Genome/gencodeVM13/gencode.vM13.annotation.3utr.fasta', 'w') as utr_3_gtf, \
            open(PATH + 'Genome/gencodeVM13/gencode.vM13.annotation.5utr.fasta', 'w') as utr_5_gtf:
        gtf_peaks_file = HTSeq.GFF_Reader(annotation_gtf)
        mouse_seq = m6a_utils.read_mouse_seq()
        transcript_to_details = defaultdict(list)
        k=0
        for feature in gtf_peaks_file:
            if 'gene_id' in feature.attr.keys():
                if feature.attr['gene_type'] == 'protein_coding':
                    transcr_id = feature.attr['gene_id']
                    #print(feature.get_gff_line())
                    transcript_to_details[transcr_id].append(feature)
                    k+=1
                    if k%10000 == 0:
                        print(k)
        no_utr = 0
        for transcr_id in transcript_to_details:
            #print(transcr_id+'\t'+str(transcript_to_details[transcr_id]))
            list_feature = transcript_to_details[transcr_id]
            list_utr = []

            for feature in list_feature:
                if feature.type == 'gene':
                    transcript = feature
                elif feature.type == 'UTR':
                    if feature.iv.length > 10 & feature.iv.length < 30000:
                        list_utr.append(feature)

            # get gene
            if len(list_utr)==0:
                print(transcr_id,'no UTR')
                #no_utr += 1
                #for feature in list_feature:
                #    print(transcript.name, transcript.attr['gene_type'], feature.type)
            else:
                utr_5_index = 1
                utr_3_index = 1
                for feature in list_feature:
                    type = 'exon'
                    diffStart = float(feature.iv.start - transcript.iv.start)
                    diffEnd = float(feature.iv.end - transcript.iv.start)
                    if transcript.iv.strand == '-':
                        diffStart = float(transcript.iv.end - feature.iv.end)
                        diffEnd = float(transcript.iv.end - feature.iv.start)

                    relative_pos_begin = diffStart / transcript.iv.length
                    relative_pos_end = diffEnd / transcript.iv.length
                    type = feature.type
                    if feature.type == 'UTR':
                        if relative_pos_begin < 0.5:
                            type = '5UTR'
                            utr_id = transcr_id + "_5UTR_"+str(utr_5_index)
                            utr_5_index += 1
                            feature.attr['utr_id'] = utr_id
                            sequence_for_motif = mouse_seq[feature.iv.chrom][feature.iv.start:feature.iv.end].seq
                            #print(len(str(sequence_for_motif)), str(sequence_for_motif))
                            if (len(str(sequence_for_motif)) > 10) and (len(str(sequence_for_motif)) < 30000):
                                utr_5_gtf.write('>'+utr_id+'\n'+str(sequence_for_motif)+'\n')

                        else:
                            type = '3UTR'
                            utr_id = transcr_id + "_3UTR_"+str(utr_3_index)
                            utr_3_index += 1
                            feature.attr['utr_id'] = utr_id
                            sequence_for_motif = mouse_seq[feature.iv.chrom][feature.iv.start:feature.iv.end].seq
                            if (len(str(sequence_for_motif)) > 10) and (len(str(sequence_for_motif)) < 30000):
                                utr_3_gtf.write('>'+utr_id+'\n'+str(sequence_for_motif)+'\n')

                        #print(type, transcr_id,feature.iv.strand,feature.type, relative_pos_begin,relative_pos_end)

        print('no UTR total: ',no_utr,len(transcript_to_details))



#
#    Prepare window file GTF
#
def create_sliding_gene_window_GTF(windowSize):
    overlap=windowSize/2
    gtf_file = HTSeq.GFF_Reader(PATH_ANNOT + "/gencodeVM13/gencode.vM13.annotation.gene.gtf", end_included=True)
    windows = HTSeq.GenomicArrayOfSets( "auto", stranded=True)
    transcriptID = 1
    with open(PATH_ANNOT + "/gencodeVM13/gencode.vM13.gene.slidingwindow.gtf", "w") as slidingGTF:
        for feature in gtf_file:
            interval=feature.iv
            transcriptID += 1
            if transcriptID%1000 == 0:
                print('Gene: '+str(transcriptID)+'/ 100000')
            windowID = 1
            if interval.strand == '+':
                begin = interval.start_d
                end = begin + windowSize
                while end<interval.end_d:
                    window=HTSeq.GenomicInterval(interval.chrom, begin, end+1, interval.strand)
                    featureWindow=HTSeq.GenomicFeature(feature.attr['gene_id']+'_window_'+
                                                       feature.attr['gene_type']+'_'+str(windowID),
                                                         'window', window)
                    windowID+=1
                    #print(featureWindow.get_gff_line())
                    begin+=overlap
                    end=begin+windowSize
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
                end = begin - windowSize
                while end > interval.end_d:
                    window = HTSeq.GenomicInterval(interval.chrom, end, begin+1, interval.strand)
                    featureWindow = HTSeq.GenomicFeature(feature.attr['gene_id'] +'_window_'+
                                                         feature.attr['gene_type']+'_'+str(windowID),
                                                         'window', window)
                    windowID += 1
                    #print(featureWindow.get_gff_line())
                    begin -= overlap
                    end = begin - windowSize
                    slidingGTF.write(featureWindow.get_gff_line())

                end = interval.end_d
                window = HTSeq.GenomicInterval(interval.chrom, end+1, begin, interval.strand)
                featureWindow = HTSeq.GenomicFeature(feature.attr['gene_id'] + '_window_' +
                                                         feature.attr['gene_type'] + '_' + str(windowID),
                                                         'window', window)
                slidingGTF.write(featureWindow.get_gff_line())


def create_sliding_exon_window_GTF(windowSize):
    overlap=windowSize/2
    gtf_file = HTSeq.GFF_Reader(PATH_ANNOT + "/gencodeVM13/gencode.vM13.annotation.exon.gtf", end_included=True)
    windows = HTSeq.GenomicArrayOfSets( "auto", stranded=True)
    transcriptID = 1
    with open(PATH_ANNOT + "/gencodeVM13/gencode.vM13.exon.slidingwindow.gtf", "w") as slidingGTF:
        for feature in gtf_file:
            if feature.type == 'exon':
                interval=feature.iv
                transcriptID += 1
                if transcriptID%1000 == 0:
                    print('Gene: '+str(transcriptID)+'/ 100000')
                windowID = 1
                if interval.strand == '+':
                    begin = interval.start_d
                    end = begin + windowSize
                    while end<interval.end_d:
                        window=HTSeq.GenomicInterval(interval.chrom, begin, end+1, interval.strand)
                        featureWindow=HTSeq.GenomicFeature(feature.attr['transcript_name']+'_'+feature.attr['exon_number']
                                                           +'_window_'+
                                                           feature.attr['gene_type']+'_'+str(windowID),
                                                             'window', window)
                        windowID+=1
                        #print(featureWindow.get_gff_line())
                        begin+=overlap
                        end=begin+windowSize
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
                    end = begin - windowSize
                    while end > interval.end_d:
                        window = HTSeq.GenomicInterval(interval.chrom, end, begin+1, interval.strand)
                        featureWindow = HTSeq.GenomicFeature(feature.attr['transcript_name'] +'_'+feature.attr['exon_number']
                                                             +'_window_'+
                                                             feature.attr['gene_type']+'_'+str(windowID),
                                                             'window', window)
                        windowID += 1
                        #print(featureWindow.get_gff_line())
                        begin -= overlap
                        end = begin - windowSize
                        slidingGTF.write(featureWindow.get_gff_line())

                    end = interval.end_d
                    window = HTSeq.GenomicInterval(interval.chrom, end+1, begin, interval.strand)
                    featureWindow = HTSeq.GenomicFeature(feature.attr['transcript_name']+'_'+feature.attr['exon_number']+
                                                         '_window_' +
                                                             feature.attr['gene_type'] + '_' + str(windowID),
                                                             'window', window)
                    if window.length > 1:
                        slidingGTF.write(featureWindow.get_gff_line())


# From DAVID website I created a list of all ENTREZID to UniprotID
# I use ti to add UniprotID in the BioMart annotation table
def parse_GTF_transcript_annotation():
    with open(PATH_ANNOT + "gencodeVM13/gencode.vM13.metadata.EntrezGene", "rU") as file_entrez, \
            open(PATH_ANNOT + "gencodeVM13/DAVID.entrezID.to.uniprotID.txt", "rU") as file_uniprot_ID, \
            open(PATH_ANNOT + "gencodeVM13/gencode.vM13.annotation.gtf", "rU") as file_transcript, \
            open(PATH_ANNOT + "gencodeVM13/gencode.vM13.annotation.txt", "w") as new_annotation_file:
        # get all UniprotIDs for each EntrezID
        file_uniprot_ID.readline()
        entrezid_to_uniprotid = dict()
        entrezid_to_description = dict()
        for row in file_uniprot_ID:
            entrezid = row.split('\t')[0]
            description = row.split('\t')[-1].strip()
            entrezid_to_description[entrezid] = description
            if not entrezid_to_uniprotid.has_key(entrezid):
                entrezid_to_uniprotid[entrezid] = [row.split('\t')[1]]
            else:
                uniprot_info = entrezid_to_uniprotid[entrezid]
                uniprot_info.append(row.split('\t')[1])
                entrezid_to_uniprotid[entrezid] = uniprot_info
        print('Add ',len(entrezid_to_uniprotid),' uniprotID')


        # get transcriptID to entrez information
        transcriptid_to_entrezid = dict()
        for row in file_entrez:
            transcriptid = row.split('\t')[0]
            entrezid = row.split('\t')[1].strip()
            entrezinformation = [entrezid]
            if entrezid in entrezid_to_uniprotid:
                entrezinformation.append(entrezid_to_uniprotid[entrezid][0])
                entrezinformation.append(';'.join(entrezid_to_uniprotid[entrezid]))
            else:
                entrezinformation.append('')
                entrezinformation.append('')

            if entrezid in entrezid_to_description:
                entrezinformation.append(entrezid_to_description[entrezid])
            else:
                entrezinformation.append('')

            #print(transcriptid,entrezinformation)
            transcriptid_to_entrezid[transcriptid] = entrezinformation
        print('Add ', len(transcriptid_to_entrezid), ' transcript id')


        # update annotation file using entrez information
        headers = ['Transcript_ID', 'begin', 'end', 'chr', 'strand', 'type', 'gene_begin', 'gene_end',
                   'transcript_name', 'gene_name', 'gene_id', 'EntrezID', 'UniprotID', 'UniprotIDs', 'Description']
        new_annotation_file.write('\t'.join(headers)+'\n')

        gtf_peaks_file = HTSeq.GFF_Reader(file_transcript)
        retained_gene = 0
        number_gene = 0
        for feature in gtf_peaks_file:
            if (feature.type == 'gene'):
                pos_gene = [str(feature.iv.start), str(feature.iv.end)]
            if ('transcript_id' in feature.attr) and (feature.type == 'transcript'):
                number_gene += 1
                transcript_id = feature.attr['transcript_id']
                row = [transcript_id, str(feature.iv.start), str(feature.iv.end), str(feature.iv.chrom),
                       str(feature.iv.strand), feature.attr['gene_type'], pos_gene[0], pos_gene[1]]

                if 'transcript_name' in feature.attr:
                    row.append(feature.attr['transcript_name'])
                else:
                    row.append('')

                if 'gene_name' in feature.attr:
                    row.append(feature.attr['gene_name'])
                    row.append(feature.attr['gene_id'])
                else:
                    row.append('')
                    row.append('')

                if transcript_id in transcriptid_to_entrezid:
                    for entry in transcriptid_to_entrezid[transcript_id]:
                        row.append(entry)

                # print(feature.attr['gene_id'],feature.attr['gene_name'])
                # print(feature.attr['gene_type'])
                new_annotation_file.write('\t'.join(row) + '\n')

                if number_gene % 10000 == 0:
                    print(feature.get_gff_line().strip())
        print(number_gene)

        print('Created annotation : gencode.vM13.annotation.txt')


def create_gene_annotation():
    with open(PATH_ANNOT + "gencodeVM13/gencode.vM13.annotation.txt", "rU") as file_annotation, \
            open(PATH_ANNOT + "gencodeVM13/gencode.vM13.gene.annotation.txt", "w") as gene_annotation:
        gene_file = set()
        csv_annot = csv.DictReader(file_annotation, delimiter='\t')
        headers = ['gene_id', 'chr', 'gene_begin', 'gene_end', 'strand', 'type', 'gene_name',
                   'EntrezID', 'UniprotID', 'UniprotIDs', 'Description']
        gene_annotation.write('\t'.join(headers)+'\n')
        for row_annot in csv_annot:
            row = ''
            for header in headers:
                #print(header,row_annot)
                if type(row_annot[header]) == str:
                    row += row_annot[header] + '\t'
            gene_file.add(row)

        for row in gene_file:
            gene_annotation.write(row.strip() + '\n')


def parse_gene_description():
    with open(PATH+'Genome/gencode.vM13.annotation.gtf','r') as annotation_gtf, \
            open(PATH + 'Genome/gencode.vM13.gene.description.txt', 'w') as transcript_description:
        gtf_peaks_file = HTSeq.GFF_Reader(annotation_gtf)
        transcript_to_details = defaultdict(list)
        k=0
        for feature in gtf_peaks_file:
            if 'gene_id' in feature.attr.keys():
                if feature.attr['gene_type'] in GENE_TYPE:
                    transcr_id = feature.attr['gene_id']
                    #print(feature.get_gff_line())
                    transcript_to_details[transcr_id].append(feature)
                    k+=1
                    if k%10000 == 0:
                        print(k)
        no_utr = 0
        for transcr_id in transcript_to_details:
            #print(transcr_id+'\t'+str(transcript_to_details[transcr_id]))
            list_feature = transcript_to_details[transcr_id]
            list_utr = []

            for feature in list_feature:
                if feature.type == 'gene':
                    transcript = feature
                elif feature.type == 'UTR':
                    list_utr.append(feature)

            # get gene
            if len(list_utr)==0:
                print(transcr_id,'no UTR')
                #no_utr += 1
                #for feature in list_feature:
                #    print(transcript.name, transcript.attr['gene_type'], feature.type)
            else:
                for feature in list_feature:
                    type = 'exon'
                    diffStart = float(feature.iv.start - transcript.iv.start)
                    diffEnd = float(feature.iv.end - transcript.iv.start)
                    if transcript.iv.strand == '-':
                        diffStart = float(transcript.iv.end - feature.iv.end)
                        diffEnd = float(transcript.iv.end - feature.iv.start)

                    relative_pos_begin = diffStart / transcript.iv.length
                    relative_pos_end = diffEnd / transcript.iv.length
                    type = feature.type
                    if feature.type == 'UTR':
                        if relative_pos_begin < 0.5:
                            type = '5UTR'
                        else:
                            type = '3UTR'
                    transcript_description.write(transcr_id+'\t'+type+'\t'+
                                                 str(feature.iv.start)+'\t'+str(feature.iv.end)+'\n')
                    print(type, transcr_id,feature.iv.strand,feature.type, relative_pos_begin,relative_pos_end)

        print('no UTR total: ',no_utr,len(transcript_to_details))


            #transcript_description.write(transcript_id+'\t'+transcript_to_details[transcript_id]+'\n')


def get_rna_fasta():
#     gtf_file = HTSeq.GFF_Reader('/Users/cbecavin/Documents/RNABindingProtein/Genome/NC_003210.gff', end_included=True)
# #    gtf_file = HTSeq.GFF_Reader(PATH_ANNOT + "/gencodeVM13/gencode.vM13.annotation.gtf", end_included=True)
#     with open('/Users/cbecavin/Documents/RNABindingProtein/Genome/NC_003210.rRNA.bed', "w") as bed_file:
#         for feature in gtf_file:
#             if feature.type == 'rRNA':
#                 print(feature)
#                 bed_file.write(feature.iv.chrom + '\t' + str(feature.iv.start) + '\t' + str(feature.iv.end) +
#                                '\t1\t1\t'+feature.iv.strand+'\n')
    # bedtools merge
    #mouse_seq = m6a_utils.read_mouse_seq()
    mouse_seq = m6a_utils.read_genome('/Users/cbecavin/Documents/RNABindingProtein/Genome/NC_003210.fna')
    with open('/Users/cbecavin/Documents/RNABindingProtein/Genome/NC_003210.rRNA.bed', "r") as bed_file, \
            open('/Users/cbecavin/Documents/RNABindingProtein/Genome/NC_003210.rRNA.fna', "w") as fasta_file:
        fasta_file.write('>NC_003210_rRNA\n')
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
                print(str(compl))



def get_type_gene():
    # gtf_file = HTSeq.GFF_Reader(PATH_ANNOT + "/gencode.vM13.annotation.gtf", end_included=True)
    # with open(PATH_ANNOT + "/gencode.vM13.gene.type.txt", "w") as bed_file:
    #     for feature in gtf_file:
    #         if feature.type == 'gene':
    #             #print(feature.get_gff_line())
    #             bed_file.write(feature.attr['gene_id'] + '\t' + feature.attr['gene_type'] + '\n')

    types = set()
    gene_types = dict()
    keep_type = ['Mt_rRNA', 'Mt_tRNA', 'rRNA', 'protein_coding', 'processed_pseudogene', 'lincRNA', 'miRNA']
    with open(PATH_ANNOT + "/gencode.vM13.gene.type.txt", "r") as bed_file, \
            open(PATH_ANNOT + "/gencode.vM13.gene.type.lite.txt", "w") as new_file:
        for row in bed_file:
            type = row.split('\t')[1].strip()
            if type in keep_type:
                gene_types[row.split('\t')[0].strip()] = row.split('\t')[1].strip()
                new_file.write(row)
        print(types)

    path_expr = PATH + 'Expression/'
    for file in glob.glob(path_expr+'*All*.txt'):
        print(file)
        type_count = defaultdict(int)
        for row in open(file,"r"):
            #print(row)
            gene_id = row.split('\t')[0].strip()
            if gene_id in gene_types.keys():
                gene_type = gene_types[gene_id]
                count_gene = int(row.split('\t')[1].strip())
                type_count[gene_type] += count_gene
        print(type_count)


# Prepare annotation files

# Download : gencode.vM13.annotation.gtf
# Download : gencode.vM13.metadata.EntrezGene

    # create gencode.vM13.annotation.txt (final annotation file)
# parse_GTF_transcript_annotation()
# create_gene_annotation()
#
# create_gene_and_exon_gtf()
#
# create_gene_and_exon_bed()
#
# #create_utr_gtf()
# #create_utr_fasta()
#
# create_sliding_gene_window_GTF(100)
# create_sliding_exon_window_GTF(100)
#
parse_gene_description()

#get_rna_fasta()

#get_type_gene()

    # Get entrez ID for all transcripts : gencode.vM13.metadata.EntrezGene
    #       90319 transcript
    #       20729 entrezID
    # Extract uniprotID using DAVID pathway website : DAVID.entrezID.to.uniprotID.txt
    #       20360 entrezID for which we have a UniprotID
    #       70317 uniprotID (multi uniprotID for one entrezID
    # Get annotation and GO information from DAVID pathway website using EntrezID: DAVID.gencode.entrez.vM13.metadata
    #       20677 EntrezID
    # Get Hallmark H gene set information from http://bioinf.wehi.edu.au/software/MSigDB/   -> HallmarkH_genesets.txt
    #       000 EntrezID
    # Regroup in one table all EntrezID information gathered
    ################################
    #
    # FINAL FILE : gencode.vM13.annotation.entrez.uniprot.txt
    #
    ################################
#parse_entrez_id()
   # Make a file where we have the relative position of UTR, exon, intron, etc. on a transcript
#parse_transcript_description()
    # Modify GTF file by adding UniprotID for HTSeq calculation
    #       20344 uniprotID counted in HTSeq

#parseGTF_HTSeq()
#test_gene_annotation()
#parse_hallmark()


