#!/usr/bin/env python
import csv
import sys, getopt
import os
import HTSeq
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import sequence_python


MOTIF_METH = {'GGACC': 6, 'CCACC': 6, 'CGACC': 6, 'GCACC': 6, 
                'GCACT': 6, 'CCACT': 6, 'GGACT': 6,'CGACT': 6,
                'GCACA': 6, 'CCACA': 6, 'GGACA': 6,'CGACA': 6}
GENE_TYPE = ['protein_coding']
#GENE_TYPE = ['protein_coding','processed_pseudogene','lincRNA','TEC','unprocessed_pseudogene']


def annotate_peaks(path, exp_design_name, mouse_seq, gene_list, annotation_filename, technique, bed_name):
    """
    Read list of windows_indexes and add :
        - annotation
        - motif presence score
        - relative position on transcript
        - overlapping refpeaks
        -

    INPUT: PATH_PEAKS + exp_design_name + '_windows_indexes.txt'
           PATH_PEAKS + exp_design_name + '_windows_list_Annot.txt'
    OUTPUT: PATH_PEAKS + exp_design_name + '_windows_indexes_Annot.txt'

    :param exp_design_name:
    :param mouse_seq:
    :param ref_peaks_list:
    :param PATH_PEAKS:
    :return:
    """
    PATH_ANNOT = path + '/Genome/'
    if technique == '' or technique == 'All':
        PATH_PEAKS = path + '/PeakDetection/Peaks'
        bed_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name +'.bed'
        annot_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '_Peaks.txt'
    else:
        PATH_PEAKS = path + '/PeakDetection/' + technique + '/'
        bed_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name +'.bed'
        annot_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '_Peaks.txt'


    with open(bed_filename, "rU") as bed_file, \
            open(annot_filename, "w") as annot_file:

        df_annot = pd.read_csv(PATH_ANNOT + annotation_filename, index_col=0, sep='\t')
        headers = ['WindowId','chromo_peak','begin_peak','end_peak','length_peak',
                   'Motif', 'Peak_Presence', 'Type_overlap', 'Ref_peak', 'Relat_pos','Gene_ID']
        headers.extend(df_annot.columns)
        annot_file.write('\t'.join(headers)+'\n')
        index_coding = 0
        index_nc = 0
        for row in bed_file:
            if technique == '':
                peak_technique = row.split('\t')[4].strip()
            else:
                peak_technique = technique
            peak = HTSeq.GenomicInterval(row.split('\t')[0], int(row.split('\t')[1]), int(row.split('\t')[2]), ".")
            gene_ids = []
            gene_ivs = []
            gene_types = []
            for gene_id in gene_list:
                gene = gene_list[gene_id]
                if peak.overlaps(gene.iv):
                    if technique == '':
                        if gene.attr['gene_type'] in GENE_TYPE:
                            # if gene.attr['gene_type'] == 'protein_coding' && gene.attr['gene_type'] in GENE_TYPE:
                            gene_ids.append(gene_id)
                            gene_ivs.append(gene.iv)
                            gene_types.append(gene.attr['gene_type'])
                    else:
                        gene_ids.append(gene_id)
                        gene_ivs.append(gene.iv)
                        gene_types.append(gene.attr['gene_type'])

            #print(len(transcript_ids),peak)
            if len(gene_ids) == 0:
                index_nc += 1
                window_id = 'Peak_NC_' + str(index_nc)
                if index_nc % 1000 == 0:
                    print(window_id)
                print('no ',peak)
            else:
                for i in range(len(gene_ids)):
                    gene_id = gene_ids[i]
                    gene_iv = gene_ivs[i]
                    gene_type = gene_types[i]
                    if gene_type == 'protein_coding':
                        index_coding += 1
                        window_id = 'Peak_' + str(index_coding)
                        if index_coding % 1000 == 0:
                            print(window_id)
                    else:
                        index_nc += 1
                        window_id = 'Peak_NC_' + str(index_nc)
                        if index_nc % 1000 == 0:
                            print(window_id)

                    # Calculate motif presence score
                    sequence_score = 0
                    sequence = mouse_seq[peak.chrom][peak.start:peak.end].seq
                    for motif in MOTIF_METH:
                        if motif in sequence:
                            sequence_score += MOTIF_METH[motif]

                    # relative position on transcript
                    relative_pos_start = 0.5
                    relative_pos_end = 0.5
                    if gene_id != '':
                        relative_pos_start = float(peak.start - gene_iv.start)
                        if gene_iv.strand == '-':
                            relative_pos_start = float(gene_iv.end - peak.end)
                        relative_pos_start = relative_pos_start / gene_iv.length

                        relative_pos_end = float(peak.end - gene_iv.start)
                        if gene_iv.strand == '-':
                            relative_pos_end = float(gene_iv.end - peak.start)
                        relative_pos_end = relative_pos_end / gene_iv.length
                    relat_pos = float(relative_pos_end - relative_pos_start)/2 + relative_pos_start

                    new_row = [window_id, peak.chrom, str(peak.start), str(peak.end), str(peak.length),
                               str(sequence_score), peak_technique, '_Over_', '_Ref_', str(relat_pos), gene_id]

                    if gene_id != '':
                        for i in range(0, len(df_annot.columns)):
                            header = df_annot.columns[i]
                            new_row.append(df_annot[header][gene_id])
                    else:
                        print('not found ' + window_id)

                    annot_file.write('\t'.join([str(i) for i in new_row]) + '\n')


def peak_relative_pos(path, exp_design_name, technique, annotation_file, bed_name):
    """
    Study the poisition in the transcript of peaks : UTR, Start_codon, CDS, exon, intron

    INPUT:  PATH + 'Genome/gencodeVM13/gencode.vM13.transcript.description.txt'
            PATH_PEAKS + exp_design_name + '_' + col_name + '_Summary_Gene.txt'

    OUTPUT: PATH_PEAKS + exp_design_name + '_' + col_name + '_Relat_pos.txt'

    :param col_name:
    :param exp_design_name:
    :param PATH_PEAKS:
    :return:
    """
    PATH_ANNOT = path + '/Genome/'
    if technique == '':
        PATH_PEAKS = path + '/PeakDetection/Peaks'
        relat_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '_Relat_pos.txt'
        annot_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '_Peaks.txt'
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '_Annot.txt'
    else:
        PATH_PEAKS = path + '/PeakDetection/' + technique
        relat_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '_Relat_pos.txt'
        annot_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '_Peaks.txt'
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '_Annot.txt'


    peak_to_overlap = defaultdict(list)
    with open(PATH_ANNOT + annotation_file, 'r') as transcript_description, \
            open(relat_filename, 'w') as relat_pos_file:
        transcript_to_details = defaultdict(list)
        headers = ['WindowId','Transcript_ID','Type','Pos']
        relat_pos_file.write('\t'.join(headers) + '\n')
        for row in transcript_description:
            transcript_id = row.split('\t')[0]
            transcript_to_details[transcript_id].append(row.strip().split('\t'))

        df_biocond = pd.read_csv(annot_filename, index_col=0, sep='\t')
        not_match = 0
        nb_start = 0
        nb_stop = 0
        nb_5utr = 0
        nb_3utr = 0
        nb_cds = 0
        nb_exon = 0
        nb_intron = 0
        set_type = set()
        for row_name, row in df_biocond.iterrows():
            if (isinstance(row['Gene_ID'], basestring)) and (row['Gene_ID'] != ''):
                transcript_id = row['Gene_ID']
                relat_pos_start = int(row['begin_peak'])
                relat_pos_end = int(row['end_peak'])
                strand = row['strand']
                chromo = row['chr']
                if transcript_id in transcript_to_details:
                    all_elements = transcript_to_details[transcript_id]
                    match_element = []
                    start_codon_present = False
                    utr_present = False
                    cds_present = False
                    exon_present = False

                    for element in all_elements:
                        begin_relat = int(element[2])
                        end_relat = int(element[3])
                        element_gene = sequence_python.SequencePython('element', begin_relat, end_relat, strand, 'mm10', chromo)
                        peak = sequence_python.SequencePython('peak', relat_pos_start, relat_pos_end, strand, 'mm10', chromo)
                        if element_gene.isoverlap(peak):
                            type_element = element[1]
                            peak_to_overlap[row_name].append(type_element)

    
    #print(peak_to_overlap)
    not_match = 0
    nb_5utr = 0
    nb_3utr = 0
    nb_cds = 0
    nb_intron = 0         
    with open(peaks_filename, 'w') as peaks_file, \
        open(annot_filename, 'rU') as annot_file:
        peaks_file.write(annot_file.readline())
        for row in annot_file:
            peak_id = row.split('\t')[0]
            overlaps = peak_to_overlap[peak_id]
            #print(overlaps)
            new_overlaps = ''
            if "3UTR" in overlaps:
                nb_3utr+=1
                new_overlaps += '3UTR' + ';'
            if "5UTR" in overlaps:
                nb_5utr+=1
                new_overlaps += '5UTR' + ';'
            if "CDS" in overlaps:
                nb_cds+=1
                new_overlaps += 'CDS' + ';'
            else:
                nb_intron+=1
                new_overlaps += 'Intron' + ';'
            peaks_file.write(row.replace('_Over_',new_overlaps))

    with open(relat_filename, 'w') as sum_overlap_file:
        new_row = 'Technique '+technique+'\n5\'UTR\t' + str(nb_5utr) + '\n3\'UTR\t' + str(nb_3utr) + \
            '\nCDS\t' + str(nb_cds) + '\nIntron \t' + str(nb_intron) + '\nNo Match \t' + str(not_match)
        #print(str(new_row))
        sum_overlap_file.write(str(new_row))



def search_ref_peak(path, exp_design_name, technique, bed_name):
    if technique == '':
        PATH_PEAKS = path + '/PeakDetection/Peaks'
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '_Annot.txt'
    else:
        PATH_PEAKS = path + '/PeakDetection/' + technique
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '_Annot.txt'

    ref_filename = path + '/Genome/' + 'pc_ep_mouse_mm10.bed'
    clip_filename = path + '/Genome/' + 'sb_m6a_mouse_mm10.bed'
    trew_filename = path + '/Genome/' + 'trew_mouse_mm10.bed'
    temp_PATH = path + '/PeakDetection/temp/'
    temp_ref_filename = temp_PATH + exp_design_name +'_' + technique + '_Peak_to_ref.txt'
    temp_clip_filename = temp_PATH + exp_design_name + '_' + technique + '_Peak_to_clip.txt'
    temp_trew_filename = temp_PATH + exp_design_name + '_' + technique + '_Peak_to_trew.txt'

    # prepare bed
    chromo_to_peak = defaultdict(list)
    peak_to_pos = dict()
    peak_to_ref = defaultdict(list)
    peak_to_clip = defaultdict(list)
    peak_to_trew = defaultdict(list)
    with open(peaks_filename, 'rU') as peak_file:
        csv_peaks = csv.DictReader(peak_file, delimiter = '\t')
        for row in csv_peaks:
            ref_gene = sequence_python.SequencePython('peak', int(row['begin_peak']), int(row['end_peak']), '+', 'mm10',\
                                                   row['chromo_peak'])
            peak_to_pos[row['WindowId']] = ref_gene
            chromo_to_peak[row['chromo_peak']].append(row['WindowId'])

    print('Search MeRIP-Seq')
    with open(ref_filename, 'rU') as ref_file:
        ref_file.readline()
        k = 0
        for element in ref_file:
            #print(element)
            begin = int(element.split('\t')[1])
            end = int(element.split('\t')[2])
            chromo = element.split('\t')[0]
            info = element.split('\t')[4]
            #if chromo == 'chr2':
            #    break
            ref_peak = sequence_python.SequencePython('element', begin, end, '+', 'mm10', chromo)
            k+=1
            if k%100000==0:
                print(k)
            if chromo_to_peak.has_key(chromo):
                for peak_id in chromo_to_peak.get(chromo):
                    peak = peak_to_pos[peak_id]
                    if (peak.begin < ref_peak.begin):
                        overlapLength = peak.end - ref_peak.begin;
                    else:
                        overlapLength = ref_peak.end - peak.begin
                    if (overlapLength > 0):
                        #print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
                        peak_to_ref[peak_id].append(info)

    with open(temp_ref_filename, 'w') as ref_peak_file:
        for peak_id, value in peak_to_ref.items():
            ref_peak_file.write(peak_id + '\t' + str(len(value))+'\t'+','.join(value)+'\n')

    print('Search CLIP-Seq')
    with open(clip_filename, 'rU') as ref_file:
        ref_file.readline()
        k = 0
        for element in ref_file:
            #print(element)
            begin = int(element.split('\t')[1])
            end = int(element.split('\t')[2])
            chromo = element.split('\t')[0]
            #if chromo == 'chr2':
            #    break
            ref_peak = sequence_python.SequencePython('element', begin, end, '+', 'mm10', chromo)
            k+=1
            if k%100000==0:
                print(k)
            if chromo_to_peak.has_key(chromo):
                for peak_id in chromo_to_peak.get(chromo):
                    peak = peak_to_pos[peak_id]
                    if (peak.begin < ref_peak.begin):
                        overlapLength = peak.end - ref_peak.begin;
                    else:
                        overlapLength = ref_peak.end - peak.begin
                    if (overlapLength > 0):
                        #print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
                        peak_to_clip[peak_id].append(1)

    with open(temp_clip_filename, 'w') as ref_peak_file:
        for peak_id, value in peak_to_clip.items():
            ref_peak_file.write(peak_id + '\t' + str(len(value))+'\n')

    print('Search TREW')
    with open(trew_filename, 'rU') as ref_file:
        ref_file.readline()
        k = 0
        for element in ref_file:
            #print(element)
            begin = int(element.split('\t')[1])
            end = int(element.split('\t')[2])
            chromo = element.split('\t')[0]
            info = element.split('\t')[4]
            ref_peak = sequence_python.SequencePython('element', begin, end, '+', 'mm10', chromo)
            k+=1
            #if chromo == 'chr2':
            #    break
            if k%100000==0:
                print(k)
            if chromo_to_peak.has_key(chromo):
                for peak_id in chromo_to_peak.get(chromo):
                    peak = peak_to_pos[peak_id]
                    if (peak.begin < ref_peak.begin):
                        overlapLength = peak.end - ref_peak.begin;
                    else:
                        overlapLength = ref_peak.end - peak.begin
                    if (overlapLength > 0):
                        #print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
                        peak_to_trew[peak_id].append(info)

    with open(temp_trew_filename, 'w') as ref_peak_file:
        for peak_id, value in peak_to_trew.items():
            ref_peak_file.write(peak_id + '\t' + str(len(value))+'\t'+','.join(value)+'\n')


def add_ref_peak(path, exp_design_name, technique, bed_name):
    if technique == '':
        PATH_PEAKS = path + '/PeakDetection/Peaks'
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '_Annot.txt'
        final_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '_Peaks.txt'
    else:
        PATH_PEAKS = path + '/PeakDetection/' + technique
        peaks_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '_Annot.txt'
        final_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '_Peaks.txt'

    ref_filename = path + '/Genome/' + 'pc_ep_mouse_mm10.bed'
    clip_filename = path + '/Genome/' + 'sb_m6a_mouse_mm10.bed'
    trew_filename = path + '/Genome/' + 'trew_mouse_mm10.bed'
    temp_PATH = path + '/PeakDetection/temp/'
    temp_ref_filename = temp_PATH + exp_design_name +'_' + technique + '_Peak_to_ref.txt'
    temp_clip_filename = temp_PATH + exp_design_name + '_' + technique + '_Peak_to_clip.txt'
    temp_trew_filename = temp_PATH + exp_design_name + '_' + technique + '_Peak_to_trew.txt'

    # prepare dict
    peak_to_ref = dict()
    peak_to_clip = dict()
    peak_to_trew = dict()
    with open(temp_ref_filename, 'rU') as ref_file:
        for row in ref_file:
            peak_to_ref[row.split('\t')[0]] = row.split('\t')[1].strip()
    with open(temp_clip_filename, 'rU') as ref_file:
        for row in ref_file:
            peak_to_clip[row.split('\t')[0]] = row.split('\t')[1].strip()
    with open(temp_trew_filename, 'rU') as ref_file:
        for row in ref_file:
            peak_to_trew[row.split('\t')[0]] = row.split('\t')[1].strip()

    nb_peak = 0
    nb_peak_meripseq = 0
    nb_peak_trew = 0
    nb_peak_clip = 0
    with open(peaks_filename, 'rU') as peak_file, \
          open(final_filename, 'w') as final_file:
        header = peak_file.readline().replace('Ref_peak','Ref_MeRIP\tCLIP\tTREW')
        final_file.write(header)
        for row in peak_file:
            nb_peak += 1
            peak_id = row.split('\t')[0]
            ref=''
            trew=''
            clip=''
            if peak_id in peak_to_ref:
                ref = peak_to_ref[peak_id]
                nb_peak_meripseq += 1
            if peak_id in peak_to_clip:
                clip = peak_to_clip[peak_id]
                nb_peak_clip += 1
            if peak_id in peak_to_trew:
                trew = peak_to_trew[peak_id]
                nb_peak_trew += 1

            new_row = row.replace('_Ref_',ref+'\t'+clip+'\t'+trew)
            final_file.write(new_row)

        print('Peaks ' + str(nb_peak) + ' MeRIP-Seq ' + str(nb_peak_meripseq) + ' CLIP-Seq ' +
            str(nb_peak_clip) + ' TREW ' + str(nb_peak_trew))

    os.remove(peaks_filename)


def create_peak_gtf(path, exp_design_name, technique, bed_name):
    """
    Read all PATH_PEAKS+'/'+exp_design_name+'_'+technique+'_'+Final.txt
    Combine peaks
    and save to GFF
    :param list_technique:
    :return:
    """
    PATH_ANNOT = path + '/Genome/'
    if technique == '' or technique == 'All':
        PATH_PEAKS = path + '/PeakDetection/Peaks'
        peak_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '_Peaks.txt'
        gtf_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '.gtf'
    else:
        PATH_PEAKS = path + '/PeakDetection/' + technique + '/'
        peak_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '_Peaks.txt'
        gtf_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '.gtf'

    with open(gtf_filename, 'w') as gtf_file, \
            open(peak_filename, 'rU') as peak_file:
        csv_peaks = csv.DictReader(peak_file, delimiter='\t')
        for row in csv_peaks:
            peak = HTSeq.GenomicInterval(row['chromo_peak'], int(row['begin_peak']), int(row['end_peak']), ".")
            peak_id = row['WindowId']
            feature = HTSeq.GenomicFeature(peak_id, 'exon', peak)
            #print(feature.get_gff_line().strip() + '; gene_id \"'+peak_id+'\"')
            gtf_file.write(feature.get_gff_line().strip() + '; gene_id \"'+peak_id+'\"'+'\n')
            


def create_peak_fasta(path, exp_design_name, mouse_seq, technique, bed_name):
    """
    Read all PATH_PEAKS+'/'+exp_design_name+'_'+technique+'_'+Final.txt
    Combine peaks
    and save to GFF
    :param list_technique:
    :return:
    """
    PATH_ANNOT = path + '/Genome/'
    if technique == '' or technique == 'All':
        PATH_PEAKS = path + '/PeakDetection/Peaks'
        peak_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '_Peaks.txt'
        fasta_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '.fasta'
    else:
        PATH_PEAKS = path + '/PeakDetection/' + technique + '/'
        peak_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '_Peaks.txt'
        fasta_filename = PATH_PEAKS + '/' + exp_design_name + '_' + technique + '_' + bed_name + '.fasta'

    with open(fasta_filename, 'w') as fasta_file, \
            open(peak_filename, 'rU') as peak_file:
        csv_peaks = csv.DictReader(peak_file, delimiter='\t')
        for row in csv_peaks:
            sequence_for_motif = mouse_seq[row['chromo_peak']][int(row['begin_peak']):int(row['end_peak'])].seq
            peak_id = row['WindowId']
            fasta_file.write('>' + peak_id + '\n' + str(sequence_for_motif) + '\n')
            


def read_mouse_seq(path, genome_file):
    with open(path + '/Genome/' + genome_file + '.fa', 'r') as mousefile:
        mouse_seq = SeqIO.to_dict(SeqIO.parse(mousefile, "fasta"))
        print('Mouse genome read')
        return mouse_seq


def read_gene_GTF(path, annotation_filename):
    """
    General function which read ref peaks file and return them
    :return: ref_peaks_list[feature.iv] = feature.get_gff_line().strip().replace('\t',';')
    """
    print(path + '/Genome/' + annotation_filename)
    gtf_peaks_file = HTSeq.GFF_Reader(path + '/Genome/' + annotation_filename )
    gene_list = dict()
    # load GTF file with all transcripts
    for feature in gtf_peaks_file:
            #print(feature.attr['transcript_id'])
            gene_list[feature.attr['gene_id']] = feature
    print('Gene GTF - ' + annotation_filename + str(len(gene_list)))
    return gene_list


def main(argv):
    '''
    Main function of peak_calling_fisher.py
    :param argv:
    -p --path Path of the working folder
    -a --annotation Annotation file in GTF format
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:e:t:a:g:b:", ["path=", "expdesign=","peak_technique=","annotation=","genome=",'bed_name='])
    except getopt.GetoptError:
        print('Cannot run command - Help: peak_annotation.py -p <path> -e <expdesign> -t <peak_technique> '
              '-a <annotation> -g <genome> -b <bed_name>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('peak_annotation.py -p <path> -e <expdesign> -t <peak_technique> -a <annotation> -g <genome> -b <bed_name>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-e", "--expdesign"):
            exp_design_name = arg
        elif opt in ("-t", "--peaktechnique"):
            peak_technique = arg
        elif opt in ("-a", "--annotation"):
            annotation_file = arg
        elif opt in ("-g", "--genome"):
            genome_file = arg
        elif opt in ("-b", "--bed_name"):
            bed_name = arg
    # Need to replace All per '' when looking at whole peaks from different techniques
    if peak_technique == 'All':
        peak_technique = ''
    
    # choose gene or transcript
    annotation_filename = annotation_file + '.annotation.gene'
    #annotation_filename = annotation_file + '.annotation.transcript.txt'
    
    # read mouse genome and list of genes
    mouse_seq = read_mouse_seq(path, genome_file)
    gene_list = read_gene_GTF(path, annotation_filename + '.gtf')
    annotate_peaks(path, exp_design_name, mouse_seq, gene_list, annotation_filename + '.txt', peak_technique, bed_name)
    create_peak_gtf(path, exp_design_name, peak_technique, bed_name)
    create_peak_fasta(path, exp_design_name, mouse_seq, peak_technique, bed_name)
    peak_relative_pos(path, exp_design_name, peak_technique, annotation_filename + '.description.txt', bed_name)
    print('Search MeRIP-Seq, CLIP-Seq, and TREW from MeT-DB V2')
    search_ref_peak(path, exp_design_name, peak_technique, bed_name)
    print('Add MeRIP-Seq, CLIP-Seq, and TREW from MeT-DB V2')
    add_ref_peak(path, exp_design_name, peak_technique, bed_name)
    
if __name__ == "__main__":
    main(sys.argv[1:])
