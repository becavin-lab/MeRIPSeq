#!/usr/bin/env python
import HTSeq
import csv
from collections import defaultdict
import sys, getopt


GENE_TYPE = ['protein_coding']
#GENE_TYPE = ['protein_coding','processed_pseudogene','lincRNA','TEC','unprocessed_pseudogene']


def create_gene_and_exon_gtf(path_annot, annotation_file):
    '''
    Keep only 'gene' in the GTF for count and median expression calculation in the Peak detection workflow
    :return:
    '''
    print('Create '+annotation_file + '.gene.gtf file for Peak detection')
    gtf_file = HTSeq.GFF_Reader(path_annot + '/' + annotation_file + '.annotation.gtf', end_included=True)
    NUMBER_ELEMENT_GTF = 0
    with open(path_annot + '/' + annotation_file + '.annotation.gene.gtf', "w") as cleanGTF:
        number_gene = 0
        for feature in gtf_file:
            NUMBER_ELEMENT_GTF += 1
            if feature.type == 'gene':
                if feature.attr['gene_type'] in GENE_TYPE:
                    #print(feature.get_gff_line())
                    number_gene+=1
                    cleanGTF.write(feature.get_gff_line())
        print('Gone through ' + str(NUMBER_ELEMENT_GTF) + ' of GTF file')
        print('Found '+str(number_gene)+' genes')

    print('Create ' + annotation_file + '.exon.gtf file for Peak detection')
    gtf_file = HTSeq.GFF_Reader(path_annot + '/' + annotation_file + '.annotation.gtf', end_included=True)
    with open(path_annot + '/' + annotation_file + '.annotation.exon.gtf', "w") as cleanGTF:
        number_exon = 0
        for feature in gtf_file:
            if feature.type == 'gene':
                if feature.attr['gene_type'] in GENE_TYPE:
                    #print(feature.get_gff_line())
                    cleanGTF.write(feature.get_gff_line())
            if feature.type == 'exon':
                if feature.attr['gene_type'] in GENE_TYPE:
                    #print(feature.get_gff_line())
                    number_exon +=1
                    cleanGTF.write(feature.get_gff_line())
        print('Found ' + str(number_exon) + ' exons')
    return NUMBER_ELEMENT_GTF

# From DAVID website I created a list of all ENTREZID to UniprotID
# I use ti to add UniprotID in the BioMart annotation table
def parse_GTF_transcript_annotation(path_annot,annotation_file):
    '''
    WARNING : Need to create gencode.vM13.DAVID.entrezID.to.uniprotID.txt
    :return:
    '''
    print('Create ' + annotation_file + '.annotation.txt file for Peak detection')
    with open(path_annot + '/' + annotation_file + '.metadata.EntrezGene', "rU") as file_entrez, \
            open(path_annot + '/' + annotation_file + '.DAVID.entrezID.to.uniprotID.txt', "rU") as file_uniprot_ID, \
            open(path_annot + '/' + annotation_file + '.annotation.gtf', "rU") as file_transcript, \
            open(path_annot + '/' + annotation_file + '.annotation.txt', "w") as new_annotation_file:
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

                #if number_gene % 10000 == 0:
                #    print(feature.get_gff_line().strip())
        print(number_gene)


def create_gene_annotation(path_annot, annotation_file):
    '''
    Create an annotation file with only the gene
    :return:
    '''
    print('Create ' + annotation_file + '.annotation.gene.txt file for Peak detection')
    with open(path_annot + '/' + annotation_file + '.annotation.txt', "rU") as file_annotation, \
            open(path_annot + '/' + annotation_file + '.annotation.gene.txt', "w") as gene_annotation:
        gene_file = set()
        csv_annot = csv.DictReader(file_annotation, delimiter='\t')
        headers = ['gene_id', 'chr', 'gene_begin', 'gene_end', 'strand', 'type', 'gene_name',
                   'EntrezID', 'UniprotID', 'UniprotIDs', 'Description']
        gene_annotation.write('\t'.join(headers)+'\n')
        for row_annot in csv_annot:
            row = ''
            if row_annot['type'] in GENE_TYPE:
                for header in headers:
                    #print(header,row_annot)
                    if type(row_annot[header]) == str:
                        row += row_annot[header] + '\t'
                gene_file.add(row)

        for row in gene_file:
            gene_annotation.write(row.strip() + '\n')


def parse_gene_description(path_annot, annotation_file, NUMBER_ELEMENT_GTF):
    print('Create ' + annotation_file + '.annotation.gene.description.txt file for Peak detection')
    print('Position of UTR, exon, and intron in each gene')
    with open(path_annot + '/' + annotation_file + '.annotation.gtf','r') as annotation_gtf, \
            open(path_annot + '/' + annotation_file + '.annotation.gene.description.txt', 'w') as gene_description:
        gtf_peaks_file = HTSeq.GFF_Reader(annotation_gtf)
        gene_to_details = defaultdict(list)
        k=0
        for feature in gtf_peaks_file:
            if 'gene_id' in feature.attr.keys():
                if feature.attr['gene_type'] in GENE_TYPE:
                    gene_id = feature.attr['gene_id']
                    #print(feature.get_gff_line())
                    gene_to_details[gene_id].append(feature)
                    k+=1
                    if k%50000 == 0:
                        print(annotation_file + '.annotation.gtf: line '+ str(k) + '/' + str(NUMBER_ELEMENT_GTF))

        no_utr = 0
        number_exon = 0
        number_3utr = 0
        number_5utr = 0
        for gene_id in gene_to_details:
            #print(transcr_id+'\t'+str(transcript_to_details[transcr_id]))
            list_feature = gene_to_details[gene_id]
            list_utr = []

            for feature in list_feature:
                if feature.type == 'gene':
                    gene = feature
                elif feature.type == 'UTR':
                    list_utr.append(feature)

            # get gene
            if len(list_utr) == 0:
                #print('no UTR in gene '+transcr_id)
                no_utr += 1
                #for feature in list_feature:
                #    print(transcript.name, transcript.attr['gene_type'], feature.type)
            else:
                for feature in list_feature:
                    type = 'exon'
                    diffStart = float(feature.iv.start - gene.iv.start)
                    if gene.iv.strand == '-':
                        diffStart = float(gene.iv.end - feature.iv.end)

                    relative_pos_begin = diffStart / gene.iv.length
                    type = feature.type
                    if feature.type == 'UTR':
                        if relative_pos_begin < 0.5:
                            type = '5UTR'
                            number_5utr += 1
                        else:
                            type = '3UTR'
                            number_3utr += 1
                    else:
                        number_exon += 1
                    gene_description.write(gene_id+'\t'+type+'\t'+str(feature.iv.start)+'\t'+str(feature.iv.end)+'\n')
                    #print(type, transcr_id,feature.iv.strand,feature.type, relative_pos_begin,relative_pos_end)
        print('Number of exon: ' + str(number_exon))
        print('Number of 5\'UTR: ' + str(number_5utr))
        print('Number of 3\'UTR: ' + str(number_3utr))
        print('Number of genes without UTR: ' + str(no_utr) + '/' + str(len(gene_to_details)))


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
        opts, args = getopt.getopt(argv,"hp:a:",["path=","annotation="])
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

    NUMBER_ELEMENT_GTF = create_gene_and_exon_gtf(path_annot, annotation_file)
    #NUMBER_ELEMENT_GTF = 1693722

    parse_GTF_transcript_annotation(path_annot, annotation_file)
    create_gene_annotation(path_annot, annotation_file)
    parse_gene_description(path_annot, annotation_file, NUMBER_ELEMENT_GTF)

if __name__ == "__main__":
    main(sys.argv[1:])


