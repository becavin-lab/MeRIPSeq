#!/usr/bin/env python

import csv
import os
import glob
import HTSeq
#import ConfigParser

PATH = '/Users/cbecavin/Documents/m6aAkker/'
PATH_Proteomics = PATH + 'Proteomics/'
PATH_ANNOT = PATH + 'Genome/'
LIST_UNIPROT = 'listes_Uniprot'
LIST_GENE_NAME = 'listes_GeneName'


def parseSARResults():
    
    ## listFolder=['AllIP','AllInput','IPZT3Seq2','IPZT3Seq345','IPZT11Seq345','IPZT11Seq1','IPZT11Seq2','IPSeq1','IPSeq2','IPSeq345','InputSeq1','InputSeq2','InputSeq345','InputZT11Seq1','InputZT11Seq2','InputZT11Seq345','InputZT3Seq2','InputZT3Seq345']
    ##listFolder=['InputSeq2345','IPSeq2345','InputSeq2345Corrected','IPSeq2345Corrected']
    listFolder = ['InputSeq789','IPSeq789', 'InputAll_Liver','IPAll_Liver']
    PATH='/Users/cbecavin/Documents/m6aAkker'
    for folder in listFolder:
        print(folder)
        pathName=PATH+'/Analysis_July_2017/'+folder+'/tables/'
        pathNameNew = PATH + '/Analysis_Sept_2017/' + folder + '/tables/'
        pathNameList = PATH+'/Analysis_Sept_2017/'+folder+'/'+LIST_UNIPROT+'/'
        pathNameListGene = PATH + '/Analysis_Sept_2017/' + folder + '/'+LIST_GENE_NAME+'/'
        try:
            os.mkdir(pathNameNew)
        except OSError:
            print('Already exists')
        try:
            os.mkdir(pathNameList)
        except OSError:
            print('Already exists')
        try:
            os.mkdir(pathNameListGene)
        except OSError:
            print('Already exists')
        print(pathName)

        #
        # Read uniprot vs GeneName
        # and Uniprot vs HallMark
        uniprot_to_genename = dict()
        uniprot_to_hallmark = dict()
        with open(PATH_ANNOT + "gencodeVM13/gencode.vM13.annotation.entrez.uniprot.txt", "r") as annotation_file, \
                open(PATH_ANNOT + 'VM9_Uniprot_to_geneName.txt', 'r') as missing_annot_file:
            missing_annot_file.readline()
            # because VM9 was used for RNASeq analysis there 1000 uniprot IDs missing in VM13 annotation
            # So we add them here
            for row in missing_annot_file:
                uniprot_to_genename[row.split('\t')[0]] = row.split('\t')[1]

            annotation_csv = csv.DictReader(annotation_file, delimiter='\t')
            for row in annotation_csv:
                uniprot_IDs = row.get('UniprotIDs')
                if uniprot_IDs is not None:
                    uniprot_IDs = uniprot_IDs.split(';')
                    gene_name = row.get('gene_name')
                    hallmark = row.get('HallMarkH_GeneSets')
                    for uniprot in uniprot_IDs:
                        uniprot_to_genename[uniprot] = gene_name
                        uniprot_to_hallmark[uniprot] = hallmark

        listFiles = glob.glob(pathName+'*.complete.txt')
        print(listFiles)
        for file in listFiles:
            print(file)
            fileUp = file.split('.')[0]+'.up.txt'
            fileUp_new = (file.split('.')[0] + '.up.txt').replace(pathName, pathNameNew)
            listfileUp = (file.split('.')[0]+'_'+folder+'.up.txt').replace(pathName, pathNameList)
            listfileUp_gene = (file.split('.')[0] + '_' + folder + '.up.txt').replace(pathName, pathNameListGene)
            allgenes=[]
            with open(fileUp, "r") as deSeqResult, \
                    open(fileUp_new, "w") as deSeq_new, \
                    open(listfileUp, "w") as listdeSeqResult, \
                    open(listfileUp_gene, "w") as listdeSeqResult_gene:
                dictreader = csv.DictReader(deSeqResult, delimiter='\t')

                headers = ''
                for header in dictreader.fieldnames:
                    headers += header + '\t'
                deSeq_new.write(headers + 'Gene_Name\tHallMark\n')

                print(fileUp,listfileUp)
                for rowProtein in dictreader:
                    gene = rowProtein.get('Id')
                    listdeSeqResult.write(gene+'\n')
                    allgenes.append(gene)
                    gene_name = ''
                    if gene in uniprot_to_genename.keys():
                        gene_name = uniprot_to_genename[gene]
                        listdeSeqResult_gene.write((gene_name+'\n'))
                    hallmark = ''
                    if gene in uniprot_to_hallmark.keys():
                        hallmark = uniprot_to_hallmark[gene]

                    new_row = ''
                    for header in dictreader.fieldnames:
                        new_row += rowProtein.get(header) + '\t'
                    new_row += gene_name + '\t' + hallmark
                    deSeq_new.write(new_row + '\n')
             
            fileDown=file.split('.')[0]+'.down.txt'
            fileDown_new = (file.split('.')[0] + '.down.txt').replace(pathName,pathNameNew)
            listfileDown=(file.split('.')[0]+'_'+folder+'.down.txt').replace(pathName, pathNameList)
            listfileDown_gene = (file.split('.')[0] + '_' + folder + '.down.txt').replace(pathName, pathNameListGene)
            with open(fileDown, "r") as deSeqResult, \
                    open(fileDown_new, "w") as deSeq_new, \
                    open(listfileDown_gene, "w") as listdeSeqResult_gene, \
                    open(listfileDown, "w") as listdeSeqResult:
                dictreader = csv.DictReader(deSeqResult, delimiter='\t')

                headers = ''
                for header in dictreader.fieldnames:
                    headers += header + '\t'
                deSeq_new.write(headers + 'Gene_Name\tHallMark\n')

                for rowProtein in dictreader:
                    gene = rowProtein.get('Id')
                    listdeSeqResult.write(gene+'\n')
                    allgenes.append(gene)
                    gene_name = ''
                    if gene in uniprot_to_genename.keys():
                        gene_name = uniprot_to_genename[gene]
                        listdeSeqResult_gene.write((gene_name+'\n'))
                    hallmark = ''
                    if gene in uniprot_to_hallmark.keys():
                        hallmark = uniprot_to_hallmark[gene]

                    new_row = ''
                    for header in dictreader.fieldnames:
                        new_row += rowProtein.get(header) + '\t'
                    new_row += gene_name + '\t' + hallmark
                    #print(new_row)
                    deSeq_new.write(new_row+'\n')

            listfileAll = (file.split('.')[0]+'_'+folder+'.all.txt').replace(pathName, pathNameList)
            listfileAll_gene = (file.split('.')[0] + '_' + folder + '.all.txt').replace(pathName, pathNameListGene)

            with open(listfileAll, "w") as allGenesFile, \
                    open(listfileAll_gene, "w") as allGenesFile_gene:
                for gene in allgenes:
                    allGenesFile.write(gene+'\n')
                    gene_name = ''
                    if gene in uniprot_to_genename.keys():
                        gene_name = uniprot_to_genename[gene]
                        allGenesFile_gene.write((gene_name+'\n'))

            fileComplete = file.split('.')[0] + '.complete.txt'
            fileComplete_new = (file.split('.')[0] + '.complete.txt').replace(pathName, pathNameNew)
            with open(fileComplete, "r") as deSeqResult, \
                    open(fileComplete_new, "w") as deSeqResult_gene:
                dictreader = csv.DictReader(deSeqResult, delimiter='\t')

                headers = ''
                for header in dictreader.fieldnames:
                    headers += header + '\t'
                deSeqResult_gene.write(headers + 'Gene_Name\tHallMark\n')

                for rowProtein in dictreader:
                    gene = rowProtein.get('Id')
                    gene_name = ''
                    if gene in uniprot_to_genename.keys():
                        gene_name = uniprot_to_genename[gene]
                    hallmark = ''
                    if gene in uniprot_to_hallmark.keys():
                        hallmark = uniprot_to_hallmark[gene]
                    new_row = ''
                    for header in dictreader.fieldnames:
                        new_row += rowProtein.get(header) + '\t'
                    new_row += gene_name + '\t' + hallmark
                    deSeqResult_gene.write(new_row + '\n')


def parse_transcriptomic_table():
    with open(PATH + 'Genome/BioMart_Mouse_annot_uniprot_06-2016.txt', 'r') as genomics_file:
        uniprot_to_allassociated_uniprot = dict()
        dictreader = csv.DictReader(genomics_file, delimiter='\t')
        for row_genome in dictreader:
            uniprot = row_genome.get('UniprotID')
            uniprot_all = row_genome.get('All_UniprotID')
            print(uniprot,uniprot_all)
            if uniprot is not None:
                uniprot_to_allassociated_uniprot[uniprot] = uniprot_all
    with open(PATH_Proteomics + 'ConvMousevsGermFree.complete.txt', 'r') as transcriptomics_file, \
            open(PATH_Proteomics + 'Transcriptomics_CM_vs_GF.txt', 'w') as transcriptomics_final_file:
        # update transcriptomics table by search All_uniprotID
        transcriptomics_final_file.write('All_UniprotID\t'+transcriptomics_file.readline())
        i=0
        for row in transcriptomics_file:
            uniprot = row.split('\t')[0]
            if uniprot in uniprot_to_allassociated_uniprot.keys():
                uniprot_all = uniprot_to_allassociated_uniprot.get(uniprot)
                print(uniprot_all)
                if uniprot_all == None:
                    uniprot_all = ''
            #print(gene_name)
            i+=1
            transcriptomics_final_file.write(uniprot_all+'\t'+row)
        print(i)


def correlate_omics():
    uniprot_to_annot = dict()
    with open(PATH + 'Genome/BioMart_Mouse_annot_uniprot_06-2016.txt', 'r') as genomics_file:
        dictreader = csv.DictReader(genomics_file, delimiter='\t')
        header_annotation = dictreader.fieldnames
        for row_genome in dictreader:
            uniprot = row_genome.get('All_UniprotID')
            if uniprot is not None:
                uniprot_all_ID = uniprot.split(';')
                for uniprot_id in uniprot_all_ID:
                    annot = []
                    for keys in dictreader.fieldnames:
                        annot.append(row_genome.get(keys))
                    uniprot_to_annot[uniprot_id] = '\t'.join(annot)


    with open(PATH_Proteomics + 'Transcriptomics_CM_vs_GF.txt', 'rU') as transcriptomics_file, \
            open(PATH_Proteomics + 'Proteomics_CM_vs_GF.txt', 'rU') as proteomics_file, \
            open(PATH_Proteomics + 'Trans_Proteo_CM_GF.txt', 'w') as omics_file:
        transcriptomics_logfc = dict()
        dictreader_transcriptomics = csv.DictReader(transcriptomics_file, delimiter='\t')
        dictreader_proteomics = csv.DictReader(proteomics_file, delimiter='\t')
        for row in dictreader_transcriptomics:
            uniprot = row.get('All_UniprotID')
            uniprot_all_ID = uniprot.split(';')
            #print(uniprot_all_ID)
            if not uniprot == '':
                for uniprot_id in uniprot_all_ID:
                    log2 = row.get('log2FoldChange')
                    padj = row.get('padj')
                    if log2 == 'NA':
                        log2 = 0.0
                    if padj == 'NA':
                        padj = 1.0
                    log2 = float(log2)
                    padj = float(padj)
                    transcriptomics_logfc[uniprot_id] = [log2,padj]
        omics_file.write('Uniprot\tTranscriptomics logFC\tTranscript p-value\tProteomics logFC\tProt p-value\t'
                         'Omics_comparison\t'+'\t'.join(header_annotation)+'\n')
        for row in dictreader_proteomics:
            uniprot = row.get('UniProt')
            annot = uniprot_to_annot.get(uniprot)
            if annot is None:
                annot=''
            log2_prot = row.get('log2')
            padj_prot = row.get('Adjusted p.value')
            if padj_prot == 'NA':
                padj_prot = 0.0
            elif padj_prot == '':
                padj_prot = 1.0
            if uniprot in transcriptomics_logfc.keys():
                log2_trans = transcriptomics_logfc.get(uniprot)[0]
                padj_trans = transcriptomics_logfc.get(uniprot)[1]
            else:
                log2_trans = 0.0
                padj_trans = 1.0

            log2_prot = float(log2_prot)
            padj_prot = float(padj_prot)

            omics_value = 0
            if padj_prot < 0.05 and padj_trans < 0.05:
                if log2_trans > 0 and log2_prot > 0:
                    omics_value = 2
                if log2_trans < 0 and log2_prot < 0:
                    omics_value = -2
                if log2_trans > 0 and log2_prot < 0:
                    omics_value = 1
                if log2_trans < 0 and log2_prot > 0:
                    omics_value = -1
            results = [str(log2_trans),str(padj_trans),str(log2_prot),str(padj_prot),str(omics_value)]
            omics_file.write(uniprot+'\t'+'\t'.join(results)+'\t'+annot+'\n')


def add_m6a():
    with open(PATH + 'Genome/BioMart_Mouse_annot_uniprot_06-2016.txt', 'r') as genomics_file:
        uniprot_to_allassociated_uniprot = dict()
        dictreader = csv.DictReader(genomics_file, delimiter='\t')
        for row_genome in dictreader:
            uniprot = row_genome.get('UniprotID')
            uniprot_all = row_genome.get('All_UniprotID')
            #print(uniprot,uniprot_all)
            if uniprot is not None:
                uniprot_to_allassociated_uniprot[uniprot] = uniprot_all


    with open(PATH_Proteomics + 'Trans_Proteo_CM_GF.txt', 'rU') as trans_proteo_file, \
            open(PATH_Proteomics + 'm6aExpDesign_CM_GF_Seq8_Filtered_Annot.txt', 'rU') as m6a_file, \
            open(PATH_Proteomics + 'Omics_CM_GF.txt', 'w') as omics_file:
        trans_proteo_dict = csv.DictReader(trans_proteo_file, delimiter='\t')
        m6a_dict = csv.DictReader(m6a_file, delimiter='\t')
        print(trans_proteo_dict.fieldnames)

        headers_to_keep = ['Peak_ID','Peak_start','Peak_end','Relat_pos','Motif_presence','Ref_Peak_nb',
                           'GermFree_Seq8_Med','ConvMouse_Seq8_Med']


        headers = headers_to_keep + ['m6a_FoldChange'] + trans_proteo_dict.fieldnames[1:]
        uniprot_to_omics = dict()
        for row_trans_prot in trans_proteo_dict:
            uniprot_id = row_trans_prot.get('Uniprot')
            all_row = []
            for header in trans_proteo_dict.fieldnames[1:]:
                value = row_trans_prot.get(header)
                if value is None:
                    value=''
                all_row.append(value)
            uniprot_to_omics[uniprot_id] = '\t'.join(all_row)

        omics_file.write('\t'.join(headers)+'\n')
        print(headers)
        for row in m6a_dict:
            uniprot = row.get('UniProt / SwissProt Accession')
            print(uniprot)
            if uniprot is not None:
                if uniprot in uniprot_to_allassociated_uniprot.keys():
                    all_uniprot_id = uniprot_to_allassociated_uniprot[uniprot]
                    all_row = []
                    for header in headers_to_keep:
                        value = row.get(header)
                        if value is None:
                            value = ''
                        all_row.append(value)

                    conv_mouse = float(row.get('ConvMouse_Seq8_Med'))
                    germ_free = float(row.get('GermFree_Seq8_Med'))
                    fold_change = 0.0
                    if conv_mouse == 0:
                        fold_change = -germ_free
                    elif germ_free == 0:
                        fold_change = conv_mouse
                    else:
                        fold_change = conv_mouse/germ_free
                    found = False
                    new_row = '\t'.join(all_row) + '\t' + str(fold_change)
                    for uniprot_id in all_uniprot_id.strip(';'):
                        if uniprot in uniprot_to_omics.keys():
                            new_row += '\t' + uniprot_to_omics[uniprot]
                            omics_file.write(new_row+'\n')
                            found = True
                            break
                    #if not found:
                    #    omics_file.write(new_row + '\n')


# extract list of genes from SARTools results
parseSARResults()

# Proteomics vs Transcriptomics
#parse_transcriptomic_table()
#correlate_omics()
#add_m6a()

