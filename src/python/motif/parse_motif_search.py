#!/usr/bin/env python
from collections import defaultdict
import os
from shutil import copyfile
import os.path
import sys, getopt

SOFTWARES = {'meme','dreme','centrimo','meme_tomtom'}
##TYPES_SEARCHES = {'All','Narrow','Const','Max'}



def parse_all_logs(path_analysis, exp_design_name):
    with open(path_analysis + 'Motif_search.log', 'w') as general_log:
        header = 'Data_name\t' + '\t'.join(SOFTWARES)
        general_log.write(header+'\n')
        for data_name in os.listdir(path_analysis + 'results/'):
            if not data_name.startswith('.') and not data_name.endswith('.log') \
                    and not data_name.endswith('.sh'): # and '_MaxMaxValues' in data_name:
                progress_filename = path_analysis + 'results/' + data_name + '/progress_log.txt'
                line = data_name + '\t'
                for software in SOFTWARES:
                    found = False
                    with open(progress_filename,'r') as file:
                        for row in file:
                            status_soft = 'name: '+software+'  status: 0'
                            if status_soft in row:
                                found = True
                    if found:
                        line += '1\t'
                    else:
                        line += '0\t'
                general_log.write(line.strip() + '\n')


def parse_dreme(path_analysis, path_motif, exp_design_name):
    with open(path_analysis + 'dreme.sh', 'w') as dreme_file:
        for data_name in os.listdir(path_analysis + 'results/'):
            #if exp_design_name in data_name:
            print(data_name)
            if not data_name.startswith('.') and not data_name.endswith('.log') \
                    and not data_name.endswith('.sh'): # and '_MaxMaxValues' in data_name:
                # print(data_name)
                dreme_xml = path_analysis + 'results/' + data_name + '/dreme_out/dreme.xml'
                if os.path.exists(dreme_xml):
                    dreme_list = list()
                    with open(dreme_xml, 'r') as file:
                        for row in file:
                            dreme_list.append(row.strip())
                            if '<command_line>' in row:
                                dreme_file.write(row)

                    dreme_txt = path_analysis + 'results/' + data_name + '/dreme_out/dreme.txt'
                    dreme_txt_list = list()
                    with open(dreme_txt, 'r') as file:
                        for row in file:
                            dreme_txt_list.append(row.strip())
                    dreme_txt_intro = ''
                    write_intro = False
                    for row in dreme_txt_list:
                        if row.startswith('MOTIF '):
                            break
                        if 'MEME version' in row:
                            write_intro = True
                        if write_intro:
                            dreme_txt_intro += row + '\n'

                    #print(dreme_txt_intro)

                    # parse dreme
                    with open(path_motif + '/logs/Dreme_' + data_name + '.log', 'w') as general_log:
                        header = 'ID\tMotif\tNb_sites\tPos_Occurence\tNeg_Occurence\tNb_Seq\t' \
                                 'Percent\tPos_Percent\tNeg_Percent\tEvalue'
                        general_log.write(header + '\n')
                        # search length fasta
                        for line in dreme_list:
                            if '<positives name=' in line:
                                number_seq = float(line.split(' ')[2].replace('\"','').replace('count=',''))
                        # search motifs
                        for i in range(1, len(dreme_list)):
                            line = dreme_list[i]
                            if '<motifs>' in line:
                                index_motifs = i
                                break

                        # get motifs info
                        for k in range(index_motifs+1, len(dreme_list)):
                            if '<motif' in dreme_list[k]:
                                new_line = dreme_list[k].replace('\"','').split(' ')
                                #print(new_line)
                                id = new_line[1].replace('id=', '')
                                #print(id,data_name)
                                print(dreme_txt_intro)
                                seq = new_line[2].replace('seq=', '')
                                occ = float(new_line[4].replace('nsites=', ''))
                                pos_occ = float(new_line[5].replace('p=', ''))
                                neg_occ = float(new_line[6].replace('n=', ''))
                                evalue = new_line[8].replace('evalue=', '')

                                # save table
                                line_to_write = [id, seq, str(occ), str(pos_occ),str(neg_occ), str(number_seq),
                                                 str(occ/number_seq), str(pos_occ/number_seq), str(neg_occ/number_seq),
                                                 str(evalue)]
                                general_log.write('\t'.join(line_to_write) + '\n')

                                # copy png file
                                for imagefile in os.listdir(path_analysis + 'results/' + data_name + '/dreme_out/'):
                                    if id in imagefile:
                                        if 'nc_' in imagefile:
                                            src = path_analysis + 'results/' + data_name + '/dreme_out/' + imagefile
                                            dst = path_motif + '/motif_figure/'+seq + '_nc.png'
                                            copyfile(src, dst)
                                        else:
                                            src = path_analysis + 'results/' + data_name + '/dreme_out/' + imagefile
                                            dst = path_motif + '/motif_figure/' + seq + '_rc.png'
                                            copyfile(src, dst)

                                # save motif
                                with open(path_motif + '/motif/' + seq + '.meme', 'w') as meme_file:
                                    meme_file.write(dreme_txt_intro+'\n')
                                    # find motif info
                                    dreme_txt_motif = ''
                                    write_motif = False
                                    for row in dreme_txt_list:
                                        if write_motif and row.startswith('MOTIF '):
                                            break
                                        if ('MOTIF '+seq+' DREME') in row:
                                            write_motif = True
                                        if write_motif:
                                            dreme_txt_motif += row + '\n'

                                    print(dreme_txt_motif)
                                    meme_file.write(dreme_txt_motif + '\n')


def regroup_motifs(path_motif, exp_design_name, pvalue_cutoff):
    motif_set = set()
    dict_dataset = defaultdict(list)
    dict_evalue = defaultdict(list)
    #data_selection = 'NoabxGF.'
    for data_name in os.listdir(path_motif + 'logs/'):
        if data_name.startswith('Dreme_') and data_name.endswith('.log') and exp_design_name in data_name:
            type_data = data_name.replace('Dreme_','').replace('.log','')
            if 'Meth' in type_data:
                #if data_selection in data_name:
                print(type_data, data_name)
                with open(path_motif + 'logs/' + data_name,'r') as log_file:
                    log_file.readline()
                    for line in log_file:
                        motif = line.split('\t')[1]
                        evalue = line.strip().split('\t')[-1]
                        if float(evalue) < pvalue_cutoff:
                            motif_set.add(motif)
                            print(evalue,motif)
                            dict_dataset[motif].append(type_data)
                            dict_evalue[motif].append(evalue)
            else:
                print(type_data, data_name)
                with open(path_motif + 'logs/' + data_name, 'r') as log_file:
                    log_file.readline()
                    for line in log_file:
                        motif = line.split('\t')[1]
                        evalue = line.strip().split('\t')[-1]
                        if float(evalue) < pvalue_cutoff:
                            motif_set.add(motif)
                            print(evalue, motif)
                            dict_dataset[motif].append(type_data)
                            dict_evalue[motif].append(evalue)


    with open(path_motif + 'Motif_'+exp_design_name+ '__raw.txt','w') as motif_log:
         motif_log.write("Motif\tData\tP-value\n")
         for motif in motif_set:
             print(motif+'\t'+';'.join(dict_dataset[motif]))
             print(motif+'\t'+';'.join(dict_dataset[motif])+'\t'+';'.join(dict_evalue[motif])+'\n')
             motif_log.write(motif+'\t'+';'.join(dict_dataset[motif])+'\t'+';'.join(dict_evalue[motif])+'\n')


def regroup_motifs_diff(path_motif, exp_design_name, pvalue_cutoff):
    '''
    Regroup all motifs coming from comparisons of biological conditions
    :param path_motif:
    :param exp_design_name:
    :param pvalue_cutoff:
    :return:
    '''
    motif_set = set()
    dict_dataset = defaultdict(list)
    dict_evalue = defaultdict(list)
    for data_name in os.listdir(path_motif + 'logs/'):
        if data_name.startswith('Dreme_') and data_name.endswith('.log') and '_MaxMaxValues' not in data_name:
            type_data = data_name.replace('Dreme_','').replace('.log','')
            print(type_data, data_name)
            with open(path_motif + 'logs/' + data_name,'r') as log_file:
                log_file.readline()
                for line in log_file:
                    motif = line.split('\t')[1]
                    evalue = line.strip().split('\t')[-1]
                    if float(evalue) < pvalue_cutoff:
                        motif_set.add(motif)
                        print(evalue,motif)
                        dict_dataset[motif].append(type_data)
                        dict_evalue[motif].append(evalue)

    with open(path_motif + 'Motif_'+exp_design_name+'_raw.txt','w') as motif_log:
        motif_log.write("Motif\tData\tP-value\n")
        for motif in motif_set:
            motif_log.write(motif+'\t'+';'.join(dict_dataset[motif])+'\t'+';'.join(dict_evalue[motif])+'\n')


def select_centrally_enriched_diff(path_analysis, path_motif, exp_design_name):
    '''
    Add a column indicating if motif are centrally enriched, and select only motif coming form a comparisons.
    Meaning : exp_design_name not in D
    :param path_analysis:
    :param path_motif:
    :param exp_design_name:
    :return:
    '''
    centrimo_motif = dict()
    for data_name in os.listdir(path_motif + 'logs/'):
        if data_name.startswith('Dreme_') and data_name.endswith('.log') and exp_design_name in data_name:
            type_data = data_name.replace('Dreme_','').replace('.log','')
            centrimo_path = path_analysis + '/results/' + type_data + '/centrimo_out/centrimo.txt'
            if os.path.exists(centrimo_path):
                with open(centrimo_path, 'r') as centrimo_file:
                    for line in centrimo_file:
                        if 'DREME' in line:
                            motif = line.split('\t')[1].split(' ')[0]
                            print(motif, type_data)
                            centrimo_motif[motif] = type_data

    print(len(centrimo_motif))

    with open(path_motif + 'Motif_' + exp_design_name + '_raw.txt', 'r') as motif_log, \
            open(path_motif + 'Motif_' + exp_design_name + '_Diff.txt', 'w') as motif_filter_log:
        motif_filter_log.write("Motif\tData\tP-value\tCentrallyEnriched\n")
        motif_log.readline()
        for line in motif_log:
            print(line.split('\t')[0])
            found = False
            datas = line.split('\t')[1].split(';')
            for data in datas:
                print(data)
            if line.split('\t')[0] in centrimo_motif:
                motif_filter_log.write(line.strip()+'\tYes\n')
            else:
                motif_filter_log.write(line.strip() + '\tNo\n')


def select_centrally_enriched(path_analysis, path_motif, exp_design_name):
    #data_selection = 'NoLP.'

    centrimo_motif = dict()
    for data_name in os.listdir(path_motif + 'logs/'):
        if data_name.startswith('Dreme_') and data_name.endswith('.log') and exp_design_name in data_name:
            type_data = data_name.replace('Dreme_','').replace('.log','')
            centrimo_path = path_analysis + '/results/' + type_data + '/centrimo_out/centrimo.txt'
            if os.path.exists(centrimo_path):
                with open(centrimo_path, 'r') as centrimo_file:
                    for line in centrimo_file:
                        if 'DREME' in line:
                            motif = line.split('\t')[1].split(' ')[0]
                            print(motif, type_data)
                            centrimo_motif[motif] = type_data

    print(len(centrimo_motif))

    with open(path_motif + 'Motif_' + exp_design_name + '_raw.txt', 'r') as motif_log, \
            open(path_motif + 'Motif_' + exp_design_name + '.txt', 'w') as motif_filter_log:
        motif_filter_log.write("Motif\tData\tP-value\tCentrallyEnriched\n")
        motif_log.readline()
        for line in motif_log:
            print(line.split('\t')[0])
            if line.split('\t')[0] in centrimo_motif:
                motif_filter_log.write(line.strip()+'\tYes\n')


def create_motif(path_motif, motif):
    motif_file = path_motif + '/motif/' + motif + '.meme'
    motif_figure_file = path_motif + '/motif_figure/' + motif + '_nc.eps'
#    motif_figure_file = path_motif + '/motif_figure/' + motif + '_nc.png'
    # Meme-suite as to be installed for using = iupac2meme
    iupac2_command = 'iupac2meme -dna ' + motif + ' > ' + motif_file
    print(iupac2_command)
    os.system(iupac2_command)
    ceqlogo_command = 'ceqlogo -i ' + motif_file + ' -m ' + motif + ' -o ' + motif_figure_file + ' -f eps'
#    ceqlogo_command = 'ceqlogo -i ' + motif_file + ' -m ' + motif + ' -o '+motif_figure_file+' -f png'
    print('Run ' + ceqlogo_command)
    os.system(ceqlogo_command)


def main(argv):
    '''
    Main function of parse_motif_search.py
    :param argv:
    -p --path Path of the working folder
    -e --expdesign name of the exp_design
    -b --bed_name suffix of the bed file
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv, "hp:e:b:",
                                   ["path=", "expdesign=", 'bed_name='])
    except getopt.GetoptError:
        print(
            'Cannot run command - Help: parse_motif_search.py -p <path> -e <expdesign> -b <bed_name>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('peak_finalize.py -p <path> -e <expdesign> -b <bed_name>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-e", "--expdesign"):
            exp_design_name = arg
        elif opt in ("-b", "--bed_name"):
            bed_name = arg

    path_analysis = path + 'PeakDiffExpression/' + exp_design_name + '_' + bed_name + '/Motif/'
    path_motif = path + 'PeakDiffExpression/Motif/'

    # Parse meme results and extract all motifs
    parse_dreme(path_analysis, path_motif, exp_design_name + '_' + bed_name)
    # regroup all motifs
    regroup_motifs(path_motif, exp_design_name + '_' + bed_name)
    pvalue_cutoff = 0.00001
    #regroup_motifs(path_motif, exp_design_name, pvalue_cutoff)
    # Select motif centrally enriched -> Could be done automatically
    #select_centrally_enriched(path_analysis, path_motif, exp_design_name)
    # regroup all motifs for comparisons
    # Remove all logs for each exp_design_name
    pvalue_cutoff = 0.01
    #regroup_motifs_diff(path_motif, exp_design_name, pvalue_cutoff)
    # Select motif centrally enriched -> Could be done automatically
    #select_centrally_enriched_diff(path_analysis, path_motif, exp_design_name)

    # Manually filter motif list and put it in path_motif
    # add new motif
    #motif = 'RRACH'
    #motif="NGGACN"
    #create_motif(path_motif, motif)
    #motif = 'NBCAN'
    #add_motif(path_motif, motif, motif_list_filename)


if __name__ == "__main__":
    main(sys.argv[1:])
