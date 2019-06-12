import pytz
#from utils import m6a_utils
#from utils import exp_design
import m6a_utils #a
import m6a_utils as utils #a
import subprocess
import csv
from collections import defaultdict
import collections
import numpy
import HTSeq
import pandas as pd
import math
import os
#from utils import sequence_python
import sequence_python #a


def create_venn_diagram(PATH_PEAKS, name, list_technique):
    list_peak_list = set()
    for technique in list_technique:
        print(technique)
        with open(PATH_PEAKS + name + '.bed', 'rU') as all_peak, \
                open(PATH_PEAKS + name + '_' + technique +'_List.txt', 'w') as list_peak:
            index = 0
            for row in all_peak:
                peak_id = 'Peak_'+str(index)
                index += 1
                #print(row)
                if technique in row.split('\t')[4]:
                    list_peak.write(peak_id+'\n')
                    if technique == 'Cecum':
                        id = row.split('\t')[5]
                        yo = id.split(',')
                        for temp in yo:
                            if 'Cecum_Peak' in temp:
                                list_peak_list.add(peak_id)
                                print(row)
                    #print(peak_id)

    print(len(list_peak_list))
    with open(PATH_PEAKS + name + '_Peaks_List_all.txt', 'w') as list_peak:
        for peak_id in list_peak_list:
            list_peak.write(peak_id+'\n')



def get_ref_peak():
    #path = '/Users/cbecavin/Documents/m6aAkker/Genome/'
    path = '/pasteur/projets/policy01/m6aAkker/Genome/' #a
    # modify sb file
    awk_command = 'awk \'BEGIN{FS=OFS="\\t"}!/^,/ && NR>1{print $2,$3,$4,"CLIP",$5,$7}\' ' \
                  ''+path+'/sb_m6a_mouse_mm10 > '+path+'/sb_m6a_mouse_mm10.bed'
    os.system(awk_command)
    # modify TREW data
    awk_command = 'awk \'BEGIN{FS=OFS="\\t"}NR>1{print $2,$3,$4,"TREW",$5,$6}\' ' \
                  ''+path+'/trew_mouse_mm10 > '+path+'/trew_mouse_mm10.bed'
    os.system(awk_command)
    # modify MeRIPSeq
    awk_command = 'awk \'BEGIN{FS=OFS="\\t"}NR>1{print $2,$3,$4,"MeRIP",$5,$7}\' ' \
                  ''+path+'/pc_ep_mouse_mm10 > '+path+'/pc_ep_mouse_mm10.bed'
    os.system(awk_command)

    list_file = ['sb_m6a_mouse_mm10','trew_mouse_mm10','pc_ep_mouse_mm10']
    for file_name in list_file:
         sort_command = 'sort -k1,1 -k2,2n ' + path + file_name+'.bed > ' + path + file_name+'_sort.bed'
         os.system(sort_command)
         merge_command = 'bedtools merge -i ' + path +file_name+'_sort.bed -c 4,5,6 -o distinct > ' + \
                         path+file_name+'_merge.bed'
         os.system(merge_command)

    wc_command = 'wc -l '+path+'/*.bed'
    os.system(wc_command)
    #print(wc_command)


def merge_tissue_peaks():
    #cecum_filename = m6a_utils.PATH_PEAKS + 'Peaks/CecAm_Peaks_Ref.txt'
    cecum_filename = m6a_utils.PATH_PEAKS + 'Peaks/Cecum_all_MaxMaxValues_3_Peaks.txt' #a
    #liver_filename = m6a_utils.PATH_PEAKS + 'Peaks/LivOld_Peaks_Ref.txt'
    liver_filename = m6a_utils.PATH_PEAKS + 'Peaks/Liver_all_MaxMaxValues_3_Peaks.txt' #a
    
    ref_filename = m6a_utils.PATH_PEAKS + 'Peaks/Ref_Peaks.txt'
    final_filename = m6a_utils.PATH_PEAKS + 'Peaks/Liver_Cecum_Peaks.txt'
    with open(final_filename, 'w') as final_file:
        df_caecum = pd.read_csv(cecum_filename, index_col=0, sep='\t')
        df_liver = pd.read_csv(liver_filename, index_col=0, sep='\t')
        df_ref = pd.read_csv(ref_filename, index_col=0, sep='\t')
        liver_peaks = dict()
        ref_peaks = dict()
        cecum_peaks = dict()
        peak_cecum_number = 0
        peak_ref_number = 0
        peak_liver_number = 0


        # Load peaks
        chrom_to_peaks_liver = defaultdict(list)
        chrom_to_peaks_ref = defaultdict(list)
        chrom_to_peaks_cecum = defaultdict(list)
        for peak_id_liver, row_liver in df_liver.iterrows():
            peak_liver_number += 1
            liver_peak = sequence_python.SequencePython('Peak', row_liver['begin_peak'], row_liver['end_peak'], '+',
                                                            'mm10',row_liver['chr'])
            liver_peaks[peak_id_liver] = liver_peak
            chrom_to_peaks_liver[row_liver['chr']].append(peak_id_liver)
        for peak_id_ref, row_ref in df_ref.iterrows():
            peak_ref_number += 1
            ref_peak = sequence_python.SequencePython('Peak', row_ref['begin_peak'], row_ref['end_peak'], '+',
                                                            'mm10',row_ref['chr'])
            ref_peaks[peak_id_ref] = ref_peak
            chrom_to_peaks_ref[row_ref['chr']].append(peak_id_ref)
        for peak_id_cecum, row_cecum in df_caecum.iterrows():
            peak_cecum_number += 1
            cecum_peak = sequence_python.SequencePython('Peak', row_cecum['begin_peak'], row_cecum['end_peak'], '+',
                                                            'mm10',row_cecum['chr'])
            cecum_peaks[peak_id_cecum] = cecum_peak
            chrom_to_peaks_cecum[row_cecum['chr']].append(peak_id_cecum)

        # create lists

        # liver_only = set()
        # ref_only() = set()
        # cecum_only() = set()
        #
        # liver_cecum = set()
        # liver_ref  = set()
        # cecum_ref = set()
        #
        # liver_cecum_ref = set()

        # cecum only
        nb_peak = 0
        index = 0

        chrom_to_2 = chrom_to_peaks_cecum
        chrom_to_3 = chrom_to_peaks_liver
        dict_peaks_1 = ref_peaks
        dict_peaks_2 = cecum_peaks
        dict_peaks_3 = liver_peaks

        print(len(dict_peaks_1))
        print(len(dict_peaks_2))
        print(len(dict_peaks_3))

        for peak_id, peak_pos in dict_peaks_1.items():
            index += 1
            if index%1000==0:
                print(index,peak_id)
            found_2 = False
            for peak_id_2 in chrom_to_2[peak_pos.chromosome]:
                peak_pos_2 = dict_peaks_2[peak_id_2]
                if (peak_pos_2.begin < peak_pos.begin):
                    overlapLength = peak_pos_2.end - peak_pos.begin;
                else:
                    overlapLength = peak_pos.end - peak_pos_2.begin
                if (overlapLength > 0):
                    # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end
                    found_2 = True
                    break
            found_3 = False
            for peak_id_3 in chrom_to_3[peak_pos.chromosome]:
                peak_pos_3 = dict_peaks_3[peak_id_3]
                if (peak_pos_3.begin < peak_pos.begin):
                    overlapLength = peak_pos_3.end - peak_pos.begin;
                else:
                    overlapLength = peak_pos.end - peak_pos_3.begin
                if (overlapLength > 0):
                    # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end
                    found_3 = True
                    break
            #if (found_3):
            if (not found_3) and (found_2):
                nb_peak += 1

        print(nb_peak)

        print(len(dict_peaks_1))
        print(len(dict_peaks_2))
        print(len(dict_peaks_3))

        # overlap_liver = set()
        # overlap_cecum_liver = 0
        # overlap_liver_cecum = 0
        # overlap_ref_cecum = 0
        # overlap_ref_cecum_liver = 0
        # overlap_ref_liver = 0
        # final_file.write('Peak_Id\t'+'\t'.join(df_caecum.columns)+'\tLiver_peak\tNb_Liver_Peak\tgeneral_id\n')
        # index_general_id = 1
        # for peak_id_caecum, row in df_caecum.iterrows():
        #     #print(peak_id_caecum)
        #     if not '_NC_' in peak_id_caecum:
        #         caecum_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10', row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_liver[row['chr']]:
        #             liver_peak = liver_peaks[peak_id_liver]
        #             if (liver_peak.begin < caecum_peak.begin):
        #                 overlapLength = liver_peak.end - caecum_peak.begin;
        #             else:
        #                 overlapLength = caecum_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_liver)
        #                 overlap_liver.add(peak_id_liver)
        #
        #         print_row = ''
        #         # for colname in df_caecum.columns:
        #         #     print_row += str(row[colname])+'\t'
        #         # general_id = 'Peak_' + str(index_general_id)
        #         # new_row = 'Caecum_'+peak_id_caecum+'\t'+print_row+';'.join(peak_overlap)+'\t'+str(len(peak_overlap))+'\t'+general_id+'\n'
        #         if len(peak_overlap) != 0:
        #             overlap_liver_cecum += 1
        #             #print('no overlap')
        #         #final_file.write(new_row)
        #         index_general_id += 1
        #
        # for peak_id_caecum, row in df_liver.iterrows():
        #     #print(peak_id_caecum)
        #     if not '_NC_' in peak_id_caecum:
        #         caecum_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10', row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_cecum[row['chr']]:
        #             liver_peak = cecum_peaks[peak_id_liver]
        #             if (liver_peak.begin < caecum_peak.begin):
        #                 overlapLength = liver_peak.end - caecum_peak.begin;
        #             else:
        #                 overlapLength = caecum_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_liver)
        #                 overlap_liver.add(peak_id_liver)
        #
        #         print_row = ''
        #         # for colname in df_caecum.columns:
        #         #     print_row += str(row[colname])+'\t'
        #         # general_id = 'Peak_' + str(index_general_id)
        #         # new_row = 'Caecum_'+peak_id_caecum+'\t'+print_row+';'.join(peak_overlap)+'\t'+str(len(peak_overlap))+'\t'+general_id+'\n'
        #         if len(peak_overlap) != 0:
        #             overlap_cecum_liver += 1
        #             #print('no overlap')
        #         #final_file.write(new_row)
        #         index_general_id += 1
        #
        #
        # for peak_id_caecum, row in df_caecum.iterrows():
        #     #print(peak_id_caecum)
        #     if not '_NC_' in peak_id_caecum:
        #         caecum_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10', row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_ref[row['chr']]:
        #             liver_peak = ref_peaks[peak_id_liver]
        #             if (liver_peak.begin < caecum_peak.begin):
        #                 overlapLength = liver_peak.end - caecum_peak.begin;
        #             else:
        #                 overlapLength = caecum_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_liver)
        #                 overlap_liver.add(peak_id_liver)
        #
        #         if len(peak_overlap) != 0:
        #             overlap_ref_cecum += 1
        #
        #
        # for peak_id_caecum, row in df_liver.iterrows():
        #     #print(peak_id_caecum)
        #     if not '_NC_' in peak_id_caecum:
        #         caecum_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10', row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_ref[row['chr']]:
        #             liver_peak = ref_peaks[peak_id_liver]
        #             if (liver_peak.begin < caecum_peak.begin):
        #                 overlapLength = liver_peak.end - caecum_peak.begin;
        #             else:
        #                 overlapLength = caecum_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_liver)
        #                 overlap_liver.add(peak_id_liver)
        #
        #         if len(peak_overlap) != 0:
        #             overlap_ref_liver += 1
        #
        # for peak_id_caecum, row in df_liver.iterrows():
        #     #print(peak_id_caecum)
        #     if not '_NC_' in peak_id_caecum:
        #         caecum_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10', row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_ref[row['chr']]:
        #             liver_peak = ref_peaks[peak_id_liver]
        #             if (liver_peak.begin < caecum_peak.begin):
        #                 overlapLength = liver_peak.end - caecum_peak.begin;
        #             else:
        #                 overlapLength = caecum_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 for peak_id_liver in chrom_to_peaks_cecum[row['chr']]:
        #                     liver_peak = cecum_peaks[peak_id_liver]
        #                     if (liver_peak.begin < caecum_peak.begin):
        #                         overlapLength = liver_peak.end - caecum_peak.begin;
        #                     else:
        #                         overlapLength = caecum_peak.end - liver_peak.begin
        #                     if (overlapLength > 0):
        #                         # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                         peak_overlap.append(peak_id_liver)
        #                         overlap_liver.add(peak_id_liver)
        #
        #         if len(peak_overlap) != 0:
        #             overlap_ref_cecum_liver += 1


        # for peak_id_ref, row in df_ref.iterrows():
        #     print('Ref',peak_id_ref)
        #     if not '_NC_' in peak_id_ref:
        #         peak_ref_number += 1
        #         ref_peak = sequence_python.SequencePython('Peak', row['begin_peak'], row['end_peak'], '+', 'mm10',
        #                                                      row['chr'])
        #         peak_overlap = list()
        #         for peak_id_liver in chrom_to_peaks_liver[row['chr']]:
        #             liver_peak = liver_peaks[peak_id_liver]
        #             if (liver_peak.begin < ref_peak.begin):
        #                 overlapLength = liver_peak.end - ref_peak.begin;
        #             else:
        #                 overlapLength = ref_peak.end - liver_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_liver)
        #         if len(peak_overlap) != 0:
        #             overlap_ref_liver += 1
        #             # print('no overlap')
        #
        #         peak_overlap = list()
        #         for peak_id_cecum in chrom_to_peaks_cecum[row['chr']]:
        #             cecum_peak = cecum_peaks[peak_id_cecum]
        #             if (cecum_peak.begin < ref_peak.begin):
        #                 overlapLength = cecum_peak.end - ref_peak.begin;
        #             else:
        #                 overlapLength = ref_peak.end - cecum_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 peak_overlap.append(peak_id_cecum)
        #         if len(peak_overlap) != 0:
        #             overlap_ref_cecum += 1
        #
        #         peak_overlap = list()
        #         for peak_id_cecum in chrom_to_peaks_cecum[row['chr']]:
        #             cecum_peak = cecum_peaks[peak_id_cecum]
        #             if (cecum_peak.begin < ref_peak.begin):
        #                 overlapLength = cecum_peak.end - ref_peak.begin;
        #             else:
        #                 overlapLength = ref_peak.end - cecum_peak.begin
        #             if (overlapLength > 0):
        #                 # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                 for peak_id_liver in chrom_to_peaks_liver[row['chr']]:
        #                     liver_peak = liver_peaks[peak_id_liver]
        #                     if (liver_peak.begin < ref_peak.begin):
        #                         overlapLength = liver_peak.end - ref_peak.begin;
        #                     else:
        #                         overlapLength = ref_peak.end - liver_peak.begin
        #                     if (overlapLength > 0):
        #                         # print(peak.begin,peak.end,ref_peak.begin,ref_peak.end)
        #                         peak_overlap.append(peak_id_liver)
        #         if len(peak_overlap) != 0:
        #             overlap_ref_cecum_liver += 1
        # for peak_id_liver, row_liver in df_liver.iterrows():
        #     peak_liver_number += 1
        #     if peak_id_liver != overlap_liver:
        #         print_row = ''
        #         for colname in df_liver.columns:
        #             print_row += str(row_liver[colname]) + '\t'
        #         general_id = 'Peak_' + str(index_general_id)
        #         new_row = 'Liver_' + peak_id_liver + '\t' + print_row + '\t\t\t'+general_id+'\n'
        #         final_file.write(new_row)
        #         index_general_id += 1
        #
        # print(general_id,index_general_id)
        # print('Overlap Caecum + Liver', overlap_liver_cecum, peak_cecum_number, float(100*overlap_liver_cecum/peak_cecum_number))
        # print('Overlap Liver + Caecum', overlap_cecum_liver, peak_liver_number, float(100*len(overlap_liver)/peak_liver_number))
        # print('Overlap Ref + Liver', overlap_ref_liver, peak_liver_number, float(100 * overlap_ref_liver / peak_liver_number))
        # print('Overlap Ref + Cecum', overlap_ref_cecum, peak_cecum_number, float(100 * overlap_ref_cecum / peak_cecum_number))
        # print('Overlap Ref (Liver - Cecum)', overlap_ref_cecum_liver, overlap_ref_liver, overlap_ref_cecum, peak_ref_number, float(100 * overlap_ref_liver / peak_ref_number), float(100 * overlap_ref_cecum / peak_ref_number))


def compare_techniques(exp_design_name):
    peaks_table = m6a_utils.PATH_PEAKS + 'Peaks/' + exp_design_name + '_Peaks_Ref.txt'
    df_peaks = pd.read_csv(peaks_table, index_col=0, sep='\t')
    file_macs2 = m6a_utils.PATH_PEAKS + 'Peaks/' + exp_design_name + '_MACS2_List.txt'
    file_fisher = m6a_utils.PATH_PEAKS + 'Peaks/' + exp_design_name + '_Fisher_List.txt'
    file_poi = m6a_utils.PATH_PEAKS + 'Peaks/' + exp_design_name + '_POI_List.txt'
    file_rpmf = m6a_utils.PATH_PEAKS + 'Peaks/' + exp_design_name + '_RPMF_List.txt'
    with open(file_macs2, 'w') as file_macs2, \
            open(file_poi, 'w') as file_poi,\
            open(file_fisher, 'w') as file_fisher, \
            open(file_rpmf, 'w') as file_rpmf:
        for peak_id, row_liver in df_peaks.iterrows():
            technique = row_liver['Peak_Presence']
            if "MACS2" in technique:
                file_macs2.write(peak_id+'\n')
            if "Fisher" in technique:
                file_fisher.write(peak_id+'\n')
            if "POI" in technique:
                file_poi.write(peak_id+'\n')
            if "RPMF" in technique:
                file_rpmf.write(peak_id+'\n')


################################################################
################################################################
#   RUN ALL SCRIPTS HERE

PATH_PEAKS = m6a_utils.PATH_PEAKS + '/Peaks/'

#list_technique = ['Cecum_Old', 'Cecum']
#name = 'Cecum_Old_vs_New'

#list_technique = ['CecAm', 'LivOld','Ref']
#name = 'CecAm_LivOld_Ref'

#list_technique = ['LivOld', 'Liver', 'LivPaper']
#name = 'LiverOld_vs_New'

list_technique = ['Cecum_all', 'Liver_all', 'Ref']
name = 'CecumAll_vs_LiverAll_Ref'

#list_technique = ['Cecum_Old', 'Liver_Old', 'Ref']
#name = 'Cecum_vs_Liver_Ref_Old'

#list_technique = ['Cecum', 'Liver', 'Ref']
#name = 'Cecum_vs_Liver_Ref'

#get_ref_peak()
#
# cat_command = 'cat '
# for technique in list_technique:
#     awk_command = 'awk \'BEGIN{OFS="\t"}NR>1{print $2,$3,$4,"+","'+technique+'",$1}\' '+PATH_PEAKS+'/'+ \
#                   technique+'_Peaks_Ref.txt > ' + PATH_PEAKS+'/'+ \
#                   technique+'.bed'
#     print(awk_command)
#     os.system(awk_command)
# for technique in list_technique:
#     cat_command += PATH_PEAKS+'/' + technique+'.bed '
# cat_command += ' > '+PATH_PEAKS+'/' + name+'_raw.bed'
# print(cat_command)
# os.system(cat_command)
# sort_command = 'sort -k1,1 -k2,2n ' +PATH_PEAKS+'/' + name+'_raw.bed > ' + \
#                PATH_PEAKS + '/' + name + '_sort.bed'
# print sort_command
# os.system(sort_command)
# merge_command = 'bedtools merge -i '+PATH_PEAKS+'/' + name+\
#                 '_sort.bed -c 4,5,6 -o distinct > '+PATH_PEAKS+'/' + name+'.bed'
# print(merge_command)
# os.system(merge_command)
# create_venn_diagram(PATH_PEAKS, name, list_technique)

merge_tissue_peaks()


#compare_techniques('LivOld')
