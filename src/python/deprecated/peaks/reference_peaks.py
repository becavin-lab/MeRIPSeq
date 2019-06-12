#!/usr/bin/env python

import csv
import HTSeq
from Bio import SeqIO

PATH = '/Users/cbecavin/Documents/m6aAkker/'

def get_peak_data():
    # modify sb file
    awk_command = 'awk \'BEGIN{FS=OFS="\\t"}!/^,/{print $2,$3,$4,"CLIP",$5,$7,$11}\' sb_m6a_mouse_mm10 > sb_m6a_mouse_mm10.txt'
    print(awk_command)
    awk_command = 'awk \'BEGIN{FS=OFS="\\t"}!/^,/ && NR>1{print $2,$3,$4,"CLIP",$5,$7}\' sb_m6a_mouse_mm10 > sb_m6a_mouse_mm10.bed'
    print(awk_command)
    # modify TREW data
    awk_command = 'awk \'BEGIN{FS=OFS="\\t"}{print $2,$3,$4,"TREW",$5,$6,$7,$13,$14,$15}\' trew_mouse_mm10 > trew_mouse_mm10.txt'
    print(awk_command)
    awk_command = 'awk \'BEGIN{FS=OFS="\\t"}NR>1{print $2,$3,$4,"TREW",$5,$6}\' trew_mouse_mm10 > trew_mouse_mm10.bed'
    print(awk_command)
    awk_command = 'awk \'BEGIN{FS=OFS="\\t"}{print $2,$3,$4,"MeRIP",$5,$7}\' pc_ep_mouse_mm10 > pc_ep_mouse_mm10.txt'
    print(awk_command)
    awk_command = 'awk \'BEGIN{FS=OFS="\\t"}NR>1{print $2,$3,$4,"MeRIP",$5,$7}\' pc_ep_mouse_mm10 > pc_ep_mouse_mm10.bed'
    print(awk_command)

    list_file = ['sb_m6a_mouse_mm10','trew_mouse_mm10','pc_ep_mouse_mm10']
    for file_name in list_file:
        sort_command = 'sort -k1,1 -k2,2n ' + PATH + 'ReferencePeaks/'+file_name+'.bed > ' + PATH + 'ReferencePeaks/'+ \
                       file_name+'_sort.bed'
        print sort_command
        merge_command = 'bedtools merge -i ' + PATH + 'ReferencePeaks/'+file_name+'_sort.bed -c 4,5,6 -o distinct > ' + \
                        PATH + 'ReferencePeaks/'+file_name+'_merge.bed'
        print(merge_command)

    wc_command = 'wc -l '+PATH+'ReferencePeaks/*.bed'
    print(wc_command)

def prepare_peak_ref_file():
    '''
    From CSV files create bed files
    :return:
    '''
    with open(PATH + 'ReferencePeaks/Listdata.txt', 'rU') as txt_file:
        for row in txt_file:
            data_name = row.split('\t')[0]
            new_data_name = row.split('\t')[1]
            with open(PATH + 'ReferencePeaks/'+data_name+'.csv', 'rU') as txt_peaks_file, \
                open(PATH + 'ReferencePeaks/' + new_data_name + '.bed', 'w') as bed_peaks_file:
                list_peaks = csv.reader(txt_peaks_file, delimiter=',')
                for elements in list_peaks:
                    chromo = elements[2]
                    start = elements[3]
                    end = elements[4]
                    enrichment = elements[5]
                    if 'Macs' in data_name:
                        strand = '+'
                        type = 'macs2'
                    else:
                        strand = elements[6]
                        type = 'exomePeak'
                    if elements[0] != 'Methylation Type':
                        new_row = [chromo,start,end,strand,enrichment,type,new_data_name]
                        bed_peaks_file.write('\t'.join(new_row) + '\n')


def regroup_peaks():
    cat_command = 'cat '+PATH+'ReferencePeaks/trew_mouse_mm10.bed '+PATH+'ReferencePeaks/sb_m6a_mouse_mm10.bed '+\
                  PATH+'ReferencePeaks/pc_ep_mouse_mm10.bed > All_Peaks_raw.bed'
    print(cat_command)
    sort_command = 'sort -k1,1 -k2,2n ' + PATH + 'ReferencePeaks/All_Peaks_raw.bed > ' + PATH + 'ReferencePeaks/All_Peaks_sort.bed'
    print sort_command
    merge_command = 'bedtools merge -i '+PATH+'ReferencePeaks/All_Peaks_sort.bed -c 4,5,6 -o distinct > '+PATH + 'ReferencePeaks/All_Peaks.bed'
    print(merge_command)
    awk_command = 'awk \'{OFS = "\\t"; print $1,$2,$3,"Ref","1","+","Ref"}\' '+PATH + 'ReferencePeaks/All_Peaks.bed > '+PATH + 'ReferencePeaks/All_Peaks_List.bed'
    print(awk_command)


##############################@
#      From MeT-DB   V2.0
#  http://180.208.58.19/metdb_v2/html/download.php
#  download :
# MERIP-Seq6 m6a sites : pc_ep_mouse_mm9
# http://180.208.58.19/metdb_v2/download/overall/pc_ep_mouse_mm9.gz
# CLIP-Seq singl base m6a site : sb_m6a_mouse_mm10
# http://180.208.58.19/metdb_v2/download/overall/sb_m6a_mouse_mm10.gz
# m6a modification target protein site : trew_mouse_mm10
# http://180.208.58.19/metdb_v2/download/overall/trew_mouse_mm10.gz

get_peak_data()
#prepare_peak_ref_file()
#regroup_peaks()
