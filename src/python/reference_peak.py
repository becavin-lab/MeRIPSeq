#!/usr/bin/env python
import os, sys, getopt


def get_peak_data(path):
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
        print sort_command
        merge_command = 'bedtools merge -i ' + path + file_name+'_sort.bed -c 4,5,6 -o distinct > ' + \
                        path +file_name+'_merge.bed'
        print(merge_command)

    wc_command = 'wc -l '+path+'/*.bed'
    os.system(wc_command)
    

def main(argv):
    '''
    Main function of reference_peak.py
    :param argv:
    -p --path Path of the working folder
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:",["path="])
    except getopt.GetoptError:
        print('Cannot run command - Help: reference_peak.py -p <path>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('reference_peak.py -p <path>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path_annot = arg

        get_peak_data(path_annot)

if __name__ == "__main__":
    main(sys.argv[1:])


