import HTSeq
import sys, getopt
from Bio import SeqIO


def find_peak_max(path, exp_design_name, size, mouse_seq, peak_technique, bed_name):
    final_bed_name = 'MaxValues'
    if bed_name == 'Raw':
        final_bed_name = 'MaxValues'
    if 'MaxValues' in bed_name:
        final_bed_name = bed_name.replace('MaxValues','MaxMaxValues')
    print(final_bed_name)
    if peak_technique == '':
        PATH_PEAKS = path + '/PeakDetection/Peaks'
        peak_init_bed_filename = PATH_PEAKS + '/' + exp_design_name + '_' + bed_name + '.bed'
        peak_txt_filename = PATH_PEAKS + '/' +  exp_design_name + '_' + final_bed_name + '.txt'
        peak_bed_filename = PATH_PEAKS + '/' +  exp_design_name + '_' + final_bed_name + '.bed'
    else:
        PATH_PEAKS = path + '/PeakDetection/' + peak_technique + '/'
        peak_init_bed_filename = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_' + bed_name + '.bed'
        peak_txt_filename = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_' + final_bed_name + '.txt'
        peak_bed_filename = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_' + final_bed_name + '.bed'

    # load wig file
    cvg = HTSeq.GenomicArray("auto", stranded=False, typecode='i')
    wig_filename = path + '/Mapping/' + exp_design_name + '_IP_Mean_clean.wig'
    k = 0
    print(wig_filename)
    with open(wig_filename, 'rU') as wig_file:
        #for i in range(1,10000000):
        for row in wig_file:
    #       row = wig_file.readline()
            k += 1
            if k % 1000000 == 0:
                print(k)
            #print(row)
            chromo = row.split('\t')[0]
            start = int(row.split('\t')[1])
            coverage = int(row.split('\t')[2])
            interval = HTSeq.GenomicInterval(chromo, start, start+1, "+")
            if (start > 1): # this is to avoid a mysterious segmentation fault when the wig data for a given chromosome starts at position 0,1,2,3
                cvg[interval] = coverage


    # load bed file with all peaks
    with open(peak_init_bed_filename, 'r') as bed_file:
        peaks = dict()
        i=0
        for row in bed_file:
            chromo = row.split('\t')[0]
            begin = int(row.split('\t')[1])
            end = int(row.split('\t')[2])
            peak_id = 'Peak_' + str(i)
            i+=1
            strand = row.split('\t')[5].strip()
            iv = HTSeq.GenomicInterval(chromo, begin, end ,strand)
            #print(feature.iv)
            peaks[peak_id] = iv

    print('Read window GTF ' + str(len(peaks)) + ' loaded')

    with open(peak_bed_filename, "w") as bed_peaks_file, \
            open(peak_txt_filename, "w") as txt_peaks_file:
        txt_peaks_file.write('PeakID\tMaxIndex\n')
        for name, peakiv in peaks.items():
            # Calculate median coverage for each window
            if len(name) != 0:
                max_value = 0
                max_index = peakiv.start + peakiv.length/2
                for iv, value in cvg[peakiv].steps():
                    # lengs.append(iv.length)
                    if value > max_value:
                        max_index = iv.start
                        max_value = value
                    # print(windowiv, medianvalue, lengs)
                txt_peaks_file.write(name+'\t'+str(max_index)+'\n')
                #sequence = mouse_seq[iv.chrom][max_index-size/2:max_index+size/2].seq
                peakiv.start = max_index - int(size)/2
                peakiv.end = max_index + int(size)/2
                #new_line = [iv.chrom, '.', 'exon', str(peakiv.start), str(peakiv.end), '.', '.', '.',
                #            'ID \"'+name+'\"; gene_id \"'+name+'\"']
                #gtf_peaks_file.write('\t'.join(new_line)+'\n')
                if peak_technique == '':
                    info = '1'
                else:
                    info = peak_technique
                bed_peaks_file.write(peakiv.chrom + '\t' + str(peakiv.start) + '\t' + str(peakiv.end) +
                               '\t' + name + '\t'+info+'\t+\n')


def read_mouse_seq(path, genome_file):
    with open(path + '/Genome/' + genome_file + '.fa', 'r') as mousefile:
        mouse_seq = SeqIO.to_dict(SeqIO.parse(mousefile, "fasta"))
        print('Mouse genome read')
        return mouse_seq


def main(argv):
    '''
    Main function of find_peak_center.py
    :param argv:
    -p --path Path of the working folder
    -a --annotation Annotation file in GTF format
    -p --peak_technique PPeak technique used
    -s --size size of the final peak
    -g --genome genome name
    :return:
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:e:t:s:g:b:", ["path=", "expdesign=","peak_technique=","size=","genome=","bed_name"])
    except getopt.GetoptError:
        print('Cannot run command - Help: find_peak_max.py -p <path> -e <expdesign> -t <peak_technique> '
              '-s <size> -g <genome> -b <bed_name>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('find_peak_max.py -p <path> -e <expdesign> -t <peak_technique> -s <size> -g <genome> -b <bed_name>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-e", "--expdesign"):
            exp_design_name = arg
        elif opt in ("-t", "--peaktechnique"):
            peak_technique = arg
        elif opt in ("-s", "--size"):
            size = arg
        elif opt in ("-g", "--genome"):
            genome_file = arg
        elif opt in ("-b", "--bed_name"):
            bed_name = arg

    if peak_technique == 'All':
        peak_technique = ''
    mouse_seq = read_mouse_seq(path, genome_file)
    find_peak_max(path, exp_design_name, size, mouse_seq, peak_technique, bed_name)


if __name__ == "__main__":
    main(sys.argv[1:])
