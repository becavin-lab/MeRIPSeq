import HTSeq
import sys, getopt
from Bio import SeqIO


MOTIF_METH = {'GGACC': 6, 'CCACC': 6, 'CGACC': 6, 'GCACC': 6, 
                'GCACT': 6, 'CCACT': 6, 'GGACT': 6,'CGACT': 6,
                'GCACA': 6, 'CCACA': 6, 'GGACA': 6,'CGACA': 6}
SEQ_CUTOFF = 1


def search_motif(path, exp_design_name, size, mouse_seq, peak_technique):
    if peak_technique == 'All':
        PATH_PEAKS = path + '/PeakDetection/Peaks'
        peak_init_bed_filename = PATH_PEAKS + '/' + exp_design_name + '_MaxValues.bed'
        peak_bed_filename = PATH_PEAKS + '/' +  exp_design_name + '_Motif.bed'
    else:
        PATH_PEAKS = path + '/PeakDetection/' + peak_technique + '/'
        peak_init_bed_filename = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_MaxValues.bed'
        peak_bed_filename = PATH_PEAKS + exp_design_name + '_' + peak_technique + '_Motif.bed'

    # load bed file with all peaks
    with open(peak_bed_filename, "w") as bed_peaks_file, \
            open(peak_init_bed_filename, 'r') as bed_file:
        i=0
        for row in bed_file:
            chromo = row.split('\t')[0]
            begin = int(row.split('\t')[1])
            end = int(row.split('\t')[2])
            name = row.split('\t')[3]
            strand = row.split('\t')[6]
            # Calculate motif presence score
            sequence_score = 0
            sequence = mouse_seq[peak.chrom][peak.start:peak.end].seq
            for motif in MOTIF_METH:
                if motif in sequence:
                    sequence_score += MOTIF_METH[motif]
            strand = row.split('\t')[6]
            if sequence_score > SEQ_CUTOFF:
                peak_id = 'Peak_' + i
                i+=1
                bed_peaks_file.write(chromo + '\t' + begin + '\t' + end +
                               '\t' + name + '\tsequence_score\tstrand\n')


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
        opts, args = getopt.getopt(argv,"hp:e:t:g:", ["path=", "expdesign=","peak_technique=","genome="])
    except getopt.GetoptError:
        print('Cannot run command - Help: search_motif.py -p <path> -e <expdesign> -t <peak_technique> '
              '-g <genome>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('find_peak_center.py -p <path> -e <expdesign> -t <peak_technique> -g <genome>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-e", "--expdesign"):
            exp_design_name = arg
        elif opt in ("-t", "--peaktechnique"):
            peak_technique = arg
        elif opt in ("-g", "--genome"):
            genome_file = arg
    mouse_seq = read_mouse_seq(path, genome_file)
    search_motif(path, exp_design_name, mouse_seq, peak_technique)


if __name__ == "__main__":
    main(sys.argv[1:])
