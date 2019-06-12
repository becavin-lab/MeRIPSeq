#!/usr/bin/env python

class SequencePython(object):
    ##
    # Create a Sequence object<br>
    # * if(begin>end) <br>
    # * this.end = begin <br>
    # * this.begin = end <br>
    # * this.strand = '-' <br>
    # * @param name
    #  @param begin
    #  @param end
    ##
    def __init__(self, name, begin, end, strand, genome, chromosome):
        self.name = name
        self.begin = begin
        self.end = end
        self.genome = genome
        self.chromosome = chromosome

        if (strand == '+' or strand == '-'):
            self.strand = strand;
        else:
            print("Wrong format for strand")

        # Be sure than from < to
        if (begin <= end):
            self.begin = begin
            self.end = end
            self.strand = strand
        elif end != -1:
            if strand == '-':
                self.strand = '-'
                self.begin = self.end
                self.end = begin;
            else:
                print(name + " has wrong strand information")
                self.strand = '-'
                self.begin = end
                self.end = begin
        else:
            self.begin, self.end, self.strand = (begin, end, strand)

        self.beginIntersect = begin
        self.endIntersect = end
        self.length = self.end - self.begin + 1

    ##
    # * Decide if two genes are overlapping <br>
    # * First, calculate the overlapping size of two genes, finding the gene upstream to the other<br>
    # * If overlap is negative it will mean that we have no overlap (gene1.end < gene2.begin) <br>
    # * @param gene1
    # * @param gene2
    # * @return
    ##
    def isoverlap(self, seq2):
        overlapLength = self.overlap(seq2)
        if (overlapLength <= seq2.length and overlapLength > -10):
            return True
        else:
            return False

    ##
    # Calculate the overlapping size of two genes <br>
    # * First, find the gene upstream to the other<br>
    # * Then, as we got gene1.begin < gene2.begin, overlap = gene1.end - gene2.begin<br>
    # * overlap can be negative, it will mean that we have no overlap (gene1.end < gene2.begin) <br>
    # * or overlap can be higher than gene2.length, in that case overlap = gene2.length
    # *
    # * @param gene1
    # * @param gene2
    # * @return
    ## */
    def overlap(self, seq2):
        begin1 = self.begin
        end1 = self.end
        begin2 = seq2.begin
        end2 = seq2.end

        ## found gene upstream to the other
        overlapLength = -1000000;
        if (self.genome == seq2.genome and self.chromosome == seq2.chromosome):
            if (begin1 < begin2):
                overlapLength = end1 - begin2;
                if (overlapLength > 0 and seq2.length <= overlapLength):
                    overlapLength = seq2.length
            else:
                overlapLength = end2 - begin1
                if (overlapLength > 0 and self.length <= overlapLength):
                    overlapLength = self.length

        return overlapLength

    ##
    # @param seq2
    ##
    def getOverlapPosition(self, seq2):

        ## found gene upstream to the other
        if self.isoverlap(seq2):
            begin = min(self.begin, seq2.begin)
            end = max(self.end, seq2.end)
            beginIntersect = max(self.beginIntersect, seq2.beginIntersect)
            endIntersect = min(self.endIntersect, seq2.endIntersect)

        return [begin, end, beginIntersect, endIntersect]

    ##
    # * Decide if two sequences have same strand
    # * @param gene1
    # * @param gene2
    # * @return
    ##
    def isSameStrand(self, seq2):
        xnor = (self.strand == seq2.strand)
        return xnor

