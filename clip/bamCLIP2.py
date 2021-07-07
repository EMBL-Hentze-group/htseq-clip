# --------------------------------------------------
# bamCLIP class
# --------------------------------------------------

import decimal
import gzip
import logging

import pysam
import numpy as np
from .output import Output

class bamCLIP(object):
        
    def __init__(self, options):
        self.fInput = options.input
        self.writeFile = True
        self.fOutput = options.output
        self.choice = options.choice
        self.mate = options.mate           
        self.data = {'maxReadLength' : options.maxReadLength,
                     'minReadLength' : options.minReadLength,
                     'primary': options.primary,
                     'maxReadIntervalLength': options.maxReadIntervalLength,
                     'minAlignmentQuality': options.minAlignmentQuality}
    
    def extract_start_sites(self):
        # @TODO: Fill me up
        pass

    def extract_middle_sites(self):
        # @TODO: fill me up
        pass

    def extract_end_sites(self):
        # @TODO Fill me up
        pass

    def extract_insertion_sites(self):
        # @TODO Fill me up
        pass

    def extract_deletion_sites(self):
        # @TODO Fill me up
        pass

def _discard_read(aln,qual = 10, max_interval_length = 10000, primary = False, mate = 2):
    '''
    Helper function, return True if any of the qc params fail
    Arugments:
        aln: pysam.AlignedSegment
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
    check qc params and return True if any of them fails
    '''
    if aln.is_duplicate or aln.is_qcfail or aln.is_unmapped  or (aln.mapping_quality < qual) or (aln.reference_length >= max_interval_length):
        # remove PCR duplicates, qc fails, unmapped ones or with low mapping quality
        return True
    if primary and aln.is_secondary:
        # remove all secondary alignments if primary flag is given
        return True
    if mate == 1 and aln.is_read2:
        # if crosslink is on mate 1 and current read is mate 2
        return True
    elif mate == 2 and aln.is_read1:
        # if crosslink is on mate 2 and current read is mate 2
        return True
    return False
    
def _start_sites(bam, chrom, outf, qual = 10, max_interval_length = 10000, primary = False, mate = 2, offset = 0 ):
    '''
    parse crosslink sites at the start positions
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    '''
    with pysam.AlignmentFile(bam,mode='rb') as bh, open(outf,'w') as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _discard_read(aln, qual = qual, max_interval_length = max_interval_length, primary = primary, mate = mate):
                continue
            pos = list()
            strand = ''
            if aln.is_reverse:
                pos.append(aln.reference_end -offset -1)
                strand = '-'
            else:
                pos.append(aln.reference_start + offset)
                strand = '+'
            if pos[0] <0:
                continue
            pos.append(pos[0]+1)
            pos = sorted(pos)
            try:
                yb  = aln.get_tag('YB')
            except KeyError:
                yb = 1
            dat = [chrom, str(pos[0]), str(pos[1]), aln.query_name+'|'+ str(aln.query_length),str(yb),strand]
            oh.write('\t'.join(dat)+'\n')

def _middle_sites(bam, chrom, outf, qual = 10, max_interval_length = 10000, primary = False, mate = 2, offset = 0 ):
    '''
    parse crosslink sites at the middle positions
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    '''
    with pysam.AlignmentFile(bam,mode='rb') as bh, open(outf,'w') as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _discard_read(aln, qual = qual, max_interval_length = max_interval_length, primary = primary, mate = mate):
                continue
            alignedDict = dict(aln.aligned_pairs)
            mid = int(decimal.Decimal(aln.query_length/2).quantize(decimal.Decimal(1),rounding=decimal.ROUND_HALF_UP))
            end = alignedDict[mid]+offset
            begin = end-1
            if aln.is_reverse:
                if aln.reference_length%2==1:
                    # adjust for reverse strand and odd read length
                    mid-=1
                end = alignedDict[mid]-offset-1
                begin = end -1
                strand = '-'
            else:
                end = alignedDict[mid] + offset
                begin = end -1
                strand = '+'
            if begin <0:
                continue
            try:
                yb  = aln.get_tag('YB')
            except KeyError:
                yb = 1
            dat = [chrom, str(begin), str(end), aln.query_name+'|'+ str(aln.query_length),str(yb),strand]
            oh.write('\t'.join(dat)+'\n')

def _end_sites(bam, chrom, outf, qual = 10, max_interval_length = 10000, primary = False, mate = 2, offset = 0 ):
    '''
    parse crosslink sites at the end positions
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    '''
    with pysam.AlignmentFile(bam,mode='rb') as bh, open(outf,'w') as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _discard_read(aln, qual = qual, max_interval_length = max_interval_length, primary = primary, mate = mate):
                continue
            pos = list()
            strand = ''
            if aln.is_reverse:
                pos.append(aln.reference_start -offset -1)
                strand = '-'
            else:
                pos.append(aln.reference_end + offset)
                strand = '+'
            if pos[0] < 0:
                continue
            pos.append(pos[0]+1)
            pos = sorted(pos)
            try:
                yb  = aln.get_tag('YB')
            except KeyError:
                yb = 1
            dat = [chrom, str(pos[0]), str(pos[1]), aln.query_name+'|'+ str(aln.query_length),str(yb),strand]
            oh.write('\t'.join(dat)+'\n')

def _insertion_sites(bam, chrom, outf, qual = 10, max_interval_length = 10000, primary = False, mate = 2, offset = 0 ):
    '''
    parse insertion sites
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    '''
    with pysam.AlignmentFile(bam,mode='rb') as bh, open(outf,'w') as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _discard_read(aln, qual = qual, max_interval_length = max_interval_length, primary = primary, mate = mate):
                continue
            if 1 not in set(map(lambda x: x[0], aln.cigartuples)):
                # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats
                # 1 indicates an insertion
                continue
            strand = '-' if aln.is_reverse else '+'
            try:
                yb  = aln.get_tag('YB')
            except KeyError:
                yb = 1
            positions = np.array(list(map(lambda x: x[1],aln.get_aligned_pairs())))
            # inserted bases have 'None' as position info
            ins_pos_fn = _insertion_positions(aln.is_reverse)
            inspositions= np.where(positions==None)[0]
            for ins in np.split(inspositions, np.where(np.diff(inspositions) != 1)[0]+1):
                begin, end = ins_pos_fn(ins)
                dat = [chrom, str(positions[begin]), str(positions[end]), aln.query_name+'|'+ str(aln.query_length),str(yb),strand]
                oh.write('\t'.join(dat)+'\n')

def _insertion_positions(is_reverse):
    '''
    Helper function, decide insertion positions based on strand info
    ''' 
    def pos_strand(insertions):
        return [np.min(insertions)-1,np.max(insertions)+1]
    
    def neg_strand(insertions):
        end = np.min(insertions)
        return [end-2, end-1]
    
    if is_reverse:
        return neg_strand
    else:
        return pos_strand

def _deletion_sites(bam, chrom, outf, qual = 10, max_interval_length = 10000, primary = False, mate = 2, offset = 0 ):
    '''
    parse deletion sites
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    '''
    with pysam.AlignmentFile(bam,mode='rb') as bh, open(outf,'w') as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _discard_read(aln, qual = qual, max_interval_length = max_interval_length, primary = primary, mate = mate):
                continue
            if 2 not in set(map(lambda x: x[0], aln.cigartuples)):
                # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats
                # 2 indicates deletion
                continue
            strand = '-' if aln.is_reverse else '+'
            try:
                yb  = aln.get_tag('YB')
            except KeyError:
                yb = 1
            positions = np.array(aln.get_aligned_pairs())
            delpositions= np.where(positions[:,0]==None)[0]
            for dels in np.split(delpositions, np.where(np.diff(delpositions) != 1)[0]+1):
                begin = np.min(dels)-1
                end = np.max(dels)
                dat = [chrom, str(positions[begin,1]), str(positions[end,1]), aln.query_name+'|'+ str(aln.query_length),str(yb),strand]
                oh.write('\t'.join(dat)+'\n')







                


