# --------------------------------------------------
# bamCLIP class
# --------------------------------------------------

import decimal
import gzip
import logging

import pysam
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


def _qc_fail(aln,qual = 10, max_interval_length = 10000, primary = False, mate = 2):
    '''
    Helper function
    Arugments:
        aln: pysam.AlignedSegment
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
    check qc params and return True if any of them fails
    '''
    if aln.is_duplicate or aln.is_qcfail or aln.is_unmapped  or (aln.mapping_quality < qual) or (aln.reference_length >= max_interval_length) or (primary and aln.is_secondary):
        # remove PCR duplicates, qc fails, unmapped ones or with low mapping quality
        return True
    if primary and aln.is_secondary:
        # remove all secondary alignments
        return True
    if mate == 1 and aln.is_read2:
        # if crosslink is on mate 1 and current read is mate 2
        return True
    elif mate == 2 and aln.is_read1:
        # if crosslink is on mate 2 and current read is mate 2
        return True
    return False
    
def get_start_sites(bam, chrom, outf, qual = 10, max_interval_length = 10000, primary = False, mate = 2, offset = 0 ):
    '''
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    parse crosslink sites at the start sites
    '''
    with pysam.AlignmentFile(bam,mode='rb') as bh, open(outf,'w') as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _qc_fail(aln, qual = qual, max_interval_length = max_interval_length, primary = primary, mate = mate):
                continue
            pos = list()
            strand = ''
            if aln.is_reverse:
                pos.append(aln.reference_end -offset -1)
                strand = '-'
            else:
                pos.append(aln.reference_start + offset)
                strand = '+'
            pos.append(pos[0]+1)
            pos = sorted(pos)
            try:
                yb  = aln.get_tag('YB')
            except KeyError:
                yb = 1
            dat = [chrom, str(pos[0]), str(pos[1]), aln.query_name+'|'+ str(aln.query_length),str(yb),strand]
            oh.write('\t'.join(dat)+'\n')

            


