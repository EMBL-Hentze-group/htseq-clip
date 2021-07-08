import decimal
import gzip
import logging
import multiprocessing as mp
import re
import sys
import tempfile
from functools import partial
from os import sched_getaffinity
from pathlib import Path
from shutil import copyfileobj, rmtree

import numpy as np
import pysam

'''
bamCLIP module:
    Process bam files and extract crosslink sites.
    Input bam files MUST be co-ordinated sorted and indexed.

Authors: Sudeep Sahadevan, sudeep.sahadevan@embl.de
          Thomas Schwarzl, thomas.schwarzl@embl.de
Instituion: EMBL Heidelberg
'''
class bamCLIP:
        
    def __init__(self, options):
        self.fInput = options.input # input bam
        self.fOutput = options.output # output file
        self.choice = options.choice # choice
        self.mate = options.mate # mate info
        self.min_len = options.minReadLength # min read len
        self.max_len = options.maxReadLength # max read len
        self.primary = options.primary # flag use primary
        self.max_interval_length = options.maxReadIntervalLength # splice length,
        self.min_qual = options.minAlignmentQuality # alignment quality
        self.chromFile = options.chromFile # count sites only for chromosomes in this file
        self.cores = options.cores # number of cores to use
        if options.tmp is None: # use parent folder of output file as tmp folder
            self.tmp = Path(self.fOutput).parent/next(tempfile._get_candidate_names())
        else: # use given dir as tmp folder
            self.tmp = Path(options.tmp).parent/next(tempfile._get_candidate_names())
        logging.info('Using {} as tmp folder'.format(self.tmp))
    
    def __enter__(self):
        self.chromes = list()
        # find all chromosomes in the given bam file
        self._bam_checker()
        self._set_cores()
        # check if the output need to be gzipped
        self.gz = False
        if self.fOutput is not None:
            if re.match(r'^.*\.gz.*$',self.fOutput,re.IGNORECASE):
                self.gz = True
        # temp folder
        try:
            self.tmp.mkdir()
        except FileExistsError:
            logging.warning('Files in {} might be re-written!'.format(str(self.tmp)))
        return self

    def __exit__(self,except_type,except_val,except_traceback):
        # clean up
        rmtree(self.tmp)
        if except_type:
            logging.exception(except_val)

    @staticmethod
    def _read_chromosomes(fin):
        '''
        Helper function
        Given a TAB seperated file with chromosome names as first column, parse them
        '''
        if fin is None:
            return None
        else:
            chromosomes = set()
            with open(fin,'r') as fh:
                for ch in fh:
                    chromosomes.add(ch.strip().split('\t')[0])
            return chromosomes
    
    def _bam_checker(self):
        '''
        Helper function
        find names of all chromosmes in the bam file,
        if a file with list of chromosomes are given, filter chromosomes based on the list
        '''
        with pysam.AlignmentFile(self.fInput,mode='rb',check_sq=True,check_header=True,require_index=True) as _bh:
            bam_header = dict(_bh.header)
            if 'SQ' not in bam_header:
                msg = 'Cannot find @SQ header lines in file {}'.format(self.fInput)
                logging.error(msg)
                raise LookupError(msg)
            chroms = set(map(lambda s: s['SN'],bam_header['SQ']))
        in_chroms = bamCLIP._read_chromosomes(self.chromFile)
        if in_chroms is not None:
            chroms = chroms & in_chroms
            diff_chroms = in_chroms - chroms
            if len(diff_chroms)>0:
                logging.warning('Cannot find reads for {} chromosomes in {}'.format(', '.join(diff_chroms), self.fInput))
        if len(chroms) == 0:
            msg = 'Cannot find common chromosmes between bam file: {} and chromosome file {}'.format(self.fInput, self.chromFile)
            logging.error(msg)
            raise ValueError(msg)
        logging.info('Paring reads in {} chromosomes from {}'.format(len(chroms),self.fInput))
        # Final list of chromosomes to work with
        self.chromes = sorted(chroms)
    
    def _set_cores(self):
        '''
        Helper function
        set number of cores
        '''
        allcores = len(sched_getaffinity(0))
        if self.cores > allcores:
            setcores =  max(allcores-1,1)
            logging.warning('Give number of cores {} > number of cores detected {}. Setting cores to {}'.format(self.cores,allcores,setcores))
            self.cores = setcores
        elif allcores==1:
            logging.warning('Available # cores: 1, resetting cores parameter from {} to 1'.format(self.cores))
            self.cores = 1
        else:
            logging.info('Using {} cores out of {}...'.format(self.cores, allcores))
    
    def _create_tmp_files(self,op='start'):
        '''
        Helper function, generate temp. file names
        '''
        ext = '.bed.gz' if self.gz else '.bed'
        tmpDict = {}
        for chrom in self.chromes:
            tmpDict[chrom] = str(self.tmp/'{}_{}{}{}'.format(chrom,op,next(tempfile._get_candidate_names()),ext))
        return tmpDict
    
    def _write_output(self,tmp_dict):
        '''
        Helper function, copy contents of tmp files to self.output   

        Sorting function credit goes to:
        https://stackoverflow.com/questions/4813061/non-alphanumeric-list-order-from-os-listdir/48030307#48030307
        '''
        # sort chromosomes first
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        chroms = sorted(tmp_dict.keys(),key = alphanum_key)
        if (self.fOutput is None) or (self.fOutput==""): # no output file, dump to sys.stdout
            for chrom in chroms:
                copyfileobj(open(tmp_dict[chrom],'r'),sys.stdout)
        else:
            with open(self.fOutput,'wb') as dest:
                for chrom in chroms:
                    copyfileobj(open(tmp_dict[chrom],'rb'),dest)

    def extract_start_sites(self, offset=0):
        start_dict = self._create_tmp_files(op='start')
        pool = mp.Pool(processes = self.cores)
        for chrom, tmp_file in start_dict.items():
            pool.apply_async(_start_sites, 
                args=(self.fInput, chrom, tmp_file, self.gz , self.min_qual, self.min_len, self.max_len, self.max_interval_length, self.primary, self.mate, offset))
        pool.close()
        pool.join()
        self._write_output(start_dict)

    def extract_middle_sites(self):
        middle_dict = self._create_tmp_files(op='middle')
        pool = mp.Pool(processes = self.cores)
        for chrom, tmp_file in middle_dict.items():
            pool.apply_async(_middle_sites, 
                args=(self.fInput, chrom, tmp_file, self.gz , self.min_qual, self.min_len, self.max_len, self.max_interval_length, self.primary, self.mate))
        pool.close()
        pool.join()
        self._write_output(middle_dict)

    def extract_end_sites(self, offset=0):
        end_dict = self._create_tmp_files(op='end')
        pool = mp.Pool(processes = self.cores)
        for chrom, tmp_file in end_dict.items():
            pool.apply_async(_end_sites, 
                args=(self.fInput, chrom, tmp_file, self.gz , self.min_qual, self.min_len, self.max_len, self.max_interval_length, self.primary, self.mate, offset))
        pool.close()
        pool.join()
        self._write_output(end_dict)

    def extract_insertion_sites(self):
        insertion_dict = self._create_tmp_files(op='insertion')
        pool = mp.Pool(processes = self.cores)
        for chrom, tmp_file in insertion_dict.items():
            pool.apply_async(_insertion_sites, 
                args=(self.fInput, chrom, tmp_file, self.gz , self.min_qual, self.min_len, self.max_len, self.max_interval_length, self.primary, self.mate))
        pool.close()
        pool.join()
        self._write_output(insertion_dict)

    def extract_deletion_sites(self):
        deletion_dict = self._create_tmp_files(op='deletion')
        pool = mp.Pool(processes = self.cores)
        for chrom, tmp_file in deletion_dict.items():
            pool.apply_async(_deletion_sites, 
                args=(self.fInput, chrom, tmp_file, self.gz , self.min_qual, self.min_len, self.max_len, self.max_interval_length, self.primary, self.mate))
        pool.close()
        pool.join()
        self._write_output(deletion_dict)

def _discard_read(aln,qual = 10, min_len = 5, max_len = 100, max_interval_length = 10000, primary = False, mate = 2):
    '''
    Helper function, return True if any of the qc params fail
    Arugments:
        aln: pysam.AlignedSegment
        qual: read mapping quality
        min_len: min. read length
        max_len: max. read length
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site to process
    check qc params and return True if any of them fails
    '''
    if aln.is_duplicate or aln.is_qcfail or aln.is_unmapped  or (aln.mapping_quality < qual) or (aln.query_length < min_len) or (aln.query_length > max_len):
        # remove PCR duplicates, qc fails, unmapped ones or with low mapping quality or too short or too long
        return True
    if (primary and aln.is_secondary) or (aln.reference_length >= max_interval_length):
        # remove all secondary alignments if primary flag is given
        return True
    if mate == 1 and aln.is_read2:
        # if crosslink is on mate 1 and current read is mate 2
        return True
    elif mate == 2 and aln.is_read1:
        # if crosslink is on mate 2 and current read is mate 2
        return True
    return False

def _get_writer_encoder(is_gzip):
    '''
    Helper function
    Return file writer and string encoder for plain text files and gzipped files
    '''
    def _to_text(s):
        '''
        if a string is given, return it as such, to write to plain text
        '''
        return s
        
    def _to_byte(s):
        '''
        if a string is given, return byte, to write to a gz file
        '''
        return s.encode('utf-8')
    
    if is_gzip:
        fwriter = partial(gzip.open,mode='w')
        return (fwriter, _to_byte)
    else:
        fwriter = partial(open,mode='w')
        return (fwriter, _to_text)

def _start_sites(bam, chrom, outf, is_gzip = False, qual = 10, min_len = 5, max_len = 100, max_interval_length = 10000, primary = False, mate = 2, offset = 0):
    '''
    parse crosslink sites at the start positions
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        is_gzip: boolean, true if output is gzipped
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    '''
    fwriter, encoder = _get_writer_encoder(is_gzip)
    with pysam.AlignmentFile(bam,mode='rb') as bh, fwriter(outf) as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _discard_read(aln, qual = qual, min_len = min_len, max_len = max_len, max_interval_length = max_interval_length, primary = primary, mate = mate):
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
            oh.write(encoder('\t'.join(dat)+'\n'))

def _middle_sites(bam, chrom, outf, is_gzip = False, qual=10, min_len = 5, max_len = 100, max_interval_length = 10000, primary = False, mate = 2):
    '''
    parse crosslink sites at the middle positions
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        is_gzip: boolean, true if output is gzipped
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    '''
    fwriter, encoder = _get_writer_encoder(is_gzip)
    with pysam.AlignmentFile(bam,mode='rb') as bh, fwriter(outf) as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _discard_read(aln, qual = qual, min_len = min_len, max_len = max_len, max_interval_length = max_interval_length, primary = primary, mate = mate):
                continue
            alignedDict = dict(aln.aligned_pairs)
            mid = int(decimal.Decimal(aln.query_length/2).quantize(decimal.Decimal(1),rounding=decimal.ROUND_HALF_UP))
            if aln.is_reverse:
                if aln.reference_length%2==1:
                    # adjust for reverse strand and odd read length
                    mid-=1
                end = alignedDict[mid]
                if end is None:
                    logging.warning('Skipping {}, middle is an inserted position'.format(aln.query_name))
                    continue
                end = end-1
                begin = end -1
                strand = '-'
            else:
                end = alignedDict[mid]
                if end is None:
                    logging.warning('Skipping {}, middle is an inserted position'.format(aln.query_name))
                    continue
                begin = end -1
                strand = '+'
            if begin <0:
                continue
            try:
                yb  = aln.get_tag('YB')
            except KeyError:
                yb = 1
            dat = [chrom, str(begin), str(end), aln.query_name+'|'+ str(aln.query_length),str(yb),strand]
            oh.write(encoder('\t'.join(dat)+'\n'))

def _end_sites(bam, chrom, outf, is_gzip = False, qual=10, min_len = 5, max_len = 100, max_interval_length = 10000, primary = False, mate = 2, offset = 0):
    '''
    parse crosslink sites at the end positions
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        is_gzip: boolean, true if output is gzipped
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    '''
    fwriter, encoder = _get_writer_encoder(is_gzip)
    with pysam.AlignmentFile(bam,mode='rb') as bh, fwriter(outf) as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _discard_read(aln, qual = qual, min_len = min_len, max_len = max_len, max_interval_length = max_interval_length, primary = primary, mate = mate):
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
            oh.write(encoder('\t'.join(dat)+'\n'))

def _insertion_sites(bam, chrom, outf, is_gzip = False, qual=10, min_len = 5, max_len = 100, max_interval_length = 10000, primary = False, mate = 2):
    '''
    parse insertion sites
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        is_gzip: boolean, true if output is gzipped
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    '''
    fwriter, encoder = _get_writer_encoder(is_gzip)
    with pysam.AlignmentFile(bam,mode='rb') as bh, fwriter(outf) as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _discard_read(aln, qual = qual, min_len = min_len, max_len = max_len, max_interval_length = max_interval_length, primary = primary, mate = mate):
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
                oh.write(encoder('\t'.join(dat)+'\n'))

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

def _deletion_sites(bam, chrom, outf, is_gzip = False, qual=10, min_len = 5, max_len = 100, max_interval_length = 10000, primary = False, mate = 2):
    '''
    parse deletion sites
    Arugments:
        bam: bam file 
        chrom: current chromsome
        outf: output file name
        is_gzip: boolean, true if output is gzipped
        qual: read mapping quality
        max_interval_length: maximum alignment length 
        primary: flag to filter off secondary alignments
        mate: mate with crosslink site
        offset: number of basepairs to use as offset
    '''
    fwriter, encoder = _get_writer_encoder(is_gzip)
    with pysam.AlignmentFile(bam,mode='rb') as bh, fwriter(outf) as oh:
        for aln in bh.fetch(chrom,multiple_iterators=True):
            if _discard_read(aln, qual = qual, min_len = min_len, max_len = max_len, max_interval_length = max_interval_length, primary = primary, mate = mate):
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
                oh.write(encoder('\t'.join(dat)+'\n'))
