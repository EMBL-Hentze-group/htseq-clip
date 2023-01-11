import gzip
import logging
import sys
from pathlib import Path

from .output import Output

class Trimmer:

    '''
    Trim large annotation files based on the output from createMatrix function
    '''

    def __init__(self, inputMatrix, inputAnn, outputAnn):
        '''
        Arguments:
         inputMatrix: output from createMatrix function
         inputAnn: output from mapToId function
         outputAnn: output  annotation function
        '''
        self.inputMatrix = inputMatrix
        self.inputAnn = inputAnn
        self.outputAnn = Output(outputAnn)
        self._uniq_ids = set()
    
    @classmethod
    def _get_reader_mode(self, fn):
        '''
        Helper function 
        Return correct parser and mode for plain text and .gz files
        '''
        with open(fn,'rb') as _rb:
            if _rb.read(2) == b'\x1f\x8b':
                logging.debug('{} is gzipped'.format(fn))
                return gzip.open, 'rt'
            else:
                logging.debug('{} is plain text'.format(fn))
                return open, 'r'
    
    def trim_annotation(self, header = True):
        '''
        trim annotation data based on read count matrix
        '''
        reader, mode = self._get_reader_mode(self.inputAnn)
        # parse matrix
        self._read_matrix()
        # write trimmed annotation
        with reader(self.inputAnn, mode) as annh:
            for l in annh:
                if header:
                    header = False
                    self.outputAnn.write(l)
                    continue
                ldat = l.split('\t')
                if ldat[0] in self._uniq_ids:
                    self.outputAnn.write(l)

    def _read_matrix(self):
        '''
        Helper function
        Read count data matrix and parse unique id columns (first column)
        '''
        reader, mode = self._get_reader_mode(self.inputMatrix)
        with reader(self.inputMatrix, mode) as _mh:
            for l in _mh:
                ldat = l.strip().split('\t')
                self._uniq_ids.add(ldat[0])
        if len(self._uniq_ids)==0:
            raise RuntimeError(f"Cannot parse unique ids from {self.inputMatrix}!")
        logging.info(f"Found {len(self._uniq_ids)} unique ids in {self.inputMatrix}")
    

