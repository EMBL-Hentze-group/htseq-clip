import gzip
import logging
import os
import re
import sys
from pathlib import Path

from .output import Output

__author__ = 'Tom'
class MatrixConverter:

    """
    Class to convert count files into count matrixes.
    """

    def __init__(self, inputDir, inputPrefix, inputPostfix, outputFilename):
        '''
        Arguments:
         inputDir: input directory
         inputPrefix: prefix for input file(s)
         inputPostfix: postfix for input file(s)
         ouputFilename: output file name
         annotation: annotation file name
        '''
        # input files
        self.inputFilenames  = self._dir_filter(inputDir, prefix = inputPrefix, postfix = inputPostfix)
        # input dir
        self.inputDir = inputDir

        # output file
        self.out  = Output(outputFilename)
        
        # dict to store the counts
        self.countDict = {}
        
        # dict to store all samples
        self.allSamples = set()
        

        # decoder for bytes in gzip
        self._decoder = None
        
    """
    Helper function
    check if a filename matches a prefix or a postfix
    """
    def _file_filter(self, filename, prefix = "", postfix = ""):
        if not filename:
            return False
        filename = filename.strip()
        return(filename.startswith(prefix) and filename.endswith(postfix))
    
    """
    Helper function
    filter filenames in directory according to prefix and postfix
    """
    def _dir_filter(self, dirname, prefix = "", postfix = ""):   
        return [filename for filename in Path(dirname).glob("*") if self._file_filter(filename.name, prefix, postfix)]
    
    def read_samples(self, colNr=3):
        """
        read in all samples
        values will be stored in 
        colNr: zero based column index
        """
    	# for each file name
        self.countDict = {}
        for file in self.inputFilenames:
            samplename = file.stem.split('.')[0]
            self.allSamples.add(samplename)
            logging.info('Reading file {}'.format(file))
            # @TODO:  use logging module for all messages
            # open file and read in the content and store the results 
            with self._file_reader(file) as f:
                for linecount, line in enumerate(f):
                    if(linecount == 0):
                        continue
                    line = self._decoder(line)
                    linesplit = line.strip().split("\t")
                    try:
                        self.countDict[ linesplit[0] ][ samplename ] = linesplit[colNr]
                    except KeyError:
                        self.countDict[ linesplit[0] ] = { samplename : linesplit[colNr] }

    def _toStr(self,line):
        '''
        helper function
        given a string return it as it is
        '''
        return line
    
    def _byteToStr(self,line):
        '''
        helper function
        given bytes decode to string
        '''
        return line.decode('utf-8')

    def _file_reader(self,fn):
        '''
        Helper function, return the correct file reade object based on file suffix
        Argument:
         fn: file name as string
        '''
        if fn.name.lower().endswith((".gz",".gzip")):
            self._decoder = self._byteToStr
            return gzip.open(fn)
        else:
            self._decoder = self._toStr
            return open(fn)
    
    """
    Helper function
    getter for header
    """
    def _get_header(self):
        return("\t".join(["unique_id"] + sorted(self.allSamples)))
    
    """
    write matrix to output file
    """
    def write_matrix(self):
        # write header
        self.out.write(self._get_header() + "\n")
        # write rows
        for uid, sample_count_dict in self.countDict.items():
            outList = [uid]
            # write column
            for sample_name in sorted(self.allSamples):
                try:
                    outList.append(sample_count_dict[sample_name])  
                except KeyError:
                    outList.append('0')
            self.out.write('\t'.join(outList) + "\n")
        self.out.close()
