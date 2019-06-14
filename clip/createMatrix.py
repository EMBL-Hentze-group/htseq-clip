__author__ = 'Tom'

import gzip
import os
import re
import sys

from output import Output

class MatrixConverter:

    """
    Class to convert count files into count matrixes.
    """

    def __init__(self, inputDir, inputPrefix, inputPostfix, outputFilename, annotation=None):
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
        
        # sample names are stored here for consistent output
        self.samplenamesList = None
        
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
        return [filename for filename in os.listdir(dirname) if self._file_filter(filename, prefix, postfix)]
    
    """
    read in all samples
    values will be stored in 
    """
    def read_samples(self):
    	# for each file name
        for file in self.inputFilenames:
            samplename = file.split(".")[0]
            self.allSamples.add(samplename)
            sys.stderr.write('Reading file {}\n'.format(file))
            # @TODO:  use logging module for all messages
            # open file and read in the content and store the results 
            with self._file_reader(os.path.join(self.inputDir, file)) as f:
                for linecount, line in enumerate(f):
                    if(linecount == 0):
                        continue
                    linesplit = line.strip().split("\t")
                    try:
                        self.countDict[ linesplit[0] ][ samplename ] = linesplit[3]
                    except KeyError:
                        self.countDict[ linesplit[0] ] = { samplename : linesplit[3] }
        self._init_samplenames_list()
    
    def _file_reader(self,fn):
        '''
        Helper function, return the correct file reade object based on file suffix
        Argument:
         fn: file name as string
        '''
        if re.match(r'.*\.gz.*',fn,re.IGNORECASE):
            return gzip.open(fn,'rt')
        else:
            return open(fn,'r')

    """
    Helper function
    getter for sample names
    """
    def _init_samplenames_list(self):
        self.samplenamesList = sorted(self.allSamples)
    
    """
    Helper function
    getter for header
    """
    def _get_header(self):
        return("\t".join(["unique_id"] + self.samplenamesList))
    
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
            for sample_name in self.samplenamesList:
                try:
                    outList.append(sample_count_dict[sample_name])  
                except KeyError:
                    outList.append('0')
            self.out.write('\t'.join(outList) + "\n")
        self.out.close()


'''
parser = argparse.ArgumentParser(description = 'Creating a count matrix from htseq-clip count files.')


parser.add_argument('--dir',
                    action   = 'store',
                    dest     = 'dir',
                    required = True, 
                    help     = 'directory of count files')
                    
parser.add_argument('--prefix',
                    action   = 'store',
                    dest     = 'prefix',
                    default  = '',
                    help     = 'prefix of count files')

parser.add_argument('--postfix',
                    action   = 'store',
                    dest     = 'postfix',
                    required = True,
                    help     = 'postfix of count files')
                    
parser.add_argument('--output',
                    action   = 'store',
                    dest     = 'output',
                    required = True, 
                    help     = 'destination of output count matrix')

args = parser.parse_args()


#print args.dir
#print args.prefix
#print args.postfix
#print args.output


m1 = MatrixConverter(args.dir, args.prefix, args.postfix, args.output)
##print("Reading in window counts")

m1.read_samples()
#print("Writing window counts") 


#print m1.output_filename
#print m1.count_dict
#print m1.allsamples
#print m1.input_dir
#print m1.samplenames_list

m1.write_matrix()
'''

