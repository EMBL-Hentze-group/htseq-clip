"""
Class for output file management
"""

import gzip, os, sys

class Output:
    # constructor 
    def __init__(self, fileName):
        # if file name is empty, redirect to stdin
        if (fileName == "") or (fileName is None):
            self.fileOut = sys.stdout
        # if file name is given, file is opened 
        else:
            if fileName.endswith(".gz"):
                self.fileOut = gzip.open(fileName, 'w')
            else:
                self.fileOut = open(fileName, 'w')

    # output
    def write(self, s):
        self.fileOut.write(s)

    # close
    def close(self):
        if not self.fileOut.closed:
            self.fileOut.close()
    
    # destructor
    def __del__(self):
        self.close()