"""
Class for output file management
"""
import gzip
import os
import sys

class Output:
    # constructor 
    def __init__(self, fileName):
        # if file name is empty, redirect to stdin
        self._writer = self._returnStr
        if (fileName == "") or (fileName is None):
            self.fileOut = sys.stdout
        # if file name is given, file is opened 
        else:
            if fileName.lower().endswith((".gz","gzip")):
                self.fileOut = gzip.open(fileName, 'w')
                self._writer = self._returnByte
            else:
                self.fileOut = open(fileName, 'w')
    
    def _returnStr(self,s):
        '''
        helper function
        if a string is given, return it as such
        '''
        return s
    
    def _returnByte(self,s):
        '''
        helper function
        if a string is given, return byte
        '''
        return s.encode('utf-8')

    # output
    def write(self, s):
        s = self._writer(s)
        self.fileOut.write(s)

    # close
    def close(self):
        if not self.fileOut.closed:
            self.fileOut.close()
    
    # destructor
    def __del__(self):
        self.close()