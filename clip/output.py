"""
Class for output file management
"""

import gzip, os

class Output:
    # constructor 
    def __init__(self, fileName):
        # if file name is empty, redirect to stdin
        if fileName == "":
            self.writeFile = False
        # if file name is given, file is opened 
        else:
            self.writeFile = True

            if fileName.endswith(".gz"):
                self.fileOut = gzip.open(fileName, 'w')
            else:
                self.fileOut = open(fileName, 'w')

    # output
    def write(self, s):
        if self.writeFile:
            self.fileOut.write(s)
        else:
            print(s)

    # close
    def close(self):
        if self.writeFile and not self.fileOut.closed:
            self.fileOut.close()
    # destructor
    def __del__(self):
        self.close()
    
