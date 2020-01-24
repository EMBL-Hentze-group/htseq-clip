import unittest
import os
import gzip

from argparse import Namespace
from clip.createMatrix import MatrixConverter

class TestCreateMatrix(unittest.TestCase):

    outDir = 'tests/testcreateMatrix/UnitTestOutput'

    @classmethod
    def setUpClass(cls):
        if not os.path.exists(cls.outDir):
            os.makedirs(cls.outDir)
    
    def test_createMatrix(self):
        inputDir = 'tests/testcreateMatrix/'
        inputPostfix = 'txt.gz'
        outputFilename = os.path.join(self.outDir,'test_out_matrix.txt.gz')
        mC = MatrixConverter(inputDir=inputDir,inputPrefix="",inputPostfix=inputPostfix,outputFilename=outputFilename)
        mC.read_samples()
        mC.write_matrix()
        outMat,expMat = set(),set()
        with open('tests/testcreateMatrix/checktestcreateMatrix/out_matrix.txt','r') as ec:
            for l in ec:
                expMat.add(l)
        with gzip.open(outputFilename,'r') as oc:
            for l in oc:
                outMat.add(l.decode('utf-8'))
        self.assertSetEqual(outMat,expMat)

if __name__ == '__main__':
    unittest.main()