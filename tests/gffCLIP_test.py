import gzip
import os
import unittest
from argparse import Namespace

# these are custom defined exceptions
from clip.gffCLIP import (EmptyFileException, FeatureOrderException,
                          NoFeaturesException, gffCLIP)


class TestGFFCLIP(unittest.TestCase):
    '''
    Unit tests for module gffCLIP
    '''
    outDir = 'tests/testgFFCLIP/UnitTestOutput/'

    @classmethod
    def setUpClass(cls):
        if not os.path.exists(cls.outDir):
            os.makedirs(cls.outDir)

    def test_empty(self):
        options = Namespace(gff='tests/testgFFCLIP/test_empty.gff3',output=os.path.join(self.outDir,'test_empty_check.bed'),
            id='gene_id',name='gene_name',type='gene_type',splitExons=True,unsorted=False)
        gffE = gffCLIP(options)
        with self.assertRaises(EmptyFileException):
            gffE.process()
    
    def test_comments(self):
        options = Namespace(gff='tests/testgFFCLIP/test_comments.gff3',output=os.path.join(self.outDir,'test_comments_check.bed'),
        id='gene_id',name='gene_name',type='gene_type',splitExons=True,unsorted=False)
        gffC = gffCLIP(options)
        with self.assertRaises(NoFeaturesException):
            gffC.process()
    
    def test_unsorted(self):
        '''
        unit test for subparser command annotation with flag --unsorted
        '''
        outFile = os.path.join(self.outDir,'test_unsorted_check.bed')
        options = Namespace(gff='tests/testgFFCLIP/test_unsorted.gff3',output=outFile,
        id='gene_id',name='gene_name',type='gene_type',splitExons=True,unsorted=False)
        gffU = gffCLIP(options)
        gffU.process(True) # unsorted gff file
        # cross check
        gffUnsorted,gffDefault = set(),set()
        with open('tests/testgFFCLIP/checkgFFCLIP/test_unsorted_check.bed','r') as so:
            for l in so:
                gffDefault.add(l)
        with open(outFile,'r') as un:
            for l in un:
                gffUnsorted.add(l)
        self.assertSetEqual(gffUnsorted,gffDefault)
    
    def test_default(self):
        '''
        unit test for subparser command annotation
        '''
        outFile = os.path.join(self.outDir,'test_process_default.bed')
        options = Namespace(gff='tests/testgFFCLIP/test.gff3',output=outFile,
        id='gene_id',name='gene_name',type='gene_type',splitExons=True,unsorted=False)
        gffD = gffCLIP(options)
        gffD.process(False) # unsorted gff file
        # cross check
        gffSorted,gffDefault = set(),set()
        with open('tests/testgFFCLIP/checkgFFCLIP/test_default_check.bed','r') as un:
            for l in un:
                gffSorted.add(l)
        with open(outFile,'r') as so:
            for l in so:
                gffDefault.add(l)
        self.assertSetEqual(gffSorted,gffDefault)
    
    def test_w50s20(self):
        '''
        unit test for subparser command createSlidingWindows
        with window size 50 & step size 20
        '''
        outFile = os.path.join(self.outDir,'test_w50S20.bed')
        options = Namespace(output=outFile,windowSize=50,windowStep=20)
        gffW50S20 = gffCLIP(options)
        gffW50S20.slidingWindow('tests/testgFFCLIP/checkgFFCLIP/test_default_check.bed')
        gffW,gffDefault = set(),set()
        with open('tests/testgFFCLIP/checkgFFCLIP/test_w50S20.bed','r') as so:
            for l in so:
                gffDefault.add(l)
        with open(outFile,'r') as sw:
            for l in sw:
                gffW.add(l)
        self.assertSetEqual(gffW,gffDefault)

    def test_w50s20_gz(self):
        '''
        unit test for subparser command createSlidingWindows
        with window size 50 & step size 20
        '''
        outFile = os.path.join(self.outDir,'test_w50S20.bed.gz')
        options = Namespace(output=outFile,windowSize=50,windowStep=20)
        gffW50S20 = gffCLIP(options)
        gffW50S20.slidingWindow('tests/testgFFCLIP/checkgFFCLIP/test_default_check.bed.gz')
        gffW,gffDefault = set(),set()
        with open('tests/testgFFCLIP/checkgFFCLIP/test_w50S20.bed','r') as so:
            for l in so:
                gffDefault.add(l)
        with gzip.open(outFile,'r') as sw:
            for l in sw:
                gffW.add(l.decode('utf-8'))
        self.assertSetEqual(gffW,gffDefault)
        
if __name__ == '__main__':
    unittest.main()
