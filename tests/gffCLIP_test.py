import os
import unittest
from argparse import Namespace

from clip.gffCLIP import gffCLIP
# these are custom defined exceptions
from clip.gffCLIP import EmptyFileException, NoFeaturesException, FeatureOrderException

class TestGFFCLIP(unittest.TestCase):
    '''
    Unit tests for module gffCLIP
    '''

    def test_empty(self):
        options = Namespace(gff='tests/testgFFCLIP/test_empty.gff3',output='tests/testgFFCLIP/checkgFFCLIP/test_empty_check.bed',
            id='gene_id',name='gene_name',type='gene_type',splitExons=True,unsorted=False)
        gffE = gffCLIP(options)
        with self.assertRaises(EmptyFileException):
            gffE.process()
    
    def test_comments(self):
        options = Namespace(gff='tests/testgFFCLIP/test_comments.gff3',output='tests/testgFFCLIP/checkgFFCLIP/test_comments_check.bed',
        id='gene_id',name='gene_name',type='gene_type',splitExons=True,unsorted=False)
        gffC = gffCLIP(options)
        with self.assertRaises(NoFeaturesException):
            gffC.process()
    
    def test_unsorted(self):
        options = Namespace(gff='tests/testgFFCLIP/test_unsorted.gff3',output='tests/testgFFCLIP/checkgFFCLIP/test_unsorted_check.bed',
        id='gene_id',name='gene_name',type='gene_type',splitExons=True,unsorted=False)
        gffU = gffCLIP(options)
        gffU.process(True) # unsorted gff file
        # cross check
        gffUnsorted,gffDefault = set(),set()
        with open('tests/testgFFCLIP/checkgFFCLIP/test_unsorted_check.bed','r') as un:
            for l in un:
                gffUnsorted.add(l)
        with open('tests/UnitTestOutput/testgFFCLIP/test_process_default.bed','r') as so:
            for l in so:
                gffDefault.add(l)
        self.assertSetEqual(gffUnsorted,gffDefault)
    
    def test_default(self):
        options = Namespace(gff='tests/testgFFCLIP/test.gff3',output='tests/testgFFCLIP/checkgFFCLIP/test_default_check.bed',
        id='gene_id',name='gene_name',type='gene_type',splitExons=True,unsorted=False)
        gffD = gffCLIP(options)
        gffD.process(False) # unsorted gff file
        # cross check
        gffSorted,gffDefault = set(),set()
        with open('tests/testgFFCLIP/checkgFFCLIP/test_default_check.bed','r') as un:
            for l in un:
                gffSorted.add(l)
        with open('tests/UnitTestOutput/testgFFCLIP/test_process_default.bed','r') as so:
            for l in so:
                gffDefault.add(l)
        self.assertSetEqual(gffSorted,gffDefault)
    
    def test_w50s20(self):
        options = Namespace(output='tests/testgFFCLIP/checkgFFCLIP/test_w50S20.bed',windowSize=50,windowStep=20)
        gffW50S20 = gffCLIP(options)
        gffW50S20.slidingWindow('tests/testgFFCLIP/checkgFFCLIP/test_default_check.bed')
        gffW,gffDefault = set(),set()
        with open('tests/testgFFCLIP/checkgFFCLIP/test_w50S20.bed','r') as un:
            for l in un:
                gffW.add(l)
        with open('tests/UnitTestOutput/testgFFCLIP/test_w50S20.bed','r') as so:
            for l in so:
                gffDefault.add(l)
        self.assertSetEqual(gffW,gffDefault)

if __name__ == '__main__':
    unittest.main()