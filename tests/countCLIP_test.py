import gzip
import os
import unittest
from argparse import Namespace

from clip.countCLIP import countCLIP


class TestcountCLIP(unittest.TestCase):
    '''
    unit tests for module countCLIP
    '''
    outDir = 'tests/testcountCLIP/UnitTestOutput'

    @classmethod
    def setUpClass(cls):
        if not os.path.exists(cls.outDir):
            os.makedirs(cls.outDir)

    def test_annotation(self):
        '''
        unit test for subparser command mapToId
        '''
        outFile = os.path.join(self.outDir,'test_annotation.txt')
        options = Namespace(annotation='tests/testcountCLIP/W50S20.bed',output=outFile)
        countC = countCLIP(options)
        countC.annotationToIDs()
        annSet, expSet = set(),set()
        with open('tests/testcountCLIP/checkcountCLIP/test_annotation.txt','r') as eann:
            for l in eann:
                expSet.add(l)
        with open(outFile,'r') as ann:
            for l in ann:
                annSet.add(l)
        self.assertSetEqual(annSet,expSet)
    
    def test_annotation_gz(self):
        '''
        unit test for subparser command mapToId
        '''
        outFile = os.path.join(self.outDir,'test_annotation.txt.gz')
        options = Namespace(annotation='tests/testcountCLIP/W50S20.bed.gz',output=outFile)
        countC = countCLIP(options)
        countC.annotationToIDs()
        annSet, expSet = set(),set()
        with open('tests/testcountCLIP/checkcountCLIP/test_annotation.txt','r') as eann:
            for l in eann:
                expSet.add(l)
        with gzip.open(outFile,'r') as ann:
            for l in ann:
                annSet.add(l.decode('utf-8'))
        self.assertSetEqual(annSet,expSet)
    
    def test_count_slidingWindows(self):
        '''
        unit test for subparser command count
        '''
        outFile = os.path.join(self.outDir,'test_slidingWindows_count.txt.gz')
        options = Namespace(annotation='tests/testcountCLIP/W50S20.bed',output=outFile,
            input='tests/testcountCLIP/end_sites.bed.gz')
        countSW = countCLIP(options)
        countSW.count()
        countSet, defaultSet = set(),set()
        with open('tests/testcountCLIP/checkcountCLIP/test_slidingWindows_count.txt','r') as ds:
            for l in ds:
                defaultSet.add(l)
        with gzip.open(outFile,'r') as cs:
            for l in cs:
                countSet.add(l.decode('utf-8'))
        self.assertSetEqual(countSet,defaultSet)
        
if __name__ == '__main__':
    unittest.main()
