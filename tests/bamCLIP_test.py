import gzip
import os
import unittest
from argparse import Namespace

from clip.bamCLIP import bamCLIP
from clip.bamCLIP2 import bamCLIP as bc2


class TestBamCLIP(unittest.TestCase):
    '''
    Unit tests for module BamCLIP
    '''
    outDir = 'tests/testBamCLIP/UnitTestOutput'
    @classmethod
    def setUpClass(cls):
        if not os.path.exists(cls.outDir):
            os.makedirs(cls.outDir)

    '''
    test01.bam is supposed to raise a ValueError in all cases
    '''
    def test01_extract_SS(self):
        optionsSS = Namespace(input = "tests/testBamCLIP/test01.bam", output = os.path.join(self.outDir,"test01_SS.bed"), 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 500, maxReadIntervalLength=10000,
            primary = False, mate=1, choice='s', cores = 1, chromFile = None)
        bamS = bamCLIP(optionsSS)
        with self.assertRaises(ValueError):
            bamS.extract_StartSites()
    
    def test01_extract_SS2(self):
        optionsSS = Namespace(input = "tests/testBamCLIP/test01.bam", output = os.path.join(self.outDir,"test01_SS.bed"), 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 500, maxReadIntervalLength=10000,
            primary = False, mate=1, choice='s', cores = 1, chromFile = None, tmp = None)
        with self.assertRaises(ValueError):
            with bc2(optionsSS) as bh:
                bh.extract_start_sites()

    def test01_extract_MS(self):
        optionsSS = Namespace(input = "tests/testBamCLIP/test01.bam", output = os.path.join(self.outDir,"test01_MS.bed"), 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 0, maxReadIntervalLength=10000,
            primary = False,mate=1,choice='m' )
        bamS = bamCLIP(optionsSS)
        with self.assertRaises(ValueError):
            bamS.extract_MiddleSites()
    
    def test01_extract_MS2(self):
        optionsSS = Namespace(input = "tests/testBamCLIP/test01.bam", output = os.path.join(self.outDir,"test01_SS.bed"), 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 500, maxReadIntervalLength=10000,
            primary = False, mate=1, choice='s', cores = 1, chromFile = None, tmp = None)
        with self.assertRaises(ValueError):
            with bc2(optionsSS) as bh:
                bh.extract_middle_sites()
    
    def test01_extract_ES(self):
        optionsSS = Namespace(input = "tests/testBamCLIP/test01.bam", output = os.path.join(self.outDir,"test01_ES.bed"), 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 0, maxReadIntervalLength=10000,
            primary = False,mate=1,choice='e' )
        bamS = bamCLIP(optionsSS)
        with self.assertRaises(ValueError):
            bamS.extract_EndSites()
    
    def test01_extract_ES2(self):
        optionsSS = Namespace(input = "tests/testBamCLIP/test01.bam", output = os.path.join(self.outDir,"test01_SS.bed"), 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 500, maxReadIntervalLength=10000,
            primary = False, mate=1, choice='s', cores = 1, chromFile = None, tmp = None)
        with self.assertRaises(ValueError):
            with bc2(optionsSS) as bh:
                bh.extract_end_sites()
    
    '''
    test03.bam is supposed to finish successfully in all cases
    '''
    def test03_extract_SS(self):
        
        outFile = os.path.join(self.outDir,"test03_SS.bed")
        optionsSS = Namespace(input = "tests/testBamCLIP/test03.bam", output = outFile, 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 0, maxReadIntervalLength=10000,
            primary = False,mate=1,choice='s' )
         
        bamS = bamCLIP(optionsSS)
        bamS.extract_StartSites()
                
        v1 = []
        v2 = []
        with open(outFile, "r") as outputTest:
            for line in outputTest:
                v1.append(line)
        
        with open("tests/testBamCLIP/checktestBamCLIP/test03_SS_check.bed", "r") as checkTest:
            for line in checkTest:
                v2.append(line)
             
        self.assertEqual(v1, v2)
    
    def test03_extract_SS2(self):

        outFile = os.path.join(self.outDir,"test03_SS2.bed")
        optionsSS = Namespace(input = "tests/testBamCLIP/test03.bam", output = outFile, 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 500, maxReadIntervalLength=10000,
            primary = False, mate=1, choice='s', cores = 1, chromFile = None, tmp = None)
        with bc2(optionsSS) as bh:
            bh.extract_start_sites()
        test_out, orig_out = set(),set()
        with open(outFile, "r") as testH:
            for line in testH:
                test_out.add(line.strip())
        with open("tests/testBamCLIP/checktestBamCLIP/test03_SS_check.bed", "r") as origH:
            for line in origH:
                orig_out.add(line.strip())
        self.assertEqual(test_out, orig_out)
    
    def test03_extract_MS(self):
        
        outFile = os.path.join(self.outDir,"test03_MS.bed")
        optionsMS = Namespace(input = "tests/testBamCLIP/test03.bam", output = outFile, 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 0, maxReadIntervalLength=10000,
            primary = False,mate=1,choice='m' )
         
        bamM = bamCLIP(optionsMS)
        bamM.extract_MiddleSites()
                
        v1 = []
        v2 = []
        with open(outFile, "r") as outputTest:
            for line in outputTest:
                v1.append(line)
        
        with open("tests/testBamCLIP/checktestBamCLIP/test03_MS_check.bed", "r") as checkTest:
            for line in checkTest:
                v2.append(line)
             
        self.assertEqual(v1, v2)
    
    def test03_extract_MS2(self):

        outFile = os.path.join(self.outDir,"test03_MS2.bed")
        optionsSS = Namespace(input = "tests/testBamCLIP/test03.bam", output = outFile, 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 500, maxReadIntervalLength=10000,
            primary = False, mate=1, choice='s', cores = 1, chromFile = None, tmp = None)
        with bc2(optionsSS) as bh:
            bh.extract_middle_sites()
        test_out, orig_out = set(),set()
        with open(outFile, "r") as testH:
            for line in testH:
                test_out.add(line.strip())
        with open("tests/testBamCLIP/checktestBamCLIP/test03_MS_check.bed", "r") as origH:
            for line in origH:
                orig_out.add(line.strip())
        self.assertEqual(test_out, orig_out)
    
    def test03_extract_ES(self):
        
        outFile = os.path.join(self.outDir,"test03_ES.bed")
        optionsES = Namespace(input = "tests/testBamCLIP/test03.bam", output = outFile, 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 0, maxReadIntervalLength=10000,
            primary = False,mate=1,choice='e' )
         
        bamE = bamCLIP(optionsES)
        bamE.extract_EndSites()
                
        v1 = []
        v2 = []
        with open(outFile, "r") as outputTest:
            for line in outputTest:
                v1.append(line)
        
        with open("tests/testBamCLIP/checktestBamCLIP/test03_ES_check.bed", "r") as checkTest:
            for line in checkTest:
                v2.append(line)
             
        self.assertEqual(v1, v2)

    def test03_extract_ES2(self):

        outFile = os.path.join(self.outDir,"test03_ES2.bed")
        optionsSS = Namespace(input = "tests/testBamCLIP/test03.bam", output = outFile, 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 500, maxReadIntervalLength=10000,
            primary = False, mate=1, choice='s', cores = 1, chromFile = None, tmp = None)
        with bc2(optionsSS) as bh:
            bh.extract_end_sites()
        test_out, orig_out = set(),set()
        with open(outFile, "r") as testH:
            for line in testH:
                test_out.add(line.strip())
        with open("tests/testBamCLIP/checktestBamCLIP/test03_ES_check.bed", "r") as origH:
            for line in origH:
                orig_out.add(line.strip())
        self.assertEqual(test_out, orig_out)
    
    def test03_extract_DEL(self):
        
        outFile = os.path.join(self.outDir,"test03_DEL.bed")
        optionsD = Namespace(input = "tests/testBamCLIP/test03.bam", output = outFile, 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 0, maxReadIntervalLength=10000,
            primary = False,mate=1,choice='e' )
        bamD = bamCLIP(optionsD)
        bamD.extract_DeletionSites()
                
        v1 = []
        v2 = []
        with open(outFile, "r") as outputTest:
            for line in outputTest:
                v1.append(line)
        
        with open("tests/testBamCLIP/checktestBamCLIP/test03_DEL_check.bed", "r") as checkTest:
            for line in checkTest:
                v2.append(line)
             
        self.assertEqual(v1, v2)
    
    def test03_extract_INS(self):

        outFile = os.path.join(self.outDir,"test03_INS.bed")
        optionsI = Namespace(input = "tests/testBamCLIP/test03.bam", output = outFile, 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 0, maxReadIntervalLength=10000,
            primary = False,mate=1,choice='e' )
        bamI= bamCLIP(optionsI)
        bamI.extract_InsertionSites()
                
        v1 = []
        v2 = []
        with open(outFile, "r") as outputTest:
            for line in outputTest:
                v1.append(line)
        
        with open("tests/testBamCLIP/checktestBamCLIP/test03_INS_check.bed", "r") as checkTest:
            for line in checkTest:
                v2.append(line)
             
        self.assertEqual(v1, v2)
    
    def test03_extract_gz(self):

        outFile = os.path.join(self.outDir,"test03_INS.bed.gz")
        optionsI = Namespace(input = "tests/testBamCLIP/test03.bam", output = outFile, 
            minAlignmentQuality = 10, minReadLength = 0, maxReadLength = 0, maxReadIntervalLength=10000,
            primary = False,mate=1,choice='e' )
        bamI= bamCLIP(optionsI)
        bamI.extract_InsertionSites()
        v1,v2 = set(),set()
        with gzip.open(outFile,"r") as outgz:
            for line in outgz:
                v1.add(line.decode('utf-8')) 
        with open("tests/testBamCLIP/checktestBamCLIP/test03_INS_check.bed", "r") as checkTest:
            for line in checkTest:
                v2.add(line)
        self.assertSetEqual(v1,v2)


if __name__ == '__main__':
    unittest.main()
