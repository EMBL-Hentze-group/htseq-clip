# --------------------------------------------------
# UNIT TEST CLASS
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: January 2016
# --------------------------------------------------

import unittest
import sys
sys.path.append("..")

from clip import bamCLIP
from clip import bedCLIP
from clip import gtfCLIP


####################################################################################
############################EXTRACTION UNIT TESTS###################################
####################################################################################
class TestExtractMethods(unittest.TestCase):
     
    def test01_extract_SS(self):
        
        options = {}
        
        bamC = bamCLIP.bamCLIP(options)  
        
        bamC.fInput = "testExtract/test01.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test01_SS.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False  
        
        with self.assertRaises(ValueError):
            bamC.extract_StartSites()
       
    def test01_extract_MS(self):
          
        options = {}
        
        bamC = bamCLIP.bamCLIP(options)  
        
        bamC.fInput = "testExtract/test01.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test01_MS.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False  
        
        with self.assertRaises(ValueError):
            bamC.extract_MiddleSites()
              
    def test01_extract_ES(self):
          
        options = {}
        
        bamC = bamCLIP.bamCLIP(options)  
        
        bamC.fInput = "testExtract/test01.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test01_ES.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False  
        
        with self.assertRaises(ValueError):
            bamC.extract_EndSites()
              
    def test01_extract_DEL(self):
              
        options = {}
        
        bamC = bamCLIP.bamCLIP(options)  
        
        bamC.fInput = "testExtract/test01.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test01_DEL.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False  
        
        with self.assertRaises(ValueError):
            bamC.extract_DeletionSites()
              
    def test01_extract_INS(self):
          
        options = {}
        
        bamC = bamCLIP.bamCLIP(options)  
        
        bamC.fInput = "testExtract/test01.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test01_INS.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False  
        
        with self.assertRaises(ValueError):
            bamC.extract_InsertionSites()
     
    def test02_extract_SS(self):
        
        options = {}
         
        bamC = bamCLIP.bamCLIP(options)
        
        bamC.fInput = "testExtract/test02.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test02_SS.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        bamC.choice = "s"
        
        bamC.extract_StartSites()
            
        outputTest = open("UnitTestOutput/testExtract/test02_SS.bed", "r")
        checkTest  = open("testExtract/checkTestExtract/test02_SS_check.bed", "r")
        
        v1 = []
        v2 = []
        
        for line in outputTest:
            v1.append(line)
        
        for line in checkTest:
            v2.append(line)
            
        outputTest.close()
        checkTest.close()
             
        self.assertEqual(v1, v2)
         
    def test02_extract_MS(self):
         
        options = {}
         
        bamC = bamCLIP.bamCLIP(options)
        
        bamC.fInput = "testExtract/test02.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test02_MS.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        
        bamC.extract_MiddleSites()
            
        outputTest = open("UnitTestOutput/testExtract/test02_MS.bed", "r")
        checkTest  = open("testExtract/checkTestExtract/test02_MS_check.bed", "r")
        
        v1 = []
        v2 = []
        
        for line in outputTest:
            v1.append(line)
        
        for line in checkTest:
            v2.append(line)
            
        outputTest.close()
        checkTest.close()
             
        self.assertEqual(v1, v2)
        
         
    def test02_extract_ES(self):
        
        options = {}
         
        bamC = bamCLIP.bamCLIP(options)
        
        bamC.fInput = "testExtract/test02.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test02_ES.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        
        bamC.extract_EndSites()
            
        outputTest = open("UnitTestOutput/testExtract/test02_ES.bed", "r")
        checkTest  = open("testExtract/checkTestExtract/test02_ES_check.bed", "r")
        
        v1 = []
        v2 = []
        
        for line in outputTest:
            v1.append(line)
        
        for line in checkTest:
            v2.append(line)
            
        outputTest.close()
        checkTest.close()
             
        self.assertEqual(v1, v2)
         
    def test02_extract_DEL(self):
        
        options = {}
         
        bamC = bamCLIP.bamCLIP(options)
        
        bamC.fInput = "testExtract/test02.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test02_DEL.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        
        bamC.extract_DeletionSites()
            
        outputTest = open("UnitTestOutput/testExtract/test02_DEL.bed", "r")
        checkTest  = open("testExtract/checkTestExtract/test02_DEL_check.bed", "r")
        
        v1 = []
        v2 = []
        
        for line in outputTest:
            v1.append(line)
        
        for line in checkTest:
            v2.append(line)
            
        outputTest.close()
        checkTest.close()
             
        self.assertEqual(v1, v2)
              
    def test02_extract_INS(self):
         
        options = {}
         
        bamC = bamCLIP.bamCLIP(options)
        
        bamC.fInput = "testExtract/test02.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test02_INS.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        
        bamC.extract_InsertionSites()
            
        outputTest = open("UnitTestOutput/testExtract/test02_INS.bed", "r")
        checkTest  = open("testExtract/checkTestExtract/test02_INS_check.bed", "r")
        
        v1 = []
        v2 = []
        
        for line in outputTest:
            v1.append(line)
        
        for line in checkTest:
            v2.append(line)
            
        outputTest.close()
        checkTest.close()
             
        self.assertEqual(v1, v2)
        
    def test03_extract_SS(self):
        
        options = {}
         
        bamC = bamCLIP.bamCLIP(options)
        
        bamC.fInput = "testExtract/test03.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test03_SS.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        bamC.choice = "s"
        
        bamC.extract_StartSites()
            
        outputTest = open("UnitTestOutput/testExtract/test03_SS.bed", "r")
        checkTest  = open("testExtract/checkTestExtract/test03_SS_check.bed", "r")
        
        v1 = []
        v2 = []
        
        for line in outputTest:
            v1.append(line)
        
        for line in checkTest:
            v2.append(line)
            
        outputTest.close()
        checkTest.close()
             
        self.assertEqual(v1, v2)
         
    def test03_extract_MS(self):
         
        options = {}
         
        bamC = bamCLIP.bamCLIP(options)
        
        bamC.fInput = "testExtract/test03.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test03_MS.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        
        bamC.extract_MiddleSites()
            
        outputTest = open("UnitTestOutput/testExtract/test03_MS.bed", "r")
        checkTest  = open("testExtract/checkTestExtract/test03_MS_check.bed", "r")
        
        v1 = []
        v2 = []
        
        for line in outputTest:
            v1.append(line)
        
        for line in checkTest:
            v2.append(line)
            
        outputTest.close()
        checkTest.close()
             
        self.assertEqual(v1, v2)
              
    def test03_extract_ES(self):
        
        options = {}
         
        bamC = bamCLIP.bamCLIP(options)
        
        bamC.fInput = "testExtract/test03.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test03_ES.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        
        bamC.extract_EndSites()
            
        outputTest = open("UnitTestOutput/testExtract/test03_ES.bed", "r")
        checkTest  = open("testExtract/checkTestExtract/test03_ES_check.bed", "r")
        
        v1 = []
        v2 = []
        
        for line in outputTest:
            v1.append(line)
        
        for line in checkTest:
            v2.append(line)
            
        outputTest.close()
        checkTest.close()
             
        self.assertEqual(v1, v2)
               
    def test03_extract_DEL(self):
        
        options = {}
         
        bamC = bamCLIP.bamCLIP(options)
        
        bamC.fInput = "testExtract/test03.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test03_DEL.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        
        bamC.extract_DeletionSites()
            
        outputTest = open("UnitTestOutput/testExtract/test03_DEL.bed", "r")
        checkTest  = open("testExtract/checkTestExtract/test03_DEL_check.bed", "r")
        
        v1 = []
        v2 = []
        
        for line in outputTest:
            v1.append(line)
        
        for line in checkTest:
            v2.append(line)
            
        outputTest.close()
        checkTest.close()
             
        self.assertEqual(v1, v2)
             
    def test03_extract_INS(self):
         
        options = {}
         
        bamC = bamCLIP.bamCLIP(options)
        
        bamC.fInput = "testExtract/test03.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test03_INS.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        
        bamC.extract_InsertionSites()
            
        outputTest = open("UnitTestOutput/testExtract/test03_INS.bed", "r")
        checkTest  = open("testExtract/checkTestExtract/test03_INS_check.bed", "r")
        
        v1 = []
        v2 = []
        
        for line in outputTest:
            v1.append(line)
        
        for line in checkTest:
            v2.append(line)
            
        outputTest.close()
        checkTest.close()
             
        self.assertEqual(v1, v2)

####################################################################################
###############################PROCESS UNIT TESTS###################################
####################################################################################
class TestProcessMethods(unittest.TestCase):
     
    def test_process(self):
        
        options = {}
        
        gtfCLIP.gtfCLIP.gtfFile = "testProcess/test.gtf"
        gtfCLIP.gtfCLIP.fOutput = "UnitTestOutput/testProcess/test.bed"
        gtfCLIP.gtfCLIP.geneType = "gene_type"
         
        gtfC = gtfCLIP.gtfCLIP(options)
        gtfC.processGTF()
        outputTest = open("UnitTestOutput/testProcess/test.bed", "r")
        checkTest  = open("testProcess/checkTestProcess/test_check.bed", "r")
         
        v1 = []
        v2 = []
         
        for line in outputTest:
            v1.append(line)
         
        for line in checkTest:
            v2.append(line)
             
        outputTest.close()
        checkTest.close()
              
        self.assertEqual(v1, v2)
         
    def test_process_empty(self):
        
        options = {}
        
        gtfCLIP.gtfCLIP.gtfFile = "testProcess/test_empty.gtf"
        gtfCLIP.gtfCLIP.fOutput = "UnitTestOutput/testProcess/test_empty.bed"
        gtfCLIP.gtfCLIP.geneType = "gene_type"
             
        gtfC = gtfCLIP.gtfCLIP(options)
        gtfC.processGTF()
        outputTest = open("UnitTestOutput/testProcess/test_empty.bed", "r")
        checkTest  = open("testProcess/checkTestProcess/test_empty_check.bed", "r")
         
        v1 = []
        v2 = []
         
        for line in outputTest:
            v1.append(line)
         
        for line in checkTest:
            v2.append(line)
             
        outputTest.close()
        checkTest.close()
              
        self.assertEqual(v1, v2)
         
####################################################################################
###############################JUNCTION UNIT TESTS##################################
####################################################################################
class TestJunctionMethods(unittest.TestCase):
      
    def test01_junction(self):
         
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testJunction/test01.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testJunction/test01.txt"
        bedCLIP.bedCLIP.fCompare = "testJunction/test01_junction.bed"
          
        bedC = bedCLIP.bedCLIP(options)
        bedC.junction()
        outputTest = open("UnitTestOutput/testJunction/test01.txt", "r")
        checkTest  = open("testJunction/checkTestJunction/test01_check.txt", "r")
          
        v1 = []
        v2 = []
          
        for line in outputTest:
            v1.append(line)
          
        for line in checkTest:
            v2.append(line)
              
        outputTest.close()
        checkTest.close()
               
        self.assertEqual(v1, v2)
          
    def test02_junction(self):
          
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testJunction/test02.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testJunction/test02.txt"
        bedCLIP.bedCLIP.fCompare = "testJunction/test02_junction.bed"
          
        bedC = bedCLIP.bedCLIP(options)
        bedC.junction()
        outputTest = open("UnitTestOutput/testJunction/test02.txt", "r")
        checkTest  = open("testJunction/checkTestJunction/test02_check.txt", "r")
          
        v1 = []
        v2 = []
          
        for line in outputTest:
            v1.append(line)
          
        for line in checkTest:
            v2.append(line)
              
        outputTest.close()
        checkTest.close()
               
        self.assertEqual(v1, v2)
          
    def test03_junction(self):
          
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testJunction/test03.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testJunction/test03.txt"
        bedCLIP.bedCLIP.fCompare = "testJunction/test03_junction.bed"
          
        bedC = bedCLIP.bedCLIP(options)
        bedC.junction()
        outputTest = open("UnitTestOutput/testJunction/test03.txt", "r")
        checkTest  = open("testJunction/checkTestJunction/test03_check.txt", "r")
          
        v1 = []
        v2 = []
          
        for line in outputTest:
            v1.append(line)
          
        for line in checkTest:
            v2.append(line)
              
        outputTest.close()
        checkTest.close()
               
        self.assertEqual(v1, v2)
                
    def test04_junction(self):
          
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testJunction/test04.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testJunction/test04.txt"
        bedCLIP.bedCLIP.fCompare = "testJunction/test04_junction.bed"
          
        bedC = bedCLIP.bedCLIP(options)
        bedC.junction()
        outputTest = open("UnitTestOutput/testJunction/test04.txt", "r")
        checkTest  = open("testJunction/checkTestJunction/test04_check.txt", "r")
          
        v1 = []
        v2 = []
          
        for line in outputTest:
            v1.append(line)
          
        for line in checkTest:
            v2.append(line)
              
        outputTest.close()
        checkTest.close()
               
        self.assertEqual(v1, v2)
          
####################################################################################
################################COUNT UNIT TESTS####################################
####################################################################################
class TestCountMethods(unittest.TestCase):
       
    def test01_count(self):
         
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testCount/test01.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test01.txt"
        bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.testCount.sorted.bed.gz"
        bedCLIP.bedCLIP.choice = "o"
           
        bedC = bedCLIP.bedCLIP(options)
        bedC.count_only()
        outputTest = open("UnitTestOutput/testCount/test01.txt", "r")
        checkTest  = open("testCount/checkTestCount/test01_check.txt", "r")
           
        v1 = []
        v2 = []
           
        for line in outputTest:
            v1.append(line)
           
        for line in checkTest:
            v2.append(line)
               
        outputTest.close()
        checkTest.close()
                
        self.assertEqual(v1, v2)
          
    def test02_count(self):
           
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testCount/test02.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test02.txt"
        bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.testCount.sorted.bed.gz"
        bedCLIP.bedCLIP.choice = "o"
           
        bedC = bedCLIP.bedCLIP(options)
        bedC.count_only()
        outputTest = open("UnitTestOutput/testCount/test02.txt", "r")
        checkTest  = open("testCount/checkTestCount/test02_check.txt", "r")
           
        v1 = []
        v2 = []
           
        for line in outputTest:
            v1.append(line)
           
        for line in checkTest:
            v2.append(line)
               
        outputTest.close()
        checkTest.close()
                
        self.assertEqual(v1, v2)
          
    def test03_count(self):
           
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testCount/test03.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test03.txt"
        bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.testCount.sorted.bed.gz"
        bedCLIP.bedCLIP.choice = "o"
           
        bedC = bedCLIP.bedCLIP(options)
        bedC.count_only()
        outputTest = open("UnitTestOutput/testCount/test03.txt", "r")
        checkTest  = open("testCount/checkTestCount/test03_check.txt", "r")
           
        v1 = []
        v2 = []
           
        for line in outputTest:
            v1.append(line)
           
        for line in checkTest:
            v2.append(line)
               
        outputTest.close()
        checkTest.close()
                
        self.assertEqual(v1, v2)
          
    def test04_count(self):
           
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testCount/test04.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test04.txt"
        bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.testCount.sorted.bed.gz"
        bedCLIP.bedCLIP.choice = "o"
           
        bedC = bedCLIP.bedCLIP(options)
        bedC.count_only()
        outputTest = open("UnitTestOutput/testCount/test04.txt", "r")
        checkTest  = open("testCount/checkTestCount/test04_check.txt", "r")
           
        v1 = []
        v2 = []
           
        for line in outputTest:
            v1.append(line)
           
        for line in checkTest:
            v2.append(line)
               
        outputTest.close()
        checkTest.close()
                
        self.assertEqual(v1, v2)
          
    def test05_count(self):
           
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testCount/test05.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test05.txt"
        bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.testCount.sorted.bed.gz"
        bedCLIP.bedCLIP.choice = "o"
           
        bedC = bedCLIP.bedCLIP(options)
        bedC.count_only()
        outputTest = open("UnitTestOutput/testCount/test05.txt", "r")
        checkTest  = open("testCount/checkTestCount/test05_check.txt", "r")
           
        v1 = []
        v2 = []
           
        for line in outputTest:
            v1.append(line)
           
        for line in checkTest:
            v2.append(line)
               
        outputTest.close()
        checkTest.close()
                
        self.assertEqual(v1, v2)
          
    def test06_count(self):
           
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testCount/test01.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test06.txt"
        bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.sorted.bed.gz"
        bedCLIP.bedCLIP.choice = "a"
           
        bedC = bedCLIP.bedCLIP(options)
        bedC.count_all()
        outputTest = open("UnitTestOutput/testCount/test06.txt", "r")
        checkTest  = open("testCount/checkTestCount/test06_check.txt", "r")
           
        v1 = []
        v2 = []
           
        for line in outputTest:
            v1.append(line)
           
        for line in checkTest:
            v2.append(line)
               
        outputTest.close()
        checkTest.close()
                
        self.assertEqual(v1, v2)
          
####################################################################################
############################SLIDING WINDOW UNIT TESTS###############################
####################################################################################
class TestSlidingWindowMethods(unittest.TestCase):
      
    def test_01(self):
         
        options = {}
         
        gtfCLIP.gtfCLIP.fInput = "testSlidingWindow/test01.bed"
        gtfCLIP.gtfCLIP.fOutput = "UnitTestOutput/testSlidingWindow/test01.bed"
        gtfCLIP.gtfCLIP.windowSize = 100
        gtfCLIP.gtfCLIP.windowStep = 50
          
        gtfC = gtfCLIP.gtfCLIP(options)
        gtfC.slidingWindow()
        outputTest = open("UnitTestOutput/testSlidingWindow/test01.bed", "r")
        checkTest  = open("testSlidingWindow/checkTestSlidingWindow/test01_check.bed", "r")
          
        v1 = []
        v2 = []
          
        for line in outputTest:
            v1.append(line)
          
        for line in checkTest:
            v2.append(line)
              
        outputTest.close()
        checkTest.close()
               
        self.assertEqual(v1, v2)
          
    def test_02(self):
          
        options = {}
         
        gtfCLIP.gtfCLIP.fInput = "testSlidingWindow/test02.bed"
        gtfCLIP.gtfCLIP.fOutput = "UnitTestOutput/testSlidingWindow/test02.bed"
        gtfCLIP.gtfCLIP.windowSize = 100
        gtfCLIP.gtfCLIP.windowStep = 50
          
        gtfC = gtfCLIP.gtfCLIP(options)
        gtfC.slidingWindow()
        outputTest = open("UnitTestOutput/testSlidingWindow/test02.bed", "r")
        checkTest  = open("testSlidingWindow/checkTestSlidingWindow/test02_check.bed", "r")
          
        v1 = []
        v2 = []
          
        for line in outputTest:
            v1.append(line)
          
        for line in checkTest:
            v2.append(line)
              
        outputTest.close()
        checkTest.close()
               
        self.assertEqual(v1, v2)
          
####################################################################################
############################COUNTSW UNIT TESTS######################################
####################################################################################
class TestCountSWMethods(unittest.TestCase):
           
    def test01(self):
         
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testCountSW/test01_cl.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCountSW/test01.txt"
        bedCLIP.bedCLIP.fCompare = "testCountSW/test01.bed"
            
        bedC = bedCLIP.bedCLIP(options)
        bedC.countSlidingWindow()
        outputTest = open("UnitTestOutput/testCountSW/test01.txt", "r")
        checkTest  = open("testCountSW/checkTestCountSW/test01_check.sw.txt", "r")
            
        v1 = []
        v2 = []
            
        for line in outputTest:
            v1.append(line)
            
        for line in checkTest:
            v2.append(line)
                
        outputTest.close()
        checkTest.close()
                 
        self.assertEqual(v1, v2)
           
    def test02(self):
            
        options = {}
          
        bedCLIP.bedCLIP.fInput = "testCountSW/test02_cl.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCountSW/test02.txt"
        bedCLIP.bedCLIP.fCompare = "testCountSW/test01.bed"
             
        bedC = bedCLIP.bedCLIP(options)
        bedC.countSlidingWindow()
        outputTest = open("UnitTestOutput/testCountSW/test02.txt", "r")
        checkTest  = open("testCountSW/checkTestCountSW/test02_check.sw.txt", "r")
             
        v1 = []
        v2 = []
             
        for line in outputTest:
            v1.append(line)
             
        for line in checkTest:
            v2.append(line)
                 
        outputTest.close()
        checkTest.close()
                  
        self.assertEqual(v1, v2)
            
    def test03(self):
            
        options = {}
          
        bedCLIP.bedCLIP.fInput = "testCountSW/test01_cl.bed"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCountSW/test03.txt"
        bedCLIP.bedCLIP.fCompare = "testCountSW/test02.bed"
             
        bedC = bedCLIP.bedCLIP(options)
        bedC.countSlidingWindow()
        outputTest = open("UnitTestOutput/testCountSW/test03.txt", "r")
        checkTest  = open("testCountSW/checkTestCountSW/test03_check.sw.txt", "r")
             
        v1 = []
        v2 = []
             
        for line in outputTest:
            v1.append(line)
             
        for line in checkTest:
            v2.append(line)
                 
        outputTest.close()
        checkTest.close()
                  
        self.assertEqual(v1, v2)
                 
####################################################################################
############################ToDexSeq WINDOW UNIT TESTS##############################
####################################################################################
class TestToDexSeqMethods(unittest.TestCase):       
       
    def test01(self):
         
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testToDexSeq/test01.txt"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testToDexSeq/test01.dex.txt"
            
        bedC = bedCLIP.bedCLIP(options)
        bedC.toDEXSeq()
        outputTest = open("UnitTestOutput/testToDexSeq/test01.dex.txt", "r")
        checkTest  = open("testToDexSeq/checkToDexSeq/test01_check.dex.txt", "r")
            
        v1 = []
        v2 = []
            
        for line in outputTest:
            v1.append(line)
            
        for line in checkTest:
            v2.append(line)
                
        outputTest.close()
        checkTest.close()
                 
        self.assertEqual(v1, v2)
          
    def test02(self):
          
        options = {}
         
        bedCLIP.bedCLIP.fInput = "testToDexSeq/test02.txt"
        bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testToDexSeq/test02.dex.txt"
            
        bedC = bedCLIP.bedCLIP(options)
        bedC.toDEXSeq()
        outputTest = open("UnitTestOutput/testToDexSeq/test02.dex.txt", "r")
        checkTest  = open("testToDexSeq/checkToDexSeq/test02_check.dex.txt", "r")
            
        v1 = []
        v2 = []
            
        for line in outputTest:
            v1.append(line)
            
        for line in checkTest:
            v2.append(line)
                
        outputTest.close()
        checkTest.close()
                 
        self.assertEqual(v1, v2)

if __name__ == '__main__':
    unittest.main()
