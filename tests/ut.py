# --------------------------------------------------
# UNIT TEST CLASS
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
#          Nadia Ashraf, nadia.ashraf@embl.de
# Institution: EMBL Heidelberg
# Date: January 2016
# --------------------------------------------------

import unittest
import sys
sys.path.append("..")
from clip.bamCLIP import bamCLIP
from clip.bedCLIP import bedCLIP
from clip.gtfCLIP import gtfCLIP
from clip.gtf import gtfClip
from clip.gffCLIP import gffClip


####################################################################################
############################EXTRACTION UNIT TESTS###################################
####################################################################################
class TestExtractMethods(unittest.TestCase):
     
    def test01_extract_SS(self):
        
        options = {}
        
        bamC = bamCLIP(options)  
        
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
        
        bamC = bamCLIP(options)  
        
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
        
        bamC = bamCLIP(options)  
        
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
        
        bamC = bamCLIP(options)  
        
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
        
        bamC = bamCLIP(options)  
        
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
         
        bamC = bamCLIP(options)
        
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
         
        bamC = bamCLIP(options)
        
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
         
        bamC = bamCLIP(options)
        
        bamC.fInput = "testExtract/test02.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test02_ES.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        
        bamC.choice = "e"

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
         
        bamC = bamCLIP(options)
        
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
         
        bamC = bamCLIP(options)
        
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
         
        bamC = bamCLIP(options)
        
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
         
        bamC = bamCLIP(options)
        
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
         
        bamC = bamCLIP(options)
        
        bamC.fInput = "testExtract/test03.bam"
        bamC.fOutput = "UnitTestOutput/testExtract/test03_ES.bed"
        bamC.minAlignmentQuality = 10
        bamC.minReadLength = 0
        bamC.maxReadLength = 0
        bamC.primary = False
        bamC.choice = "e"
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
         
        bamC = bamCLIP(options)
        
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
         
        bamC = bamCLIP(options)
        
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

# ####################################################################################
# ###############################PROCESS UNIT TESTS###################################
# ####################################################################################
# class TestProcessMethods(unittest.TestCase):
     
#     def test01_process(self):
        
#         options = {}
        
#         gtfCLIP.gtfFile = "testProcess/test.gtf"
#         gtfCLIP.fOutput = "UnitTestOutput/testProcess/test01.bed"
#         gtfCLIP.geneType = "gene_type"
         
#         gtfC = gtfCLIP(options)
#         gtfC.processGTF()
#         outputTest = open("UnitTestOutput/testProcess/test01.bed", "r")
#         checkTest  = open("testProcess/checkTestProcess/test_check.bed", "r")
         
#         v1 = []
#         v2 = []
         
#         for line in outputTest:
#             v1.append(line)
         
#         for line in checkTest:
#             v2.append(line)
             
#         outputTest.close()
#         checkTest.close()
              
#         self.assertEqual(v1, v2)

#     def test02_process(self):

#         options = {}

#         gtfCLIP.gtfFile = "testProcess/gtftry.gtf"
#         gtfCLIP.fOutput = "UnitTestOutput/testProcess/test02.bed"
#         gtfCLIP.geneType = "gene_type"
#         gtfCLIP.geneName = 'gene_name'

#         gtfC = gtfCLIP(options)
#         gtfC.processGTF()
#         outputTest = open("UnitTestOutput/testProcess/test02.bed", "r")
#         checkTest  = open("testProcess/checkTestProcess/test02_check.bed", "r")

#         v1 = []
#         v2 = []

#         for line in outputTest:
#             v1.append(line)

#         for line in checkTest:
#             v2.append(line)

#         outputTest.close()
#         checkTest.close()

#         self.assertEqual(v1, v2)

#     def test03_process(self):

#         options = {}

#         gtfClip.gtfFile = "testProcess/gtftry.gtf"
#         gtfClip.fOutput = "UnitTestOutput/testProcess/test03.bed"
#         gtfClip.geneType = "gene_type"

#         gtfC = gtfClip(options)
#         gtfC.processGTF()
#         outputTest = open("UnitTestOutput/testProcess/test03.bed", "r")
#         checkTest  = open("testProcess/checkTestProcess/test03_check.bed", "r")

#         v1 = []
#         v2 = []

#         for line in outputTest:
#             v1.append(line)

#         for line in checkTest:
#             v2.append(line)

#         outputTest.close()
#         checkTest.close()

#         self.assertEqual(v1, v2)

#     def test04_process(self):

#         options = {}

#         gtfClip.gtfFile = "testProcess/gtftry.gtf"
#         gtfClip.fOutput = "UnitTestOutput/testProcess/test04.bed"
#         gtfClip.geneType = "gene_type"
#         gtfClip.geneName = 'gene_name'

#         gtfC = gtfClip(options)
#         gtfC.processGTF()
#         outputTest = open("UnitTestOutput/testProcess/test04.bed", "r")
#         checkTest  = open("testProcess/checkTestProcess/test04_check.bed", "r")

#         v1 = []
#         v2 = []

#         for line in outputTest:
#             v1.append(line)

#         for line in checkTest:
#             v2.append(line)

#         outputTest.close()
#         checkTest.close()

#         self.assertEqual(v1, v2)

#     def test05_process(self):

#         options = {}

#         gtfCLIP.gtfFile = "testProcess/test.gff3"
#         gtfCLIP.fOutput = "UnitTestOutput/testProcess/test05.bed"
#         gtfCLIP.geneType = "biotype"

#         gtfC = gtfCLIP(options)
#         gtfC.processGTF()
#         outputTest = open("UnitTestOutput/testProcess/test05.bed", "r")
#         checkTest  = open("testProcess/checkTestProcess/test05_check.bed", "r")

#         v1 = []
#         v2 = []

#         for line in outputTest:
#             v1.append(line)

#         for line in checkTest:
#             v2.append(line)

#         outputTest.close()
#         checkTest.close()

#         self.assertEqual(v1, v2)

#     def test06_process(self):

#         options = {}

#         gtfCLIP.gtfFile = "testProcess/test.gff3"
#         gtfCLIP.fOutput = "UnitTestOutput/testProcess/test06.bed"
#         gtfCLIP.geneType = "biotype"
#         gtfCLIP.geneName = 'Name'

#         gtfC = gtfCLIP(options)
#         gtfC.processGTF()
#         outputTest = open("UnitTestOutput/testProcess/test06.bed", "r")
#         checkTest  = open("testProcess/checkTestProcess/test06_check.bed", "r")

#         v1 = []
#         v2 = []

#         for line in outputTest:
#             v1.append(line)

#         for line in checkTest:
#             v2.append(line)

#         outputTest.close()
#         checkTest.close()

#         self.assertEqual(v1, v2)

#     def test07_process(self):

#         options = {}

#         gffClip.gtfFile = "testProcess/test.gff3"
#         gffClip.fOutput = "UnitTestOutput/testProcess/test07.bed"
#         gffClip.geneType = "biotype"

#         gtfC = gffClip(options)
#         gtfC.processGFF()
#         outputTest = open("UnitTestOutput/testProcess/test07.bed", "r")
#         checkTest  = open("testProcess/checkTestProcess/test07_check.bed", "r")

#         v1 = []
#         v2 = []

#         for line in outputTest:
#             v1.append(line)

#         for line in checkTest:
#             v2.append(line)

#         outputTest.close()
#         checkTest.close()

#         self.assertEqual(v1, v2)

#     def test08_process(self):

#         options = {}

#         gffClip.gtfFile = "testProcess/test.gff3"
#         gffClip.fOutput = "UnitTestOutput/testProcess/test08.bed"
#         gffClip.geneType = "biotype"
#         gffClip.geneName = 'Name'

#         gtfC = gffClip(options)
#         gtfC.processGFF()
#         outputTest = open("UnitTestOutput/testProcess/test08.bed", "r")
#         checkTest  = open("testProcess/checkTestProcess/test08_check.bed", "r")

#         v1 = []
#         v2 = []

#         for line in outputTest:
#             v1.append(line)

#         for line in checkTest:
#             v2.append(line)

#         outputTest.close()
#         checkTest.close()

#         self.assertEqual(v1, v2)

#     def test_process_empty01(self):
        
#         options = {}
        
#         gtfCLIP.gtfFile = "testProcess/test_empty.gtf"
#         gtfCLIP.fOutput = "UnitTestOutput/testProcess/test01_empty.bed"
#         gtfCLIP.geneType = "gene_type"
             
#         gtfC = gtfCLIP(options)
#         gtfC.processGTF()
#         outputTest = open("UnitTestOutput/testProcess/test01_empty.bed", "r")
#         checkTest  = open("testProcess/checkTestProcess/test01_empty_check.bed", "r")
         
#         v1 = []
#         v2 = []
         
#         for line in outputTest:
#             v1.append(line)
         
#         for line in checkTest:
#             v2.append(line)
             
#         outputTest.close()
#         checkTest.close()
              
#         self.assertEqual(v1, v2)

#     def test_process_empty02(self):

#         options = {}

#         gtfCLIP.gtfFile = "testProcess/test_empty.gff3"
#         gtfCLIP.fOutput = "UnitTestOutput/testProcess/test02_empty.bed"
#         gtfCLIP.geneType = "biotype"

#         gtfC = gtfCLIP(options)
#         gtfC.processGTF()
#         outputTest = open("UnitTestOutput/testProcess/test02_empty.bed", "r")
#         checkTest  = open("testProcess/checkTestProcess/test02_empty_check.bed", "r")

#         v1 = []
#         v2 = []

#         for line in outputTest:
#             v1.append(line)

#         for line in checkTest:
#             v2.append(line)

#         outputTest.close()
#         checkTest.close()

#         self.assertEqual(v1, v2)
         
# ####################################################################################
# ###############################JUNCTION UNIT TESTS##################################
# ####################################################################################
# class TestJunctionMethods(unittest.TestCase):
      
#     def test01_junction(self):
         
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testJunction/test01.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testJunction/test01.txt"
#         bedCLIP.bedCLIP.fCompare = "testJunction/test01_junction.bed"
          
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.junction()
#         outputTest = open("UnitTestOutput/testJunction/test01.txt", "r")
#         checkTest  = open("testJunction/checkTestJunction/test01_check.txt", "r")
          
#         v1 = []
#         v2 = []
          
#         for line in outputTest:
#             v1.append(line)
          
#         for line in checkTest:
#             v2.append(line)
              
#         outputTest.close()
#         checkTest.close()
               
#         self.assertEqual(v1, v2)
          
#     def test02_junction(self):
          
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testJunction/test02.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testJunction/test02.txt"
#         bedCLIP.bedCLIP.fCompare = "testJunction/test02_junction.bed"
          
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.junction()
#         outputTest = open("UnitTestOutput/testJunction/test02.txt", "r")
#         checkTest  = open("testJunction/checkTestJunction/test02_check.txt", "r")
          
#         v1 = []
#         v2 = []
          
#         for line in outputTest:
#             v1.append(line)
          
#         for line in checkTest:
#             v2.append(line)
              
#         outputTest.close()
#         checkTest.close()
               
#         self.assertEqual(v1, v2)
          
#     def test03_junction(self):
          
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testJunction/test03.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testJunction/test03.txt"
#         bedCLIP.bedCLIP.fCompare = "testJunction/test03_junction.bed"
          
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.junction()
#         outputTest = open("UnitTestOutput/testJunction/test03.txt", "r")
#         checkTest  = open("testJunction/checkTestJunction/test03_check.txt", "r")
          
#         v1 = []
#         v2 = []
          
#         for line in outputTest:
#             v1.append(line)
          
#         for line in checkTest:
#             v2.append(line)
              
#         outputTest.close()
#         checkTest.close()
               
#         self.assertEqual(v1, v2)
                
#     def test04_junction(self):
          
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testJunction/test04.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testJunction/test04.txt"
#         bedCLIP.bedCLIP.fCompare = "testJunction/test04_junction.bed"
          
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.junction()
#         outputTest = open("UnitTestOutput/testJunction/test04.txt", "r")
#         checkTest  = open("testJunction/checkTestJunction/test04_check.txt", "r")
          
#         v1 = []
#         v2 = []
          
#         for line in outputTest:
#             v1.append(line)
          
#         for line in checkTest:
#             v2.append(line)
              
#         outputTest.close()
#         checkTest.close()
               
#         self.assertEqual(v1, v2)
          
# ####################################################################################
# ################################COUNT UNIT TESTS####################################
# ####################################################################################
# class TestCountMethods(unittest.TestCase):
       
#     def test01_count(self):
         
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testCount/test01.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test01.txt"
#         bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.testCount.sorted.bed.gz"
#         bedCLIP.bedCLIP.choice = "o"
           
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.count_only()
#         outputTest = open("UnitTestOutput/testCount/test01.txt", "r")
#         checkTest  = open("testCount/checkTestCount/test01_check.txt", "r")
           
#         v1 = []
#         v2 = []
           
#         for line in outputTest:
#             v1.append(line)
           
#         for line in checkTest:
#             v2.append(line)
               
#         outputTest.close()
#         checkTest.close()
                
#         self.assertEqual(v1, v2)
          
#     def test02_count(self):
           
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testCount/test02.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test02.txt"
#         bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.testCount.sorted.bed.gz"
#         bedCLIP.bedCLIP.choice = "o"
           
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.count_only()
#         outputTest = open("UnitTestOutput/testCount/test02.txt", "r")
#         checkTest  = open("testCount/checkTestCount/test02_check.txt", "r")
           
#         v1 = []
#         v2 = []
           
#         for line in outputTest:
#             v1.append(line)
           
#         for line in checkTest:
#             v2.append(line)
               
#         outputTest.close()
#         checkTest.close()
                
#         self.assertEqual(v1, v2)
          
#     def test03_count(self):
           
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testCount/test03.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test03.txt"
#         bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.testCount.sorted.bed.gz"
#         bedCLIP.bedCLIP.choice = "o"
           
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.count_only()
#         outputTest = open("UnitTestOutput/testCount/test03.txt", "r")
#         checkTest  = open("testCount/checkTestCount/test03_check.txt", "r")
           
#         v1 = []
#         v2 = []
           
#         for line in outputTest:
#             v1.append(line)
           
#         for line in checkTest:
#             v2.append(line)
               
#         outputTest.close()
#         checkTest.close()
                
#         self.assertEqual(v1, v2)
          
#     def test04_count(self):
           
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testCount/test04.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test04.txt"
#         bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.testCount.sorted.bed.gz"
#         bedCLIP.bedCLIP.choice = "o"
           
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.count_only()
#         outputTest = open("UnitTestOutput/testCount/test04.txt", "r")
#         checkTest  = open("testCount/checkTestCount/test04_check.txt", "r")
           
#         v1 = []
#         v2 = []
           
#         for line in outputTest:
#             v1.append(line)
           
#         for line in checkTest:
#             v2.append(line)
               
#         outputTest.close()
#         checkTest.close()
                
#         self.assertEqual(v1, v2)
          
#     def test05_count(self):
           
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testCount/test05.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test05.txt"
#         bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.testCount.sorted.bed.gz"
#         bedCLIP.bedCLIP.choice = "o"
           
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.count_only()
#         outputTest = open("UnitTestOutput/testCount/test05.txt", "r")
#         checkTest  = open("testCount/checkTestCount/test05_check.txt", "r")
           
#         v1 = []
#         v2 = []
           
#         for line in outputTest:
#             v1.append(line)
           
#         for line in checkTest:
#             v2.append(line)
               
#         outputTest.close()
#         checkTest.close()
                
#         self.assertEqual(v1, v2)
          
#     def test06_count(self):
           
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testCount/test01.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCount/test06.txt"
#         bedCLIP.bedCLIP.fCompare = "testCount/Homo_sapiens.GRCh37.82.sorted.bed.gz"
#         bedCLIP.bedCLIP.choice = "a"
           
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.count_all()
#         outputTest = open("UnitTestOutput/testCount/test06.txt", "r")
#         checkTest  = open("testCount/checkTestCount/test06_check.txt", "r")
           
#         v1 = []
#         v2 = []
           
#         for line in outputTest:
#             v1.append(line)
           
#         for line in checkTest:
#             v2.append(line)
               
#         outputTest.close()
#         checkTest.close()
                
#         self.assertEqual(v1, v2)
          
# ####################################################################################
# ############################SLIDING WINDOW UNIT TESTS###############################
# ####################################################################################
# class TestSlidingWindowMethods(unittest.TestCase):
      
#     def test_01(self):
         
#         options = {}
         
#         gtfCLIP.gtfCLIP.fInput = "testSlidingWindow/test01.bed"
#         gtfCLIP.gtfCLIP.fOutput = "UnitTestOutput/testSlidingWindow/test01.bed"
#         gtfCLIP.gtfCLIP.windowSize = 100
#         gtfCLIP.gtfCLIP.windowStep = 50
          
#         gtfC = gtfCLIP.gtfCLIP(options)
#         gtfC.slidingWindow()
#         outputTest = open("UnitTestOutput/testSlidingWindow/test01.bed", "r")
#         checkTest  = open("testSlidingWindow/checkTestSlidingWindow/test01_check.bed", "r")
          
#         v1 = []
#         v2 = []
          
#         for line in outputTest:
#             v1.append(line)
          
#         for line in checkTest:
#             v2.append(line)
              
#         outputTest.close()
#         checkTest.close()
               
#         self.assertEqual(v1, v2)
          
#     def test_02(self):
          
#         options = {}
         
#         gtfCLIP.gtfCLIP.fInput = "testSlidingWindow/test02.bed"
#         gtfCLIP.gtfCLIP.fOutput = "UnitTestOutput/testSlidingWindow/test02.bed"
#         gtfCLIP.gtfCLIP.windowSize = 100
#         gtfCLIP.gtfCLIP.windowStep = 50
          
#         gtfC = gtfCLIP.gtfCLIP(options)
#         gtfC.slidingWindow()
#         outputTest = open("UnitTestOutput/testSlidingWindow/test02.bed", "r")
#         checkTest  = open("testSlidingWindow/checkTestSlidingWindow/test02_check.bed", "r")
          
#         v1 = []
#         v2 = []
          
#         for line in outputTest:
#             v1.append(line)
          
#         for line in checkTest:
#             v2.append(line)
              
#         outputTest.close()
#         checkTest.close()
               
#         self.assertEqual(v1, v2)
          
# ####################################################################################
# ############################COUNTSW UNIT TESTS######################################
# ####################################################################################
# class TestCountSWMethods(unittest.TestCase):
           
#     def test01(self):
         
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testCountSW/test01_cl.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCountSW/test01.txt"
#         bedCLIP.bedCLIP.fCompare = "testCountSW/test01.bed"
            
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.countSlidingWindow()
#         outputTest = open("UnitTestOutput/testCountSW/test01.txt", "r")
#         checkTest  = open("testCountSW/checkTestCountSW/test01_check.sw.txt", "r")
            
#         v1 = []
#         v2 = []
            
#         for line in outputTest:
#             v1.append(line)
            
#         for line in checkTest:
#             v2.append(line)
                
#         outputTest.close()
#         checkTest.close()
                 
#         self.assertEqual(v1, v2)
           
#     def test02(self):
            
#         options = {}
          
#         bedCLIP.bedCLIP.fInput = "testCountSW/test02_cl.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCountSW/test02.txt"
#         bedCLIP.bedCLIP.fCompare = "testCountSW/test01.bed"
             
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.countSlidingWindow()
#         outputTest = open("UnitTestOutput/testCountSW/test02.txt", "r")
#         checkTest  = open("testCountSW/checkTestCountSW/test02_check.sw.txt", "r")
             
#         v1 = []
#         v2 = []
             
#         for line in outputTest:
#             v1.append(line)
             
#         for line in checkTest:
#             v2.append(line)
                 
#         outputTest.close()
#         checkTest.close()
                  
#         self.assertEqual(v1, v2)
            
#     def test03(self):
            
#         options = {}
          
#         bedCLIP.bedCLIP.fInput = "testCountSW/test01_cl.bed"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testCountSW/test03.txt"
#         bedCLIP.bedCLIP.fCompare = "testCountSW/test02.bed"
             
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.countSlidingWindow()
#         outputTest = open("UnitTestOutput/testCountSW/test03.txt", "r")
#         checkTest  = open("testCountSW/checkTestCountSW/test03_check.sw.txt", "r")
             
#         v1 = []
#         v2 = []
             
#         for line in outputTest:
#             v1.append(line)
             
#         for line in checkTest:
#             v2.append(line)
                 
#         outputTest.close()
#         checkTest.close()
                  
#         self.assertEqual(v1, v2)
                 
# ####################################################################################
# ############################ToDexSeq WINDOW UNIT TESTS##############################
# ####################################################################################
# class TestToDexSeqMethods(unittest.TestCase):       
       
#     def test01(self):
         
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testToDexSeq/test01.txt"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testToDexSeq/test01.dex.txt"
            
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.toDEXSeq()
#         outputTest = open("UnitTestOutput/testToDexSeq/test01.dex.txt", "r")
#         checkTest  = open("testToDexSeq/checkToDexSeq/test01_check.dex.txt", "r")
            
#         v1 = []
#         v2 = []
            
#         for line in outputTest:
#             v1.append(line)
            
#         for line in checkTest:
#             v2.append(line)
                
#         outputTest.close()
#         checkTest.close()
                 
#         self.assertEqual(v1, v2)
          
#     def test02(self):
          
#         options = {}
         
#         bedCLIP.bedCLIP.fInput = "testToDexSeq/test02.txt"
#         bedCLIP.bedCLIP.fOutput = "UnitTestOutput/testToDexSeq/test02.dex.txt"
            
#         bedC = bedCLIP.bedCLIP(options)
#         bedC.toDEXSeq()
#         outputTest = open("UnitTestOutput/testToDexSeq/test02.dex.txt", "r")
#         checkTest  = open("testToDexSeq/checkToDexSeq/test02_check.dex.txt", "r")
            
#         v1 = []
#         v2 = []
            
#         for line in outputTest:
#             v1.append(line)
            
#         for line in checkTest:
#             v2.append(line)
                
#         outputTest.close()
#         checkTest.close()
                 
#         self.assertEqual(v1, v2)

if __name__ == '__main__':
    unittest.main()
