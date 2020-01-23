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
        pass
        
