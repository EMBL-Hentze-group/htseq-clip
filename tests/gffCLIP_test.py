import os
import unittest
from argparse import Namespace

from clip.gffCLIP import gffCLIP

class TestGFFCLIP(unittest.TestCase):
    '''
    Unit tests for module gffCLIP
    '''

    def test_empty(self):
        options = Namespace(gff='tests/testgFFCLIP/test_empty.gff3',output='tests/testgFFCLIP/checkgFFCLIP/test01_empty_check.bed',
            id='gene_id',name='gene_name',type='gene_type',splitExons=True,unsorted=False)
        gffE = gffCLIP(options)
        with self.assertRaises(OSError):
            gffE.process()
    
    def test_comments(self):
        options = Namespace(gff='tests/testgFFCLIP/test_comments.gff3',output='tests/testgFFCLIP/checkgFFCLIP/test01_comments_check.bed',
        id='gene_id',name='gene_name',type='gene_type',splitExons=True,unsorted=False)
        gffC = gffCLIP(options)
        with self.assertRaises(ValueError):
            gffC.process()

if __name__ == '__main__':
    unittest.main()