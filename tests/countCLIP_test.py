import unittest
from argparse import Namespace
from clip.countCLIP import countCLIP

class TestcountCLIP(unittest.TestCase):
    '''
    unit tests for module countCLIP
    '''
    def test_annotation(self):
        '''
        unit test for subparser command mapToId
        '''
        options = Namespace(annotation='tests/testcountCLIP/W50S20.bed',output='tests/testcountCLIP/checkcountCLIP/test_annotation.txt')
        countC = countCLIP(options)
        countC.annotationToIDs()
        annSet, expSet = set(),set()
        with open('tests/testcountCLIP/checkcountCLIP/test_annotation.txt','r') as ann:
            for l in ann:
                annSet.add(l)
        with open('tests/UnitTestOutput/testcountCLIP/test_annotation.txt','r') as eann:
            for l in eann:
                expSet.add(l)
        self.assertSetEqual(annSet,expSet)
    
    # def test_count(self):
    #     '''
    #     unit test for subparser command count
    #     '''
    #     options = Namespace(annotation='tests/testcountCLIP/W50S20.bed',output='tests/testcountCLIP/checkcountCLIP/test_slidingWindows_count.txt',
    #         input='tests/testcountCLIP/end_sites.bed.gz')
    #     countSW = countCLIP(options)
    #     countSW.count()


if __name__ == '__main__':
    unittest.main()