# --------------------------------------------------
# GeneLengthSummary class
# Authors: Thomas Schwarzl, schwarzl@embl.de
# ---------------------------------------------------

import gzip
import logging
from collections import OrderedDict, defaultdict

import HTSeq
from .output import Output


"""
The Class 'GeneLengthSummary' manages global length summary for genes. 
The length of the gene regions on each chromosome as well as length per
gene type is stored.
"""

class GeneLengthSummary:
    splitExons = True
    
    def __init__(self):
        self.chromosomes = defaultdict(int)
        self.genetypes = defaultdict(int)

    # adds a length to a chromosome 
    def addToChromosome(self, chrom, length):
        self.chromosomes[ chrom ] += int(length)
    
    # adds a length to a gene type
    def addToGeneType(self, type, length):
        self.genetypes[ type ] += int(length)
        
    # add a gene
    def add(self, gene):
        if gene is not None:
            # region list which is used for adding to the summary
            countList = None
            
            if self.splitExons:
                countList = gene.detailedRegionList
            else:
                countList = gene.regionList
            
            if countList is not None:
                for region in countList:
                    self.addToGeneType(region.type, region.length())
                    self.addToChromosome(region.interval.chrom, region.length())
    
    # print    
    def __str__(self):
        out = ""
        for chrom,length in self.chromosomes.items():
            out += "track chr %s %s\n" % (str(chrom), str(length))
        for type,length in self.genetypes.items():
            out += "track type %s %s\n" % (str(type), str(length))
        return str(out)
