# --------------------------------------------------
# GTxFeature class
# Authors: Thomas Schwarzl, schwarzl@embl.de
# --------------------------------------------------

from builtins import str
from builtins import object

import gzip
import logging
from collections import OrderedDict, defaultdict

import HTSeq


"""
The Class 'GTxFeature'  
"""
class GTxFeature(object):
    # String which is returned if attribute for identifier,
    # e.g. geneTypeAttrib for geneType, or geneNameAttrib for 
    # geneName is not available
    __NA__ = "NA"

    #TODO: definitions should come from an option
    geneDefinitions      = ("tRNA", "gene", "tRNAscan", "tRNA_gene")
    exonDefinitions      = ("exon", "exonic")
# @TODO review this list!!!
# "C_gene_segment",
# "D_gene_segment",
# "D_loop",
# "J_gene_segment",
# "RNase_MRP_RNA",
# "RNase_P_RNA",
# "SRP_RNA",
# "V_gene_segment",
# "Y_RNA",
# "antisense_RNA",
# "cDNA_match",
# "lnc_RNA",
# "mRNA",
# "match",
# "miRNA",
# "ncRNA",
# "primary_transcript",
# "rRNA",
# "sequence_feature",
# "snRNA",
# "snoRNA",
# "tRNA",
# "telomerase_RNA",
# "vault_RNA") # "transcript",
    cdsDefinitions       = ("CDS")
    utr3Definitions      = ("three_prime_UTR", "UTR3", "3UTR")
    utr5Definitions      = ("five_prime_UTR", "UTR5", "5UTR")
    utrDefinitions       = ("UTR") # undefined if 5 or 3 prime
    utrNameSplit         = ":"
    utrNameSplitPosition = 0
    
    # Identifier for the attribute field containing the gene type 
    # in the feature
    geneTypeAttrib = None # geneTypeId
    
    # Identifier for the attribute field containing the gene name  
    # in the feature 
    geneNameAttrib = None # geneNameAttrib
    geneIdAttrib = None
    
    def __init__(self, feature):
        self.feature    = feature
        self.geneType   = None
        self.geneName   = None
        self.geneId = None
    
    def getInterval(self):
        return self.feature.iv

    def getStrand(self):
        return self.feature.iv.strand

    def getName(self):
        return self.feature.name

    """
    Returns the geneType of the feature.
    If the identifier 'self.geneTypeAttrib' is not set, it will return the not available 
    string self.__NA__.
    """
    def getGeneType(self):
        if self.geneType is None:
            self.setGeneType()
         
        return self.geneType

    """
    Returns the geneName of the feature, if the identifier 'self.geneTypeName' is set. 
    """
    def getGeneName(self):
        if self.geneName is None:
            self.setGeneName()
        
        return self.geneName
    
    def getGeneId(self):
        if self.geneId is None:
            self.setGeneId()
        return self.geneId
    
    def setGeneType(self):
        if self.geneTypeAttrib is not None:
            if self.geneTypeAttrib in self.feature.attr:
                self.geneType = str(self.feature.attr[str(self.geneTypeAttrib)])
            else:
                # @TODO reformat this message
                error = "Gene type attribute -t '" + self.geneTypeAttrib + "' was not found. Please use -t to specifiy the identifier used in your annotation file to describe the gene type. Alternatively, locate and change the record with the missing identifier. Feature: " + str(self.feature)
                raise KeyError(error)
        else:
            self.geneType = self.__NA__

    def setGeneName(self):
        if self.geneNameAttrib is not None:
            if self.geneNameAttrib in self.feature.attr:
                self.geneName = str(self.feature.attr[str(self.geneNameAttrib)])
            else:
                # @TODO: reformat this message as well
                error = "Gene name attribute -n '" + self.geneNameAttrib + "' was not found. Please use -n to specify the idenfier used in your annotation file to describe the gene name. Alternatively, locate and change the record with the missing identifier. Feature: " + str(self.feature)
                raise KeyError(error)
        else:
            self.geneName = self.__NA__
    
    def setGeneId(self):
        if self.geneIdAttrib is not None:
            if self.geneIdAttrib in self.feature.attr:
                self.geneId = str(self.feature.attr[str(self.geneIdAttrib)])
            else:
                # @TODO: reformat this message and change the silly attribute letter
                error = "Gene id attribute -a '" + self.geneIdAttrib + "' was not found. Please use -n to specify the idenfier used in your annotation file to describe the gene name. Alternatively, locate and change the record with the missing identifier. Feature: " + str(self.feature)
                raise KeyError(error)
        else:
            self.geneId = self.__NA__

    def isGene(self):
        return self.feature.type in self.geneDefinitions

    def isExon(self):
        return self.feature.type in self.exonDefinitions
    
    def isCDS(self):
        return self.feature.type in self.cdsDefinitions
    
    def is5UTR(self):
        return self.feature.type in self.utr5Definitions or ( self.feature.type == self.utrDefinitions and self.feature.name.split(self.utrNameSplit)[ self.utrNameSplitPosition ] in self.utr5Definitions )

    def is3UTR(self):
        return self.feature.type in self.utr3Definitions or ( self.feature.type == self.utrDefinitions and self.feature.name.split(self.utrNameSplit)[ self.utrNameSplitPosition ] in self.utr3Definitions )
