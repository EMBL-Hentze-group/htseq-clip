# --------------------------------------------------
# GeneRegion class
# Authors: Thomas Schwarzl, schwarzl@embl.de
# --------------------------------------------------

import gzip
import logging
from collections import OrderedDict, defaultdict
import HTSeq

from .output import Output


"""
The Class 'GeneRegion' manages a genomic region of a gene.
"""
class GeneRegion:

    def __init__(self, gene, interval):
        logging.debug("Init GeneRegion")

        # corresponding gene
        self.gene           = gene
        # genomic interval
        self.interval       = interval

        self.type           = None
        self.index          = None
        self.total          = None
        self.upstreamFlag   = None
        self.downstreamFlag = None
        

    """
    Creates a boolean encoded flag
    F - F =  0, both ends cannot be trusted
    T - F =  1, only upstream end can be trusted
    F - T =  2, only downstream end can be trusted
    T - T =  3, both ends can be trusted
    trust = no alternative start sites or ends, respectively, were found,
    therefore all given 
    # TODO: Quick hack, could be done with bit conversion.
    """
    def getFlag(self):
        if self.downstreamFlag is None:
            raise Exception("downstream flag does not exit")
        if self.upstreamFlag is None:
            raise Exception("upstream flag does not exist")

        flag = -1

        if self.upstreamFlag == False and self.downstreamFlag == False:
            flag = 0
        elif self.upstreamFlag == False and self.downstreamFlag == True:
            flag = 2 
        elif self.upstreamFlag == True and self.downstreamFlag == False:
            flag = 1
        else: # both True
            flag = 3 

        logging.debug("Calculated flag: (%s, %s) -> %s" % (self.upstreamFlag, self.downstreamFlag, flag))
        return flag

    def isExon(self):
        return self.type == self.gene.__EXON__

    def isIntron(self):
        return self.type == self.gene.__INTRON__
    
    def length(self):
        return(self.interval.end - self.interval.start)

    """
    returns a .bed format style output string
    """
    def toBed(self):
        '''
        ouput format:
        chr begin   end Id@Name@GeneType@FeatureType@FeatureNumber/TotalFeatures@Id:FeatureType00FeatureNumber  confidence  strand
        '''
        return("\t".join( ( str(self.interval.chrom),
                   str(self.interval.start),
                   str(self.interval.end),
                   "@".join((self.gene.getId(), self.gene.getGeneName(),
                             self.gene.getGeneType(), self.type, "{}/{}".format(self.index,self.total),
                             "{0}:{1}{2:0>4}".format(self.gene.getId(),self.type,self.index))),
                    str(self.getFlag()), self.interval.strand)))

    def __str__(self):
        return self.toBed()

    """
    Splits exon regions into CDS and UTRs
    """
    def split(self):
        # only split if the region is an exon
        if self.isExon():
            return self.doSplit()
        # return if region is an intron
        elif self.isIntron():
            return [ self ]
        else:
            raise Exception("Region is neither Exon or Intron. In any case," + 
                            "biology is wrong, not the program. Please contact the developer.")
    
    """
    Splits an exon region into CDS and UTRs
    """
    def doSplit(self):
        out = []
        
        logging.debug("splitting exon region %s" % self)
        
        for (regionInterval, regionType) in self.gene.details[ self.interval ].steps():
            region = GeneRegion(self.gene, regionInterval)
            region.type = regionType
            region.index = self.index
            region.total = self.total

            #IMPORTANT TODO
            #TODO flag propagation is not done right, neglected for CDS/intron
            region.upstreamFlag = self.upstreamFlag
            region.downstreamFlag = self.downstreamFlag
            out.append( region )
        return(out) 
