# --------------------------------------------------
# Gene class
# Authors: Thomas Schwarzl, schwarzl@embl.de
# --------------------------------------------------


import gzip, HTSeq
from collections import OrderedDict
from output import Output
import logging
from collections import defaultdict

#logging.basicConfig(level = logging.DEBUG)
logging.basicConfig(level = logging.ERROR)

from GTxFeature import GTxFeature
from GeneRegion import GeneRegion
from HTSeq import GenomicFeature

"""
The Class 'Gene' stores genomic gene location information
and flattens gene info. 

It sores the Exon information of different transcripts
and calculates the Introns from the inbetween spaces. 

Also, the CDS and UTR information is stored and flatten.
"""
class Gene:
    logger = logging.getLogger(__name__)  
    __CDS__    = "CDS"
    __3UTR__   = "3UTR"
    __5UTR__   = "5UTR"
    __EXON__   = "exon"
    __INTRON__ = "intron"

    """
    'Gene': init
    """
    def __init__(self, feature):
        self.logger.debug("Initializing new gene")

        self.stranded = True
        self.forwardSymbol = "+"
        self.reverseSymbol = "-"
        self.regionPriorityOrder = (self.__CDS__,
                                    self.__3UTR__,
                                    self.__5UTR__)

        self.features = {
            self.__CDS__  : HTSeq.GenomicArrayOfSets('auto', stranded = self.stranded),
            self.__3UTR__ : HTSeq.GenomicArrayOfSets('auto', stranded = self.stranded),
            self.__5UTR__ : HTSeq.GenomicArrayOfSets('auto', stranded = self.stranded),
            self.__EXON__ : HTSeq.GenomicArrayOfSets('auto', stranded = self.stranded)
        }

        self.details = HTSeq.GenomicArrayOfSets('auto', stranded = self.stranded)

        self.startSites = {
            self.__CDS__  : HTSeq.GenomicArray('auto', stranded = self.stranded),
            self.__3UTR__ : HTSeq.GenomicArray('auto', stranded = self.stranded),
            self.__5UTR__ : HTSeq.GenomicArray('auto', stranded = self.stranded),
            self.__EXON__ : HTSeq.GenomicArray('auto', stranded = self.stranded)
        }

        self.endSites = {
            self.__CDS__  : HTSeq.GenomicArray('auto', stranded = self.stranded),
            self.__3UTR__ : HTSeq.GenomicArray('auto', stranded = self.stranded),
            self.__5UTR__ : HTSeq.GenomicArray('auto', stranded = self.stranded),
            self.__EXON__ : HTSeq.GenomicArray('auto', stranded = self.stranded)
        }

        self.storage = { 
            self.__CDS__  : [],
            self.__3UTR__ : [],
            self.__5UTR__ : []
        }

        self.exonTotal   = 0
        self.intronTotal = 0

        # List of unprocessed gene regions
        self.rawRegionList    = None
        
        # List of GeneRegions containing the processed exons and introns
        self.regionList    = list()

        # List of GeneRegions containing the processed CDS, UTRs, exons and introns
        self.detailedRegionList = None

        # gene symbol and gene type
        self.symbol   = "NA"
        self.geneType = "NA"
        
        # Determines if Exons should be split into CDS and UTR regions
        self.splitExons    = True
        
        self.feature = feature
        
        
    """
    'Gene': Getters for convenience
    """
    def getGeneInterval(self):
        return self.feature.getInterval()

    def getGeneType(self):
        return self.feature.getGeneType()

    def getGeneName(self):
        return self.feature.getGeneName()

    def getStrand(self):
        return self.feature.getStrand()

    def getName(self):
        return self.feature.getName()

    """
    'Gene': Basic functions
    """
    def isForwardStrand(self):
        return self.getStrand() == self.forwardSymbol

    def isReverseStrand(self):
        return self.getStrand() == self.reverseSymbol

    """
    Returns true if the object was processed, which is required for returning output
    """
    def isProcessed(self):
        return len(self.regionList) > 0
    
    
    """
    Calculating total exon and total intron number
    """
    def calculateTotalExonAndIntronCount(self):
        self.logger.debug("processing gene")
        
        self.exonTotal  = (len(self.rawRegionList) + 1) / 2
        self.intronTotal = self.exonTotal - 1
        
        self.logger.debug("exon total: %s, intron total: %s" % (str(self.exonTotal), str(self.intronTotal)))
    
    
    """
    Merge/Flatten the exons and store into regionList
    """
    def mergeExons(self):
        self.rawRegionList = list(self.features[ self.__EXON__ ][ self.getGeneInterval() ].steps())
    
    
    def exonsWereAdded(self):
        return len(list(self.features[ self.__EXON__ ].steps())) > 1
    """
    Gene annotation:
    Calculates the exon, intron regions and their corresponding flags
    The processed annotation is stored in the variable self.regionList
    """
    def process(self):
        self.logger.debug("processing gene")
        
        if self.exonsWereAdded():
        
            self.mergeExons()
                    
            # in this function the flattened gene definition is created, also all start and end sites
                # of exons and introns are stored so alternative isoforms of exons and introns can be assigned.
                # those are used to calculate flags, (no alternative isoform, or 5' or 3' isoform, or 5' and 3'
                # isoform variants
            self.calculateExonsAndIntrons()
            self.processStoredRegions()
            
            if self.splitExons:
                self.splitExonsIntoUTRandCDSRegions()
        else:
            if self.processGeneOnlyInformation:
                
                self.logger.debug("adding an exon for a gene without exon information.")
                
                feature = GTxFeature(GenomicFeature("name", self.__EXON__, self.feature.feature.iv))
                
                self.addRegion(feature, self.__EXON__)
                
                self.process()
            else:
                raise Exception("The gene annotation file provides a gene without any exon information. " +
                                "Either add exon information to the annotation or use the processGeneOnlyInformation " +
                                "of htseq-clip.")
                
        
    """
    Gene.getGeneLength:
    returns length of gene (sum of all gene regions).
    """
    def length(self):
        length = 0
        for region in self.regionList:
            length += region.length()
        return length
    
    """
    Gene: Checks if the preprocessing was done, which is an essential step
    before providing output 
    """
    def checkIfProcessed(self):
        if not self.isProcessed():
            self.process()

    """
    Gene: Calculates exon and intron regions and their corresponding flags 
    and number, as well as the total count (e.g. exon 1/10).
    """
    def calculateExonsAndIntrons(self):
        self.calculateTotalExonAndIntronCount()
        
        self.logger.debug("calculating exons and introns and their corresponding flags, number, and total count")
        
        exonIndex   = 1 if self.isForwardStrand() else self.exonTotal 
        intronIndex = 1 if self.isForwardStrand() else self.intronTotal 

        regionIndex = 0
        
        for (regionInterval, regionInfo) in self.rawRegionList:
            
            self.logger.debug("processing region %s - %s" % (str(regionInterval), str(regionInfo)))
            
            region = GeneRegion(self, regionInterval)

            # if exon
            if len(regionInfo) > 0: # == "exon":
                self.logger.debug("processing exon")
                
                upstreamFlag   = self.getExonUpstreamFlag(regionInterval)
                downstreamFlag = self.getExonDownstreamFlag(regionInterval)

                region.type = self.__EXON__
                region.index = exonIndex
                region.total = self.exonTotal
                region.upstreamFlag = upstreamFlag
                region.downstreamFlag = downstreamFlag
                
                exonIndex = self.nextIndex(exonIndex)

            # else intron
            else:
                self.logger.debug("processing intron")
                
                region.type = self.__INTRON__
                region.index = intronIndex
                region.total = self.intronTotal
                # intron flags will be determined after all exon flags have been assigned
                
                intronIndex = self.nextIndex(intronIndex)

            # update regionList
            self.regionList.append(region)

            regionIndex += 1
        
        # calculate all intron flags
        self.calculateIntronFlags()

    """
    Calculates flags for Introns by retrieving flags from the neighboring regions.
    """
    def calculateIntronFlags(self):
        self.logger.debug("calculating intron flags directional")

        regionIndex = 1
        
        while regionIndex < len(self.regionList):
            self.regionList[ regionIndex ].upstreamFlag   = self.regionList[ self.previousIndex(regionIndex) ].downstreamFlag
            self.regionList[ regionIndex ].downstreamFlag = self.regionList[ self.previousIndex(regionIndex) ].upstreamFlag

            regionIndex += 2
    
   
    """
    Splits exon regions into UTR and CDS regions and stores all regions
    to 'detailedRegionList'
    """
    def splitExonsIntoUTRandCDSRegions(self):
        self.logger.debug("calculate UTR and CDS regions")
        
        self.detailedRegionList = []
        
        for region in self.regionList:
            for newRegion in region.split():
                self.detailedRegionList.append(newRegion)

    """
    Gene: These functions increments or decrements the index,
    depending on strandness of Gene
    """
    def nextIndex(self, index):
        return self.indexStep(index, 1)
    
    def previousIndex(self, index):
        return self.indexStep(index, -1)

    def indexStep(self, index, step):
        if self.isForwardStrand():
            index += step
        elif self.isReverseStrand():
            index -= step
        else:
            raise Exception("Sorry, but htseq-clip cannot work with unstranded data yet.")
        return index

    """
    Gene: Get the strand specific exon upstream flag
    """
    def getExonUpstreamFlag(self, interval):
        length = len(list(self.startSites[ self.__EXON__ ][ interval ].steps()))
        self.logger.debug("Get exon start sites in interval: %s" % length)
        return length == 1

    """
    Gene: Get the strand specific exon downstream flag
    """
    def getExonDownstreamFlag(self, interval):
        length = len(list(self.endSites[ self.__EXON__ ][ interval ].steps()))
        self.logger.debug("Get exon end sites in interval: %s" % length)
        return length == 1

    """
    Gene: Adding a gff/gtf feature to the gene.
    """
    def add(self, feature):
        if self.isProcessed():
            raise Exception("Gene already has been processed, you cannot add additional regions.")
        
        self.logger.debug("adding info %s " % feature)

        self.assertFeatureBelongsToGene(feature)
        
        if feature.isExon():
            self.logger.debug("invoking addRegion %s %s" % (feature, self.__EXON__))
            self.addRegion(feature, self.__EXON__)

        elif feature.isCDS():
            self.storeRegion(feature, self.__CDS__)
       
        elif feature.is5UTR():
            self.storeRegion(feature, self.__5UTR__)
        
        elif feature.is3UTR():
            self.storeRegion(feature, self.__3UTR__)
        
        else:
            self.logger.debug("ignoring feature %s" % feature)
    
    """
    Assert that the feature belongs to the gene
    """
    def assertFeatureBelongsToGene(self, feature):
        if not self.feature.getGeneName() == feature.getGeneName():
            raise SyntaxError("The order of gene and gene features in the input file are incorrect. The current feature does not belong to the gene being processed.")

    """
    Gene: Adds a region to genes and adds the corresponding
    start and end sites. Start and end sites are later used to 
    figure out alternative isoforms for a given exon or intron
    
    """
    def addRegion(self, feature, name):
        self.logger.debug("adding {} {} {}".format(name, feature, feature.getInterval()))
        
        self.features[ name ][ feature.getInterval() ] = name 
        self.startSites[ name ][ HTSeq.GenomicPosition( feature.getInterval().chrom,
                                                        feature.getInterval().start_d,
                                                        strand = feature.getStrand() ) ] = True 
        self.endSites[ name ][ HTSeq.GenomicPosition(   feature.getInterval().chrom,
                                                        feature.getInterval().end_d,
                                                        strand = feature.getStrand() ) ] = True 
        self.details[ feature.getInterval() ] = name
        
        self.logger.debug("finished adding %s %s" % (name, feature))


    """
    Gene: stores the region for processing.
    Processing can be only started once all the regions are stored.
    """
    def storeRegion(self, feature, name):
        self.logger.debug("storing %s %s" % (str(name), str(feature)))
        self.storage[ name ].append(feature)

    
    """ 
    This function adds regions according to the region priority order
    """
    def processStoredRegions(self):
        self.logger.debug("processing stored regions")
        
        for regionName in self.regionPriorityOrder:
            for feature in self.storage[ regionName ]:
                self.addRegion(feature, regionName)


    """"
    Gene: Outputs the genomic location info of 
    exons and introns to a bed format 
    """
    def toBed(self):
        self.checkIfProcessed()

        for region in self.regionList:
            yield region.toBed()

    """"
    Gene: Outputs the genomic location info of
    CDS, UTRs, introns and remaining exons to 
    a bed format 
    """
    def toBedDetailed(self):
        # calculates all the coordinates if not calculated
        self.checkIfProcessed()
        
        regionList = None
        
        if self.splitExons:
            regionList = self.detailedRegionList
        else:
            regionList = self.regionList
        
        # write the individual regions to the output
        for region in regionList:
            yield region.toBed()
