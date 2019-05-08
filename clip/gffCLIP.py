# --------------------------------------------------
# --------------------------------------------------
# gffCLIP class
# Authors: Thomas Schwarzl, schwarzl@embl.de
#          Nadia Ashraf, nadia.ashraf@embl.de
#          Marko Fritz, marko.fritz@embl.de
# EMBL Heidelberg
# --------------------------------------------------


import gzip
import logging
import sys
from collections import OrderedDict, defaultdict

import HTSeq
from Gene import Gene
from GeneLengthSummary import GeneLengthSummary
from GeneRegion import GeneRegion
from GTxFeature import GTxFeature
from output import Output

"""
The class gffCLIP reads in genomic location,
flattens the annotation and outputs the 
results as Bed file
"""

class GeneInfo(object):

    '''
    class for all gene feature related info, such as gene, intron, exon...
    '''

    def __init__(self, geneDefinitions = ["tRNA", "gene", "tRNAscan", "tRNA_gene"]):
        self._featList = []
        self._geneDefinitions = set(geneDefinitions)
        self._typeMap = {} # feature type to list index map
    
    def add_feature(self,f):
        '''
        add gene feature
        Arguments:
         f: HTSeq.GenomicFeature object
        '''
        self._featList.append(f)
        try:
            self._typeMap[f.type].append(len(self._featList)-1)
        except KeyError:
            self._typeMap[f.type]=[len(self._featList)-1]

    @property
    def gene(self):
        '''
        Get gene object
        @TODO: convert sys.stderr calls to logging.warning
        '''
        geneDef = list(self._geneDefinitions & set(self._typeMap.keys()))
        if len(geneDef)==0:
            sys.stderr.write('No gene definitions found!...Skipping\n')
            return None
        elif len(geneDef)>1:
            sys.stderr.write('Multiple gene definitions found!...Skipping\n')
            return None
        else:
            geneInd = self._typeMap[geneDef[0]]
            if len(geneInd)>1:
                sys.stderr.write('Multiple gene features found!...Skipping\n')
                return None
            else:
                return self._featList[geneInd[0]]
    
    @property
    def features(self):
        '''
        Return all stored features
        '''
        return self._featList



class gffCLIP:
    features = {}
    inputFile  = ''
    gtfFile    = ''
    outputFile = ''
    geneTypeAttrib = ''
    geneNameAttrib = ''
    geneIdAttrib = ''
    windowSize = 50
    windowStep = 20
    detailed   = True
    splitExons = True
    processGeneOnlyInformation = False
    geneDefinitions = ["tRNA", "gene", "tRNAscan", "tRNA_gene"] 
    
    logger    = logging.getLogger(__name__) 

    """
    gffCLIP: Init
    """
    def __init__(self, options):
        # @TODO: add gene ID option
        if hasattr(options, 'gtf'):
            self.gtfFile = options.gtf
            self.logger.debug("set gtf")

        if hasattr(options, 'windowSize'):
            self.windowSize = options.windowSize

        if hasattr(options, 'windowStep'):
            self.windowStep = options.windowStep

        if hasattr(options, 'type'):
            self.geneType = options.type

        if hasattr(options,'name'):
            self.geneName = options.name
        
        if hasattr(options,'id'):
            self.geneId = options.id
        
        if hasattr(options,'geneOnly'):
            self.processGeneOnlyInformation = options.geneOnly
        
        #TODO:
        if hasattr(options, 'detailed'):
            self.detailed = options.detailed
        self.fOutput = Output(options.output)
        self._geneMap = None # for unsorted GFF files
        self.summary = None # GeneSummary

    def process(self,unsorted = False):
        """
        This method goes through the gtf file and determines 
        the positions of exons and introns.

        Attention. At the moment, this is done iteratively. 
        All gene info of a gene has to be provided in 
        consecutive lines. An exception will be raised if
        the annotation is out of order.
        """
        self.logger.debug("Reading from '%s'" % self.gtfFile)
        
        # initializing gff reader
        gtf = HTSeq.GFF_Reader(self.gtfFile, end_included=True)
        
        
        # initializing gene length summary, which records total length of
        # chromosomes, gene types etc
        self.summary = GeneLengthSummary()
        self.summary.splitExons = self.splitExons
        
        # initializing field identifiers
        # @TODO: change ID to Attrib and add geneIdAttrib
        GTxFeature.geneTypeAttrib = self.geneType
        GTxFeature.geneNameAttrib = self.geneName
        GTxFeature.geneIdAttrib = self.geneId
        if unsorted:
            self._process_unsorted(gtf)
        else:
            gene = self._process_sorted(gtf)
        self._write_summary()
        self.fOutput.close()
        # for f in gtf:
        #     if self.geneName not in f.attr:
        #         continue
        #     try:
        #         self._geneMap[f.attr[self.geneName]].add_feature(f)
        #     except KeyError:
        #         self._geneMap[f.attr[self.geneName]] = GeneInfo()
        #         self._geneMap[f.attr[self.geneName]].add_feature(f)
        # # sanity check
        # if len(self._geneMap)==0:
        #     raise ValueError('Cannot parse gene features from {}! Please check this input file'.format(self.gtfFile))
        
        # sys.stderr.write('Found info for{} genes\n'.format(len(self._geneMap)))
        # # for each feature in gtf file
        # # for each feature in geneMap 
        # for _, geneObj in self._geneMap.items():
        #     if geneObj.gene is None:
        #         continue
        #     gene = Gene(GTxFeature(geneObj.gene))
        #     gene.splitExons = self.splitExons
        #     gene.processGeneOnlyInformation = self.processGeneOnlyInformation
        #     for f in geneObj.features:
        #         if f.type in set(self.geneDefinitions):
        #             continue
        #         gene.add(GTxFeature(f))
        #     self.processGene(gene,summary)
        # print(summary)
            
        # for f in gtf:
        #     # initialize a new feature object
        #     feature = GTxFeature(f)
            
        #     # if this feature is a gene (which is specified as certain ids in the gff field)
        #     if feature.isGene():
        #         # if new gene info is found, annotation the gene info
        #         self.processGene(gene, summary)
                
        #         # then create a new Gene object
        #         gene = Gene(feature)
        #         gene.splitExons = self.splitExons
        #         gene.processGeneOnlyInformation = self.processGeneOnlyInformation
                
        #     # else add the gene region info  
        #     else:
        #         # if there was no gene definition provided, raise an exception
        #         if gene is None:
        #             raise Exception("GTF/GFF file provides gene feature info before the actual gene definition.")
        #         # else add the gene info to the genes
        #         else:
        #             gene.add(feature)
        
        # # annotation the last gene
        # self.processGene(gene)
        
        # print(summary)
        # return gene
    
    def _process_sorted(self,gtf):
        '''
        Helper function, process chromosome co-ordinate sorted GFF files
        '''
        gene = None
        for f in gtf:
            # initialize a new feature object
            feature = GTxFeature(f)
            
            # if this feature is a gene (which is specified as certain ids in the gff field)
            if feature.isGene():
                # if new gene info is found, annotation the gene info
                self.processGene(gene)
                
                # then create a new Gene object
                gene = Gene(feature)
                gene.splitExons = self.splitExons
                gene.processGeneOnlyInformation = self.processGeneOnlyInformation
                
            # else add the gene region info  
            else:
                # if there was no gene definition provided, raise an exception
                if gene is None:
                    raise Exception("GTF/GFF file provides gene feature info before the actual gene definition.")
                # else add the gene info to the genes
                else:
                    gene.add(feature)
        
        # annotation the last gene
        self.processGene(gene)
    
    def _process_unsorted(self,gtf):
        '''
        Helper function, process unsorted GFF files
        Arguments:
         gtf: HTSeq.GFFReader object
        '''
        self._geneMap = {}
        for f in gtf:
            if self.geneName not in f.attr:
                continue
            try:
                self._geneMap[f.attr[self.geneName]].add_feature(f)
            except KeyError:
                self._geneMap[f.attr[self.geneName]] = GeneInfo()
                self._geneMap[f.attr[self.geneName]].add_feature(f)
        # sanity check
        if len(self._geneMap)==0:
            raise ValueError('Cannot parse gene features from {}! Please check this input file'.format(self.gtfFile))
        
        # sys.stderr.write('Found info for {} gene(s)\n'.format(len(self._geneMap)))
        # for each feature in geneMap 
        for _, geneObj in self._geneMap.items():
            if geneObj.gene is None:
                continue
            gene = Gene(GTxFeature(geneObj.gene))
            gene.splitExons = self.splitExons
            gene.processGeneOnlyInformation = self.processGeneOnlyInformation
            for f in geneObj.features:
                if f.type in set(self.geneDefinitions):
                    continue
                gene.add(GTxFeature(f))
            self.processGene(gene)

    """
    Processing gene
    """
    def processGene(self, gene):
        self.logger.debug("annotation gene function")
        if gene is not None:
            #gene.toBed()
            # write the gene as bed entries to output
            for be in gene.toBedDetailed():
                self.fOutput.write(be+'\n')
            # store the nucleotide length of the regions in the summary
            self.summary.add(gene)
    
    def _write_summary(self):
        '''
        Helper function: Write chromosome and gene type summary to the end of bed file
        '''
        for chrom,length in self.summary.chromosomes.items():
            self.fOutput.write('track chr {} {}\n'.format(chrom,length))
        for gtype,length in self.summary.genetypes.items():
            self.fOutput.write('track type {} {}\n'.format(gtype,length))


    """
    This functions calculates the sliding window positions
    """
    def slidingWindow(self,inputFile):
        if inputFile.endswith('.gz'):
            almnt_file  =  gzip.open(inputFile,'r')
        else:
            almnt_file = open(inputFile,'r')
        windowidMap = {}

        for line in almnt_file:
            if line.startswith("track"):
                continue
            line = line.split('\n')
            line = line[0].split('\t')

            name = line[3].split('@')

            # this will creaet duplicate windowCount if there are multiple genes annotated to the same chromosomal locations
            # this needs to be reimplemented
            # if currentName == None or currentName != name[0]:
            #     currentName = name[0]
            #     windowCount = 1
            try:
                windowCount = windowidMap[name[0]]
                windowCount +=1
            except KeyError:
                windowidMap[name[0]] = 1
                windowCount = 1
            strand = line[5]

            start = int(line[1])
            end = int(line[2])

            pos1 = start
            pos2 = start + self.windowSize

            #if length shorter than given windowsize then the whole feature is one window
            #else split the current feature up by the given options
            if pos2 >= end:
                seq = (line[0], str(start), str(end), name[0]+"@"+str(windowCount)+"@"+name[2]+"@"+name[3], line[4], strand)
                self.fOutput.write(str('\t').join(seq) + "\n")

                windowCount+=1
            else:
                while pos2 < end:

                    seq = (line[0], str(pos1), str(pos2), name[0]+"@"+str(windowCount)+"@"+name[2]+"@"+name[3], line[4], strand)
                    self.fOutput.write(str('\t').join(seq) + "\n")

                    pos1 = pos1 + self.windowStep
                    pos2 = pos2 + self.windowStep

                    windowCount+=1

                    if pos2 > end:
                        pos2 = end
                        seq = (line[0], str(pos1), str(pos2), name[0]+"@"+str(windowCount)+"@"+name[2]+"@"+name[3], line[4], strand)
                        self.fOutput.write(str('\t').join(seq) + "\n")

                        windowCount+=1
            windowidMap[name[0]] = windowCount
        almnt_file.close()
        self.fOutput.close()
