
# --------------------------------------------------
# gffCLIP class
# Authors: Thomas Schwarzl, schwarzl@embl.de
#          Nadia Ashraf, nadia.ashraf@embl.de
#          Marko Fritz, marko.fritz@embl.de
# Modified by: Sudeep Sahadevan, sudeep.sahadevan@embl.de
# EMBL Heidelberg
# --------------------------------------------------

from builtins import str
from builtins import object

import gzip
import logging
import os
import sys
from collections import OrderedDict, defaultdict

import HTSeq
from .Gene import Gene
from .GeneLengthSummary import GeneLengthSummary
from .GeneRegion import GeneRegion
from .GTxFeature import GTxFeature
from .output import Output

"""
The class gffCLIP reads in genomic location,
flattens the annotation and outputs the 
results as Bed file
"""

class EmptyFileException(FileExistsError):
    pass

class NoFeaturesException(ValueError):
    pass

class FeatureOrderException(LookupError):
    pass

class GeneInfo:

    '''
    class for all gene feature related info, such as gene, intron, exon...
    '''

    def __init__(self, id, geneDefinitions = ["tRNA", "gene", "tRNAscan", "tRNA_gene"]):
        self.id  = id
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
            logging.warning('{}: No gene definitions found!...Skipping'.format(self.id))
            return None
        else:
            if(len(geneDef)>1):
                logging.warning('{}: Multiple gene definitions found! using one at random'.format(self.id))
            geneInd = self._typeMap[geneDef[0]]
            if len(geneInd)>1:
                if self._checkGeneFeature(geneInd):
                    logging.warning('{}: Found multiple attributes for the same gene co-ordinates, using one at random'.format(self.id))
                    return self._featList[geneInd[0]]
                else:
                    logging.warning('{}: Multiple gene features found!...Skipping'.format(self.id))
                    return None
            else:
                return self._featList[geneInd[0]]
    
    def _checkGeneFeature(self,geneInd):
        '''
        Helper function, if multiple gene features are found, 
        return true all the co-oridnates are same
        '''
        start_pos,end_pos,strand = set(),set(),set()
        for gi in geneInd:
            start_pos.add(self._featList[gi].iv.start)
            end_pos.add(self._featList[gi].iv.end)
            strand.add(self._featList[gi].iv.strand)
        if len(start_pos)==1 and len(end_pos)==1 and len(strand)==1:
            return True
        else:
            return False
    
    @property
    def features(self):
        '''
        Return all stored features
        '''
        return self._featList

class gffCLIP(object):
    # @TODO: probably remove these
    # features = {}
    # inputFile  = ''
    # gtfFile    = ''
    # outputFile = ''
    # geneTypeAttrib = ''
    # geneNameAttrib = ''
    # geneIdAttrib = ''
    windowSize = 50
    windowStep = 20
    detailed   = True
    splitExons = True
    processGeneOnlyInformation = True
    geneDefinitions = ["tRNA", "gene", "tRNAscan", "tRNA_gene"] 
    

    """
    gffCLIP: Init
    """
    def __init__(self, options):
        if hasattr(options, 'gff'):
            self.gtfFile = options.gff

        if hasattr(options, 'windowSize'):
            self.windowSize = options.windowSize

        if hasattr(options, 'windowStep'):
            self.windowStep = options.windowStep

        if hasattr(options, 'type'):
            self.geneType = options.type

        if hasattr(options,'name'):
            self.geneName = options.name

        if hasattr(options,'splitExons'):
            self.splitExons =  options.splitExons
        
        if hasattr(options,'id'):
            self.geneId = options.id
        
        #@TODO: do these options complicate things
        # if hasattr(options,'geneOnly'):
        #     self.processGeneOnlyInformation = options.geneOnly
        # if hasattr(options, 'detailed'):
        #     self.detailed = options.detailed
        self.fOutput = Output(options.output)
        self._geneMap = None # for unsorted GFF files
        self.summary = None # GeneSummary
        self._decoder = None # decoder for gz files

    def process(self,unsorted = False):
        """
        This method goes through the gtf file and determines 
        the positions of exons and introns.

        Attention. At the moment, this is done iteratively. 
        All gene info of a gene has to be provided in 
        consecutive lines. An exception will be raised if
        the annotation is out of order.
        """
        logging.debug("Reading from {}".format(self.gtfFile))
        # file size check
        if os.path.getsize(self.gtfFile)==0:
            raise EmptyFileException('Input file {} is empty!'.format(self.gtfFile))
        # initializing gff reader
        gtf = HTSeq.GFF_Reader(self.gtfFile, end_included=True)
        # initializing gene length summary, which records total length of
        # chromosomes, gene types etc
        self.summary = GeneLengthSummary()
        self.summary.splitExons = self.splitExons
        
        # initializing field identifiers
        GTxFeature.geneTypeAttrib = self.geneType
        GTxFeature.geneNameAttrib = self.geneName
        GTxFeature.geneIdAttrib = self.geneId
        if unsorted:
            self._process_unsorted(gtf)
        else:
            self._process_sorted(gtf)
        self._write_summary()
        self.fOutput.close()
    
    def _process_sorted(self,gtf):
        '''
        Helper function, process chromosome co-ordinate sorted GFF files
        '''
        gene = None
        featCount = 0
        noGeneFeature = False
        for f in gtf:
            # initialize a new feature object
            feature = GTxFeature(f)
            # if this feature is a gene (which is specified as certain ids in the gff field)
            if feature.isGene():
                # if new gene info is found, annotation the gene info
                self.processGene(gene)
                
                # then create a new Gene object
                gene = Gene(feature,self.splitExons,self.processGeneOnlyInformation)
                featCount+=1
            # else add the gene region info  
            else:
                # if there was no gene definition provided, raise an exception
                if gene is None:
                    noGeneFeature = True
                    break
                # else add the gene info to the genes
                else:
                    gene.add(feature)
                    featCount+=1
        if noGeneFeature:
            raise FeatureOrderException("GTF/GFF file should provide gene feature info before the actual gene definition.")
        if featCount==0:
            raise NoFeaturesException('Cannot parse features from file {}'.format(self.gtfFile))
        
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
            if self.geneId not in f.attr:
                continue
            try:
                self._geneMap[f.attr[self.geneId]].add_feature(f)
            except KeyError:
                self._geneMap[f.attr[self.geneId]] = GeneInfo(f.attr[self.geneId])
                self._geneMap[f.attr[self.geneId]].add_feature(f)
        # sanity check
        if len(self._geneMap)==0:
            raise ValueError('Cannot parse gene features from {}! Please check this input file'.format(self.gtfFile))
        # for each feature in geneMap 
        for _, geneObj in self._geneMap.items():
            _gene = geneObj.gene
            if _gene is None:
                continue
            gene = Gene(GTxFeature(_gene),self.splitExons,self.processGeneOnlyInformation)
            for f in geneObj.features:
                if f.type in set(self.geneDefinitions):
                    continue
                gene.add(GTxFeature(f))
            self.processGene(gene)

    """
    Processing gene
    """
    def processGene(self, gene):
        logging.debug("annotation gene function")
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

    def _toStr(self,line):
        '''
        helper function
        given a string return it as it is
        '''
        return line
    
    def _byteToStr(self,line):
        '''
        helper function
        given bytes decode to string
        '''
        return line.decode('utf-8')

    def _fo(self,fileName):
        '''
        helper function
        return file handle based on file name suffix
        '''
        if fileName.lower().endswith(('.gz','gzip')):
            self._decoder = self._byteToStr
            return gzip.open(fileName,'r')
        else:
            self._decoder = self._toStr
            return open(fileName,'r')
    
    """
    This functions calculates the sliding window positions
    """
    def slidingWindow(self,inputFile):
        windowidMap = {}
        with self._fo(inputFile) as almnt_file:
            for line in almnt_file:
                line = self._decoder(line)
                if line.startswith("track"):
                    continue
                line = line.strip().split('\t')
                name = line[3].split('@')
                featureNumber, _ = name[4].split('/')
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
                # @TODO: add gene name, window id, UID to last
                
                UID = "{0}:{1}{2:0>4}W{3:0>5}".format(name[0],name[3],featureNumber,windowCount)
                while pos2 < end:
                    seq = (line[0], str(pos1), str(pos2), "@".join(name[:5]+[UID,str(windowCount)]), line[4], strand)
                    self.fOutput.write("\t".join(seq) + "\n")
                    pos1 +=  self.windowStep
                    pos2 +=  self.windowStep
                    windowCount+=1
                flen = end - start +1
                if (pos2 >= end) and (flen < self.windowSize) :
                    seq = (line[0], str(start), str(end), "@".join(name[:5]+[UID,str(windowCount)]), line[4], strand)
                    self.fOutput.write("\t".join(seq) + "\n")
                    windowCount+=1
                elif (pos2 >= end):
                    seq = (line[0], str(pos1), str(end), "@".join(name[:5]+[UID,str(windowCount)]), line[4], strand)
                    self.fOutput.write("\t".join(seq) + "\n")
                    windowCount+=1
                windowidMap[name[0]] = windowCount
        self.fOutput.close()
