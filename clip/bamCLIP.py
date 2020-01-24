# --------------------------------------------------
# bamCLIP class
# Authors: Marko Fritz, marko.fritz@embl.de
#          Thomas Schwarzl, schwarzl@embl.de
# Institution: EMBL Heidelberg
# Date: October 2015
# modified by Sudeep Sahadevan, sudeep.sahadevan@embl.de
# --------------------------------------------------

import decimal
import gzip
import logging

from HTSeq import BAM_Reader, GenomicPosition
from .output import Output


class bamCLIP(object):
    
    # default parameters
    data = {}
    fInput = ""
    fOutput = ""
    maxReadLength = 0
    minReadLength = 0
    maxReadIntervalLength = 10000
    minAlignmentQuality = 10
    primary = False
    choice = ''
    mate = 1
    
    count = 0
    
    def __init__(self, options):
        self.fInput = options.input
        self.writeFile = True
        self.fOutput = options.output
        self.choice = options.choice
        self.mate = options.mate           
        self.data = {'maxReadLength' : options.maxReadLength,
                     'minReadLength' : options.minReadLength,
                     'primary': options.primary,
                     'maxReadIntervalLength': options.maxReadIntervalLength,
                     'minAlignmentQuality': options.minAlignmentQuality}
              
    #================================================================================ 
    '''
    This method determines if a read fullfills the critera to be included in the analysis.
    '''
    def readFullfillsQualityCriteria(self, almnt):
        if self.mate == 2: #select the second read of a pair to extract
            if almnt.paired_end and almnt.pe_which == "second":
                return(almnt.aligned and almnt.iv.length <= self.data['maxReadIntervalLength'] and
                       almnt.aQual >= self.data['minAlignmentQuality'] and not almnt.failed_platform_qc
                       and self.primaryFilter(almnt))
            elif almnt.paired_end and almnt.pe_which == "first":
                return False
            elif not almnt.paired_end:
                raise Exception("Alignment is not paired end.")
        elif self.mate == 1: #select the first read of a pair to extract
            if ( not almnt.paired_end ) or (almnt.paired_end and almnt.pe_which =="first"):
                return(almnt.aligned and
                       almnt.iv.length <= self.data['maxReadIntervalLength'] and
                       almnt.aQual >= self.data['minAlignmentQuality'] and
                       not almnt.failed_platform_qc and #SAM flag 0x0200
                       self.primaryFilter(almnt))
            elif almnt.paired_end and almnt.pe_which == "second":
                return False
        else:
            raise ValueError("Mate argument can only be 1 for first read or 2 for second")

    '''
    If the primary filter is activated (option primary is set),
    for multimapping reads it will filter out the best location
    using the SAM flag 0x0100
    '''
    def primaryFilter(self, almnt):
        if self.data['primary']:
            return(not almnt.not_primary_alignment)
        else:
            return(True) 

    '''
    Calculates the read length and stores the min and max length
    '''
    def getSequenceLength(self, almnt):
        return(len(almnt.read.seq))   
       
    def calcMinMax(self, almnt):
        length = self.getSequenceLength(almnt)
      
        if not (self.data['maxReadLength'] == 0 and self.data['minReadLength'] == 0):
            if length > self.data['maxReadLength']:
                self.data['maxReadLength'] = length
            if length < self.data['minReadLength']:
                self.data['minReadLength'] = length
        return(length)

    '''
    Returns a list of Genomic intervals for the positions of deletions/insertions/mutations
    '''
    def parseCigar(self, almnt):
        variations = {'deletions': list(),
                     'insertions': list(),
                     'hits': list() }
      
        for i in almnt.cigar:
            if i.type == 'D':
                variations['deletions'].append(i)
            elif i.type == 'I':
                variations['insertions'].append(i)
            elif i.type == 'M':
                variations['hits'].append(i)
         
        return(variations)

    # def posCalcStartSite(self, pos_d, strand, x):
    #     if strand == "+":
    #         return(pos_d + x)
    #     elif strand == "-":
    #         return(pos_d + x)
    #     else:
    #         raise Exception("Strand not known {}".format(strand))

    def posCalcMiddleSite(self, pos_d, strand, x):
        if strand == "+":
            return(pos_d + x)
        elif strand == "-":
            return(pos_d - x)
        else:
            raise Exception("Strand not known {}".format(strand))

    # def posCalcEndSite(self, pos_d, strand, x):
    #     if strand == "+":
    #         return(pos_d + x)
    #     elif strand == "-":
    #         return(pos_d + x)
    #     else:
    #         raise Exception("Strand not known {}".format(strand))

    '''
    extractOptions
    extracts the options ignore and offset from the choice parameter
    '''
    # def extractOptions(self, option):
    #     # option gets ['', '1i'], for eg
    #     ignore = False
    #     if len(option) == 1:
    #         # print option
    #         print self.choice
    #     if len(option) <= 0:
    #         return(ignore, 0)
    #     else:
    #         if option[1] == '':
    #             offset = 0
    #         else:
    #             if "i" in option[1]:
    #                 ignore = True     
    #                 offset = option[1].split("i")
    #                 offset = int(offset[0])
    #             else:
    #                 offset = int(option[1])
            
    #         return((ignore, offset))
 
    '''
    Returns GenomicPosition for middle site and supports Gapped alignments!
    '''
    def determineMiddleSite(self, almnt):

        pos = decimal.Decimal(decimal.Decimal(len(almnt.read.seq))/2).quantize(decimal.Decimal(1),rounding=decimal.ROUND_HALF_UP)
        cigarList = almnt.cigar

        if almnt.iv.strand == "-":
            cigarList = reversed(cigarList)

        for cig in cigarList:

            if cig.type == "M" or cig.type == "I":
                if cig.size - pos > 0:
                    pos = self.posCalcMiddleSite(cig.ref_iv.start_d, cig.ref_iv.strand, pos)
                    break
                else:
                    pos = pos - cig.size
        return(GenomicPosition(almnt.iv.chrom, pos, almnt.iv.strand))

    def getMiddleSiteAsBed(self, almnt):
        pos = self.determineMiddleSite(almnt)
        x = pos.start_d
        y = x + 1

        try:
            yb = almnt.optional_field('YB')
        except Exception:
            yb = 1
            pass
        
        seq = (almnt.iv.chrom, str(min(x,y)), str(max(x,y)), almnt.read.name + "|" + str(len(almnt.read.seq)), str(yb), almnt.iv.strand)
        
        return("\t".join(seq))
 
    '''
    Returns GenomicPosition for end site
    '''
    def determineEndSite(self, iv):
        return(GenomicPosition(iv.chrom, iv.end_d, iv.strand))

    '''
    Writes deletions of an alignment parsed by CIGAR string  as bed line
    '''
    def writeDeletionSiteAsBed(self, almnt, fOutput):
           
        seq = ()
        deletion = []
        
        #Goes through cigar string and only searches for deletions
        #if no deletion is found the read has none
        for i in almnt.cigar:
            if i.type == 'D':
                deletion.append(i)
        
        #if the read has a deletion it will be written out in the bed file
        #for deletions query_from and query_to are always the same
        if deletion:
            
            for i in range(len(deletion)):
                x = deletion[i].ref_iv.start_d
                y = deletion[i].ref_iv.end_d
                if almnt.iv.strand == '-':
                    b = x
                    x = y
                    y = b
                seq = (almnt.iv.chrom, str(x), str(y), almnt.read.name + "|" + str(len(almnt.read.seq)), str(deletion[i].query_from), almnt.iv.strand)
                seq = ("\t".join(seq)) 
                fOutput.write(seq + "\n")    
   
    '''
    Writes insertions of an alignment parsed by CIGAR string as bed line
    '''
    def writeInsertionSiteAsBed(self, almnt, fOutput):
        seq = ()
        insertion = []
        
        #Goes through cigar string and only searches for insertions
        #if no insertion is found the read has none
        for i in almnt.cigar:
            if i.type == 'I':
                insertion.append(i)
        
        #if the read has an insertion it will be written out in the bed file
        if insertion:
            
            for i in range(len(insertion)):
                x = insertion[i].ref_iv.start_d
                y = insertion[i].ref_iv.end_d
                if almnt.iv.strand == '-':
                    b = x
                    x = y
                    y = b
                seq = (almnt.iv.chrom, str(x), str(y), almnt.read.name+"|"+str(len(almnt.read.seq)), str(insertion[i].query_from), almnt.iv.strand)
                seq = ("\t".join(seq)) 
                fOutput.write(seq + "\n")     

    '''
    Returns GenomicPosition for desired site with offset
    '''
    def getOffsetPosition(self, almnt, position, offset = 0, ignore = True):
        # @TODO: trace error and rewrite
        if almnt.iv.strand == "+":
            x = position + offset
        elif almnt.iv.strand == "-":
            x = position - offset
        else:
            raise ValueError("Strand not known {}".format(almnt.iv.strand))

        if x < 0:
            msg = 'Start position cannot be less than zero. Alignment:{} , Read: {}'.format(str(almnt.iv),almnt.read.name)
            if ignore:
                logging.warning('Skipping {}'.format(almnt.read.name))
                logging.warning(msg)
            else:
                raise ValueError(msg)
            return None
        else:
        
            y = x+1
            
            try:
                yb = almnt.optional_field('YB')
            except Exception:
                yb = 1
                pass
            
            seq = (almnt.iv.chrom, str(min(x,y)), str(max(x,y)), almnt.read.name + "|" + str(len(almnt.read.seq)), str(yb), almnt.iv.strand)
    
            return("\t".join(seq))
	            
    '''
    Returns start site (with offset) as bed line 
    '''
    def getStartSiteAsBed(self, almnt, offset, ignore):
        return(self.getOffsetPosition(almnt, almnt.iv.start_d, offset, ignore))

    '''
    Returns end site (with offset) as bed line 
    '''
    def getEndSiteAsBed(self, almnt, offset, ignore):
        return(self.getOffsetPosition(almnt, almnt.iv.end_d, offset, ignore))

    '''
    Extract start sites
    '''
    def extract_StartSites(self, offset = 0, ignore = False):
        almnt_file = BAM_Reader(self.fInput)
        fOutput = Output(self.fOutput)
        for almnt in almnt_file:
            if not self.readFullfillsQualityCriteria(almnt):
                continue
            out = self.getStartSiteAsBed(almnt = almnt, ignore = ignore, offset = offset)
            if out is None:
                continue
            fOutput.write(out + "\n")

        fOutput.close()

    '''
    Extract middle sites
    '''   
    def extract_MiddleSites(self):
        almnt_file = BAM_Reader(self.fInput)
        fOutput = Output(self.fOutput)
    
        for almnt in almnt_file:
            if not self.readFullfillsQualityCriteria(almnt):
                continue
            out = self.getMiddleSiteAsBed(almnt)
            if out is None:
                continue
            fOutput.write(out + "\n")
       
        fOutput.close() 
           
    '''
    Extract end sites. 
    '''   
    def extract_EndSites(self, offset = 0, ignore = False):
        almnt_file = BAM_Reader(self.fInput)
        fOutput = Output(self.fOutput)
        for almnt in almnt_file:
            if not self.readFullfillsQualityCriteria(almnt):
                continue
            out = self.getEndSiteAsBed(almnt = almnt, ignore = ignore, offset = offset)
            if out is None:
                continue
            fOutput.write(out + "\n")

        fOutput.close()
        
    '''
    Extract deletion sites.
    Deletion sites are determined by parsing the CIGAR string.
    ''' 
    def extract_DeletionSites(self):
        almnt_file = BAM_Reader(self.fInput)
        fOutput = Output(self.fOutput)
        for almnt in almnt_file:
            if not self.readFullfillsQualityCriteria(almnt):
                continue
            self.writeDeletionSiteAsBed(almnt, fOutput)
       
        fOutput.close()
        
    '''
    Extract insertion sites.
    Insertion sites are determined by parsing the CIGAR string.
    ''' 
    def extract_InsertionSites(self):
        almnt_file = BAM_Reader(self.fInput)
        fOutput = Output(self.fOutput)
        for almnt in almnt_file:
            if not self.readFullfillsQualityCriteria(almnt):
                continue
            self.writeInsertionSiteAsBed(almnt, fOutput)
       
        fOutput.close()      
